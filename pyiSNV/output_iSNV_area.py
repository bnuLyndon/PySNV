# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import time
import numpy as np
from pyiSNV.utils import decode


class VariantLocator:
    def __init__(self, config):
        self.kmer_length = config.kmer_length
        self.snv_limit = config.snv_limit
        self.indel_limit = config.indel_limit

    def output_iSNV_area(
        self,
        high_frequency_kmer_dict,
        connection_mat,
        unique_long_kmers,
        long_kmer_counts,
        reads_abun,
        ref_kmer_builder,
    ):

        hf_kmer_set = set()

        ref_kmer_f_dict = ref_kmer_builder.get_ref_kmer("forward")
        seq = ref_kmer_builder.get_genome_seq()

        kmer_length = self.kmer_length
        snv_limit = self.snv_limit
        snv_length_limt = self.indel_limit

        t0 = time.time()

        # snp_list=get_snp_kmer(high_frequency_kmer_dict, ref_kmer_f_dict, reads_abun, kmer_length, snv_limit)

        connection_dict = dict()
        for i in range(len(connection_mat)):
            if connection_mat[i, 0] not in connection_dict.keys():
                connection_dict[connection_mat[i, 0]] = []
            connection_dict[connection_mat[i, 0]].append(connection_mat[i, 1])
        print("connection dict built", time.time() - t0)

        t0 = time.time()
        long_kmer_dict = dict()
        for i in range(len(unique_long_kmers)):
            long_kmer = unique_long_kmers[i, :]
            key = (long_kmer[0], long_kmer[1])
            count = long_kmer_counts[i]
            if key not in long_kmer_dict.keys():
                long_kmer_dict[key] = 0
            long_kmer_dict[key] += count
        # long_kmer_dict = set(unique_long_kmers[:,0]+unique_long_kmers[:,1])
        print("long kmer set built", time.time() - t0)

        # find starting kmer
        iSNV_list = []
        for kmer in connection_dict.keys():
            connection = connection_dict[kmer]
            if kmer != -1 and ref_kmer_f_dict.get(kmer):
                # print(kmers[0], kmers)
                loc = ref_kmer_f_dict.get(kmer)
                # if loc<50:
                #    print(loc, kmers)
                for next_kmer in connection:
                    count = high_frequency_kmer_dict.get(next_kmer, -1)
                    if count > reads_abun[loc] * snv_limit * 0.9:
                        iSNV_list.append([kmer, next_kmer])

        print("number of anchor:", len(iSNV_list))

        # find starting area
        iSNV_starting_list = []
        output_kmer_set = set()
        k = 0
        for temp in iSNV_list:
            start_area_kmer = [temp[0], temp[1]]
            output_kmer_set.update(set([temp[0], temp[1]]))
            k += 2
            start_idx = ref_kmer_f_dict.get(temp[0])

            stack = [start_area_kmer]

            threshold = reads_abun[start_idx] * snv_limit * 0.75

            kk = 0
            while len(stack):
                temp_area_kmer = stack.pop()
                last_kmer = temp_area_kmer[-1]
                kk += 1

                end_idx = ref_kmer_f_dict.get(last_kmer)
                connections = connection_dict.get(last_kmer)

                if len(temp_area_kmer) == kmer_length + 1:
                    check_seq = (temp_area_kmer[0], temp_area_kmer[kmer_length])
                    # check_seq = temp_area_kmer[0]+temp_area_kmer[kmer_length]

                    # if check_seq not in long_kmer_dict.keys():
                    if long_kmer_dict.get(check_seq, 0) < 5:
                        # print('ignoring a psueduo starting at', start_idx)
                        continue
                    else:
                        iSNV_starting_list.append(temp_area_kmer)

                elif end_idx:

                    if start_idx <= end_idx and temp_area_kmer[0] != temp_area_kmer[-1]:
                        iSNV_starting_list.append(temp_area_kmer)

                elif connections:

                    for next_kmer in connections:
                        # print(high_frequency_kmer_dict.get(kmers[1],-1), threshold)
                        if (
                            next_kmer != last_kmer
                            and high_frequency_kmer_dict.get(next_kmer, -1) >= threshold
                            or ref_kmer_f_dict.get(next_kmer)
                        ):

                            new_branch = temp_area_kmer.copy()
                            # print('last kmer: ',last_kmer)
                            new_branch.append(next_kmer)
                            stack.append(new_branch)

        # filtering adapters
        adapter_dict = dict()
        for item in iSNV_starting_list:
            # key=item[-1][:-5]
            key = item[-1] >> 2 * 5
            if key not in adapter_dict.keys():
                adapter_dict[key] = 0
            adapter_dict[key] += 1
        adapter_suffixs = [
            key for key in adapter_dict.keys() if adapter_dict[key] >= 10
        ]
        adapter_suffixs = set(adapter_suffixs)

        f_kmers = np.array(list(connection_dict.keys()))
        f_kmers_suffix = f_kmers >> 2 * 5

        for suffix in adapter_suffixs:
            kmers = f_kmers[f_kmers_suffix == suffix]
            for kmer in kmers:
                high_frequency_kmer_dict[kmer] = -1
                # inverted_table[kmer] = []

        filtered_iSNV_list = []
        for item in iSNV_starting_list:
            if adapter_dict[item[-1] >> 2 * 5] < 10:
                filtered_iSNV_list.append(item)
        print("number of filtered anchor:", len(filtered_iSNV_list))

        mark = False
        # connect kmers
        iSNV_area_list = []
        nb_iSNV = 0
        change_mark = -1
        for temp in filtered_iSNV_list:
            if nb_iSNV % 50 == 0 and nb_iSNV != change_mark:
                change_mark = nb_iSNV
                print("variant region found:", nb_iSNV)
                # print(nb_iSNV, temp[-1])

            start_idx = ref_kmer_f_dict.get(temp[0])

            stack = [temp]

            median_count = np.median(
                [high_frequency_kmer_dict.get(item, -1) for item in temp[1:-2]]
            )
            threshold = max(reads_abun[start_idx], median_count) * snv_limit * 0.75
            # threshold = reads_abun[start_idx]*snv_limit*0.75
            threshold = max(threshold, 1.1)
            long_kmer_threshold = max(threshold / 2, 5)

            if mark:
                break

            kk = 0
            while len(stack):

                temp_area_kmer = stack.pop()
                hf_kmer_set.update(temp_area_kmer)
                if kk > snv_length_limt and len(temp_area_kmer) > snv_length_limt // 2:
                    continue
                if (
                    kk > snv_length_limt * 10
                    and len(temp_area_kmer) > snv_length_limt // 10
                ):
                    continue
                last_kmer = temp_area_kmer[-1]
                kk += 1

                end_idx = ref_kmer_f_dict.get(last_kmer)
                connections = connection_dict.get(last_kmer)

                if len(temp_area_kmer) > kmer_length:

                    check_seq = (
                        temp_area_kmer[len(temp_area_kmer) - kmer_length - 1],
                        temp_area_kmer[len(temp_area_kmer) - 1],
                    )
                    # check_seq = temp_area_kmer[len(temp_area_kmer)-kmer_length-1]+temp_area_kmer[len(temp_area_kmer)-1]
                    # if check_seq not in long_kmer_dict.keys():
                    if long_kmer_dict.get(check_seq, 0) < long_kmer_threshold:
                        # print('ignoring a FP connection at', start_idx)
                        continue

                if end_idx:

                    if end_idx <= start_idx:  # or end_idx - start_idx>300 :
                        continue

                    # var_seq=''
                    for kmer in temp_area_kmer:
                        if kmer == temp_area_kmer[0]:
                            var_seq = decode(kmer, kmer_length, return_base=False)
                        else:
                            var_seq = "{}{}".format(
                                var_seq,
                                decode(kmer, kmer_length, return_base=False)[-1],
                            )

                    if len(var_seq) >= 2 * kmer_length:
                        length_diff = (len(var_seq) - 2 * kmer_length) // 2

                        check_seq_head = (
                            temp_area_kmer[0],
                            temp_area_kmer[kmer_length],
                        )

                        check_seq_mid = (
                            temp_area_kmer[length_diff],
                            temp_area_kmer[length_diff + kmer_length],
                        )
                        check_seq_tail = (
                            temp_area_kmer[len(var_seq) - 2 * kmer_length],
                            temp_area_kmer[len(var_seq) - kmer_length],
                        )

                        # return temp_area_kmer, last_kmer

                        # check_seq_head = temp_area_kmer[0]+temp_area_kmer[kmer_length]
                        # check_seq_mid = temp_area_kmer[length_diff]+temp_area_kmer[length_diff+kmer_length]
                        # check_seq_tail = temp_area_kmer[len(var_seq)-2*kmer_length] + \
                        #    temp_area_kmer[len(var_seq)-kmer_length]
                        # if check_seq_head not in long_kmer_dict or\
                        #    check_seq_mid not in long_kmer_dict or\
                        #    check_seq_tail not in long_kmer_dict:
                        #    continue
                        if (
                            long_kmer_dict.get(check_seq_head, 0) < long_kmer_threshold
                            or long_kmer_dict.get(check_seq_mid, 0)
                            < long_kmer_threshold
                            or long_kmer_dict.get(check_seq_tail, 0)
                            < long_kmer_threshold
                        ):
                            print("ignoring a psueduo connection at", start_idx)
                            continue

                    if start_idx <= end_idx and temp_area_kmer[0] != temp_area_kmer[-1]:

                        area_ref = seq[start_idx : end_idx + kmer_length]
                        area_var_kmers = []
                        area_var_kmer_counts = []
                        for index, kmer in enumerate(temp_area_kmer):
                            area_var_kmers.append(decode(kmer, kmer_length))
                            area_var_kmer_counts.append(
                                high_frequency_kmer_dict.get(kmer, -1)
                            )

                        iSNV_area_list.append(
                            [
                                start_idx,
                                area_ref,
                                area_var_kmers,
                                temp_area_kmer,
                                area_var_kmer_counts,
                            ]
                        )
                        nb_iSNV += 1
                        # print(start_idx,threshold,nb_iSNV,len(stack))

                elif connections:

                    for next_kmer in connections:

                        if next_kmer == last_kmer:  # or kmers[1] in temp_area_kmer:
                            # print(len(stack))
                            drop_list = temp_area_kmer[1 - kmer_length // 2 :]
                            # inverted_table[kmers[0]].remove(kmers[0])
                            for item in stack:
                                if item[-1] in drop_list:
                                    stack.remove(item)
                            continue

                        # print(high_frequency_kmer_dict.get(kmers[1],-1), threshold)
                        if high_frequency_kmer_dict.get(
                            next_kmer, -1
                        ) >= threshold or ref_kmer_f_dict.get(next_kmer):
                            new_branch = temp_area_kmer.copy()
                            # print('last kmer: ',last_kmer)
                            new_branch.append(next_kmer)
                            # connection_dict[connection][count_ind]=-1
                            stack.append(new_branch)

        print(len(iSNV_area_list))

        k = 0
        for kmer in hf_kmer_set:
            k += high_frequency_kmer_dict.get(kmer, 0)

        self.hf_kmer_count = k

        hf_kmer_set = set()
        k = 0
        for area in iSNV_area_list:
            hf_kmer_set.update(area[3])
        for kmer in hf_kmer_set:
            k += high_frequency_kmer_dict.get(kmer, 0)
        self.connected_kmer_count = k
        """
        filterd_list = []
        for region in iSNV_area_list:
            count = seq.count(region[1])
            seq.count(region[1]) > snv_limit * reads_abun[region[0]]
            region.append(count)
            filterd_list.append(region)
        """
        return iSNV_area_list


def output_iSNV_area_bk(
    seq,
    high_frequency_kmer_dict,
    ref_kmer_dict,
    ref_kmer_f_dict,
    connection_mat,
    unique_long_kmers,
    long_kmer_counts,
    reads_abun,
    kmer_length,
    snv_limit=0.04,
    snv_length_limt=300,
):
    t0 = time.time()

    # snp_list=get_snp_kmer(high_frequency_kmer_dict, ref_kmer_f_dict, reads_abun, kmer_length, snv_limit)

    connection_dict = dict()
    for i in range(len(connection_mat)):
        if connection_mat[i, 0] not in connection_dict.keys():
            connection_dict[connection_mat[i, 0]] = []
        connection_dict[connection_mat[i, 0]].append(connection_mat[i, 1])
    print("connection dict built", time.time() - t0)

    t0 = time.time()
    long_kmer_dict = dict()
    for i in range(len(unique_long_kmers)):
        long_kmer = unique_long_kmers[i, :]
        key = (long_kmer[0], long_kmer[1])
        count = long_kmer_counts[i]
        if key not in long_kmer_dict.keys():
            long_kmer_dict[key] = 0
        long_kmer_dict[key] += count
    # long_kmer_dict = set(unique_long_kmers[:,0]+unique_long_kmers[:,1])
    print("long kmer set built", time.time() - t0)

    # iSNV_list = [[2397565332595, 794168308175]]

    # find starting kmer
    iSNV_list = []
    for kmer in connection_dict.keys():
        connection = connection_dict[kmer]
        if kmer != -1 and ref_kmer_f_dict.get(kmer):
            # print(kmers[0], kmers)
            loc = ref_kmer_f_dict.get(kmer)
            # if loc<50:
            #    print(loc, kmers)
            for next_kmer in connection:
                count = high_frequency_kmer_dict.get(next_kmer, -1)
                if count > reads_abun[loc] * snv_limit * 0.9:
                    iSNV_list.append([kmer, next_kmer])

    print("number of anchor:", len(iSNV_list))

    # find starting area
    iSNV_starting_list = []
    output_kmer_set = set()
    k = 0
    for temp in iSNV_list:
        start_area_kmer = [temp[0], temp[1]]
        output_kmer_set.update(set([temp[0], temp[1]]))
        k += 2
        start_idx = ref_kmer_f_dict.get(temp[0])

        stack = [start_area_kmer]

        threshold = reads_abun[start_idx] * snv_limit * 0.75

        kk = 0
        while len(stack):
            temp_area_kmer = stack.pop()
            last_kmer = temp_area_kmer[-1]
            kk += 1

            end_idx = ref_kmer_f_dict.get(last_kmer)
            connections = connection_dict.get(last_kmer)

            if len(temp_area_kmer) == kmer_length + 1:
                check_seq = (temp_area_kmer[0], temp_area_kmer[kmer_length])
                # check_seq = temp_area_kmer[0]+temp_area_kmer[kmer_length]

                # if check_seq not in long_kmer_dict.keys():
                if long_kmer_dict.get(check_seq, 0) < 5:
                    # print('ignoring a psueduo starting at', start_idx)
                    continue
                else:
                    iSNV_starting_list.append(temp_area_kmer)

            elif end_idx:

                if start_idx <= end_idx and temp_area_kmer[0] != temp_area_kmer[-1]:
                    iSNV_starting_list.append(temp_area_kmer)

            elif connections:

                for next_kmer in connections:
                    # print(high_frequency_kmer_dict.get(kmers[1],-1), threshold)
                    if (
                        next_kmer != last_kmer
                        and high_frequency_kmer_dict.get(next_kmer, -1) >= threshold
                        or ref_kmer_f_dict.get(next_kmer)
                    ):

                        new_branch = temp_area_kmer.copy()
                        # print('last kmer: ',last_kmer)
                        new_branch.append(next_kmer)
                        stack.append(new_branch)

    # filtering adapters
    adapter_dict = dict()
    for item in iSNV_starting_list:
        # key=item[-1][:-5]
        key = item[-1] >> 2 * 5
        if key not in adapter_dict.keys():
            adapter_dict[key] = 0
        adapter_dict[key] += 1
    adapter_suffixs = [key for key in adapter_dict.keys() if adapter_dict[key] >= 10]
    adapter_suffixs = set(adapter_suffixs)

    f_kmers = np.array(list(connection_dict.keys()))
    f_kmers_suffix = f_kmers >> 2 * 5

    for suffix in adapter_suffixs:
        kmers = f_kmers[f_kmers_suffix == suffix]
        for kmer in kmers:
            high_frequency_kmer_dict[kmer] = -1
            # inverted_table[kmer] = []

    filtered_iSNV_list = []
    for item in iSNV_starting_list:
        if adapter_dict[item[-1] >> 2 * 5] < 10:
            filtered_iSNV_list.append(item)
    print("number of filtered anchor:", len(filtered_iSNV_list))

    mark = False
    # connect kmers
    iSNV_area_list = []
    nb_iSNV = 0
    change_mark = -1
    for temp in filtered_iSNV_list:
        if nb_iSNV % 50 == 0 and nb_iSNV != change_mark:
            change_mark = nb_iSNV
            print("iSNV found:", nb_iSNV)
            # print(nb_iSNV, temp[-1])

        start_idx = ref_kmer_f_dict.get(temp[0])

        stack = [temp]

        median_count = np.median(
            [high_frequency_kmer_dict.get(item, -1) for item in temp[1:-2]]
        )
        threshold = max(reads_abun[start_idx], median_count) * snv_limit * 0.75
        # threshold = reads_abun[start_idx]*snv_limit*0.75
        threshold = max(threshold, 1.1)
        long_kmer_threshold = max(threshold / 2, 5)

        if mark:
            break

        kk = 0
        while len(stack):

            temp_area_kmer = stack.pop()
            if kk > snv_length_limt and len(temp_area_kmer) > snv_length_limt // 2:
                continue
            if (
                kk > snv_length_limt * 10
                and len(temp_area_kmer) > snv_length_limt // 10
            ):
                continue
            last_kmer = temp_area_kmer[-1]
            kk += 1

            end_idx = ref_kmer_f_dict.get(last_kmer)
            connections = connection_dict.get(last_kmer)

            if len(temp_area_kmer) > kmer_length:

                check_seq = (
                    temp_area_kmer[len(temp_area_kmer) - kmer_length - 1],
                    temp_area_kmer[len(temp_area_kmer) - 1],
                )
                # check_seq = temp_area_kmer[len(temp_area_kmer)-kmer_length-1]+temp_area_kmer[len(temp_area_kmer)-1]
                # if check_seq not in long_kmer_dict.keys():
                if long_kmer_dict.get(check_seq, 0) < long_kmer_threshold:
                    # print('ignoring a FP connection at', start_idx)
                    continue

            if end_idx:

                if end_idx <= start_idx:  # or end_idx - start_idx>300 :
                    continue

                # var_seq=''
                for kmer in temp_area_kmer:
                    if kmer == temp_area_kmer[0]:
                        var_seq = decode(kmer, kmer_length, return_base=False)
                    else:
                        var_seq = "{}{}".format(
                            var_seq, decode(kmer, kmer_length, return_base=False)[-1]
                        )

                if len(var_seq) >= 2 * kmer_length:
                    length_diff = (len(var_seq) - 2 * kmer_length) // 2

                    check_seq_head = (temp_area_kmer[0], temp_area_kmer[kmer_length])

                    check_seq_mid = (
                        temp_area_kmer[length_diff],
                        temp_area_kmer[length_diff + kmer_length],
                    )
                    check_seq_tail = (
                        temp_area_kmer[len(var_seq) - 2 * kmer_length],
                        temp_area_kmer[len(var_seq) - kmer_length],
                    )

                    # return temp_area_kmer, last_kmer

                    # check_seq_head = temp_area_kmer[0]+temp_area_kmer[kmer_length]
                    # check_seq_mid = temp_area_kmer[length_diff]+temp_area_kmer[length_diff+kmer_length]
                    # check_seq_tail = temp_area_kmer[len(var_seq)-2*kmer_length] + \
                    #    temp_area_kmer[len(var_seq)-kmer_length]
                    # if check_seq_head not in long_kmer_dict or\
                    #    check_seq_mid not in long_kmer_dict or\
                    #    check_seq_tail not in long_kmer_dict:
                    #    continue
                    if (
                        long_kmer_dict.get(check_seq_head, 0) < long_kmer_threshold
                        or long_kmer_dict.get(check_seq_mid, 0) < long_kmer_threshold
                        or long_kmer_dict.get(check_seq_tail, 0) < long_kmer_threshold
                    ):
                        print("ignoring a psueduo connection at", start_idx)
                        continue

                if start_idx <= end_idx and temp_area_kmer[0] != temp_area_kmer[-1]:

                    area_ref = seq[start_idx : end_idx + kmer_length]
                    area_var_kmers = []
                    area_var_kmer_counts = []
                    for index, kmer in enumerate(temp_area_kmer):
                        area_var_kmers.append(decode(kmer, kmer_length))
                        area_var_kmer_counts.append(
                            high_frequency_kmer_dict.get(kmer, -1)
                        )

                    iSNV_area_list.append(
                        [
                            start_idx,
                            area_ref,
                            area_var_kmers,
                            temp_area_kmer,
                            area_var_kmer_counts,
                        ]
                    )
                    nb_iSNV += 1
                    # print(start_idx,threshold,nb_iSNV,len(stack))

            elif connections:

                for next_kmer in connections:

                    if next_kmer == last_kmer:  # or kmers[1] in temp_area_kmer:
                        # print(len(stack))
                        drop_list = temp_area_kmer[1 - kmer_length // 2 :]
                        # inverted_table[kmers[0]].remove(kmers[0])
                        for item in stack:
                            if item[-1] in drop_list:
                                stack.remove(item)
                        continue

                    # print(high_frequency_kmer_dict.get(kmers[1],-1), threshold)
                    if high_frequency_kmer_dict.get(
                        next_kmer, -1
                    ) >= threshold or ref_kmer_f_dict.get(next_kmer):
                        new_branch = temp_area_kmer.copy()
                        # print('last kmer: ',last_kmer)
                        new_branch.append(next_kmer)
                        # connection_dict[connection][count_ind]=-1
                        stack.append(new_branch)

    print(len(iSNV_area_list))

    """
    filterd_list = []
    for region in iSNV_area_list:
        count = seq.count(region[1])
        seq.count(region[1]) > snv_limit * reads_abun[region[0]]
        region.append(count)
        filterd_list.append(region)
    """
    return iSNV_area_list
