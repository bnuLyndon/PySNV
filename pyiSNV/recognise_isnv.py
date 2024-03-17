# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:58:54 2022

@author: lilab
"""

import numpy as np
from scipy.stats import binom
from pyiSNV.utils import decode, encode, TRAN_TABLE


class VariantCaller:
    def __init__(self, config):

        self.count_threshold = config.count_limit
        self.error_rate = config.error_rate
        self.vcf_threshold = config.p_threshold
        self.correct_factor = config.correct_factor
        self.indel_limit = config.indel_limit

        self.kmer_dim = config.kmer_length
        self.snp_limit = config.snv_limit
        self.downsample = config.downsample
        self.average_reads_length = config.average_reads_length

        self.verbose = config.verbose

    def output_isnv(
        self,
        iSNV_area_list,
        high_frequency_kmer_dict,
        reads_abun,
        output_file_name,
        ref_kmer_builder,
        sample_mapper,
    ):

        kmer_dim = self.kmer_dim
        snp_limit = self.snp_limit
        downsample = self.downsample

        count_threshold = self.count_threshold
        error_rate = self.error_rate
        vcf_threshold = self.vcf_threshold
        indel_length_limt = self.indel_limit
        correct_factor = self.correct_factor
        average_reads_length = self.average_reads_length

        verbose = self.verbose

        seq = ref_kmer_builder.get_genome_seq()
        snp_list = sample_mapper.snp_list
        border_list = sample_mapper.border_snv_list

        iSNV_dict = {}

        iSNV_area_list = sort_by_count(iSNV_area_list, -1)

        kernel = 4 ** np.array(range(kmer_dim), dtype=np.int64)

        iSNV_dict = recognise_isnv(
            iSNV_area_list, iSNV_dict, kmer_dim, indel_length_limt
        )

        border_SNV_dict = recognise_simple_snv(border_list, iSNV_dict, kmer_dim)

        simple_SNV_dict = recognise_simple_snv(snp_list, iSNV_dict, kmer_dim)

        for key in border_SNV_dict.keys():
            if not iSNV_dict.get(key):
                iSNV_dict[key] = border_SNV_dict[key]

        for key in simple_SNV_dict.keys():
            if not iSNV_dict.get(key):
                iSNV_dict[key] = simple_SNV_dict[key]

        index_dict = {}
        pos_snv_dict = {}
        for key in iSNV_dict.keys():
            if not key:
                continue
            loc = int(iSNV_dict[key].loc) - 1
            ref_length = len(iSNV_dict[key].ref_base)

            if iSNV_dict[key].type == "snp":
                pos1 = loc
                pos2 = loc + 1
            elif (
                iSNV_dict[key].type == "deletion"
                and len(iSNV_dict[key].ref_base) <= average_reads_length
            ):
                pos1 = loc
                pos2 = loc + ref_length
            else:
                pos1 = loc
                pos2 = loc + 1

            for j in range(pos1, pos2):

                if not pos_snv_dict.get(j):
                    pos_snv_dict[j] = set()
                pos_snv_dict[j].add(key)

            if not index_dict.get(loc):
                index_dict[loc] = []
            index_dict[loc].append(key)

            # if loc == 23054:
            #    print(1)
            count = iSNV_dict[key].get_count(high_frequency_kmer_dict)

        snv_depth = np.zeros(reads_abun.shape)
        for pos in pos_snv_dict.keys():
            iSNV_keys = pos_snv_dict.get(pos)
            for key in iSNV_keys:
                # if iSNV_dict[key].type == 'snp' or iSNV_dict[key].type == 'deletion':
                snv_depth[pos] += iSNV_dict[key].count

        base_depth = max_filter(snv_depth + reads_abun, kmer_dim)

        # filter isnv cases
        filtered_isnv_dict = dict()
        filtered_pos_snv_dict = dict()
        key_list = list(index_dict.keys())
        key_list.sort()
        k = 0
        for key in key_list:

            for item in index_dict[key]:
                elems = item.split()
                loc = int(elems[0]) - 1

                if iSNV_dict[item].count <= 0:
                    continue

                if base_depth[loc] * snp_limit <= count_threshold:
                    if (
                        binom.cdf(iSNV_dict[item].count, base_depth[loc], error_rate)
                        < vcf_threshold
                    ):
                        # if binom.pmf(iSNV_dict[item].count, base_depth[loc], 0.01)>0.001:
                        if verbose:
                            print("ignored low cdf", loc)
                        continue
                else:

                    if iSNV_dict[item].count < count_threshold:
                        if verbose:
                            print("ignored low count isnv at", loc)
                        continue

                if elems[1] == "deletion":
                    freq1 = iSNV_dict[item].counts[0] / base_depth[loc]
                    freq2 = (
                        iSNV_dict[item].counts[1] / base_depth[loc + len(elems[2]) - 1]
                    )
                    freq_limit = snp_limit * 0.5
                else:
                    freq1 = iSNV_dict[item].counts[0] / base_depth[loc]
                    freq2 = iSNV_dict[item].counts[1] / base_depth[loc]
                    freq_limit = snp_limit * 0.9

                if freq1 < freq_limit or freq2 < freq_limit:

                    if verbose and freq1 > 0.009 and freq2 > 0.009:
                        print("ignored low freq isnv at", loc, freq1, freq2)
                    continue

                if elems[1] == "deletion":
                    freq = iSNV_dict[item].count / max(
                        base_depth[loc], base_depth[loc + len(elems[2]) - 1]
                    )
                    # freq_limit = snp_limit*0.5
                else:
                    # freq_limit = snp_limit*0.9
                    freq = iSNV_dict[item].count / base_depth[loc]

                # print(freq)
                if freq < freq_limit:

                    if verbose and freq > 0.009:
                        print("ignored low freq isnv at", loc, freq)
                    continue

                filtered_isnv_dict[item] = iSNV_dict[item].count
                if not filtered_pos_snv_dict.get(loc):
                    filtered_pos_snv_dict[loc] = set()
                filtered_pos_snv_dict[loc].add(item)

        # compensite reads depth
        kmer_depth = reads_abun.copy()
        for key in filtered_isnv_dict.keys():

            elems = key.split()
            loc = int(elems[0]) - 1
            ref_length = len(elems[2])

            neighbor_snv = []
            for i in range(loc + 1, loc + kmer_dim):
                if filtered_pos_snv_dict.get(i) and i != loc:
                    keys = filtered_pos_snv_dict[i]
                    for pos_key in keys:
                        if pos_key.split()[1] == "snp":
                            neighbor_snv.append(pos_key)

            if len(neighbor_snv) > 0:
                ref_kmer = seq[loc : loc + kmer_dim]
                neighbor_combo_list = powerSetsBinary(neighbor_snv)
                for neighbor_combo in neighbor_combo_list:
                    if neighbor_combo == []:
                        continue
                    var_kmer = list(ref_kmer)
                    # test_kmer = list(ref_kmer)
                    for neighbor_key in neighbor_combo:
                        var_loc, _, _, var_base = neighbor_key.split()
                        loc_diff = int(var_loc) - loc - 1
                        var_kmer[loc_diff] = var_base
                        # test_kmer[loc_diff] = var_base.lower()

                    var_kmer_code = encode("".join(var_kmer), kernel, TRAN_TABLE)
                    count = high_frequency_kmer_dict.get(var_kmer_code, 0)

                    if count > kmer_depth[loc] * 2:
                        kmer_depth[loc] += count

        updated_base_depth = max_filter(snv_depth + kmer_depth, kmer_dim)

        # avoid false compensation
        for i in range(len(base_depth)):
            if base_depth[i] > updated_base_depth[i] * correct_factor:
                base_depth[i] = updated_base_depth[i]

        # calculate freq for simple snp
        snp_depth = reads_abun.copy()
        simple_snp_set = set()
        for key in filtered_isnv_dict.keys():

            elems = key.split()
            loc = int(elems[0]) - 1
            if elems[1] != "snp" or len(filtered_pos_snv_dict[loc]) != 1:
                continue

            mark = False
            for key in pos_snv_dict.get(loc):
                if filtered_isnv_dict.get(key):
                    mark = True
                    break
            if mark:
                continue

            shift_size = kmer_dim // 2

            neighbor_snv = []
            type_set = set()
            for i in range(loc + 1, loc + kmer_dim):
                if filtered_pos_snv_dict.get(i) and i != loc:
                    keys = filtered_pos_snv_dict[i]
                    for pos_key in keys:
                        neighbor_snv.append(pos_key)
                        type_set.add(pos_key.split()[1])

            if "insertion" in type_set or "deletion" in type_set:
                continue

            simple_snp_set.add(key)
            if len(neighbor_snv) == 0:
                continue

            ref_kmer = seq[loc : loc + kmer_dim]
            neighbor_combo_list = powerSetsBinary(neighbor_snv)
            for neighbor_combo in neighbor_combo_list:
                if neighbor_combo == []:
                    continue
                var_kmer = list(ref_kmer)
                test_kmer = list(ref_kmer)
                for neighbor_key in neighbor_combo:
                    var_loc, _, _, var_base = neighbor_key.split()
                    loc_diff = int(var_loc) - loc - 1
                    var_kmer[loc_diff] = var_base
                    test_kmer[loc_diff] = var_base.lower()

                var_kmer_code = encode("".join(var_kmer), kernel, TRAN_TABLE)
                count = high_frequency_kmer_dict.get(var_kmer_code, 0)

                snp_depth[loc] += count

        # ouput
        k = 0
        # t0 = time.time()
        rows = []
        output_dict = dict()
        new_dict = dict()

        for key in filtered_isnv_dict.keys():

            elems = key.split()
            loc = int(elems[0]) - 1
            depth = int(base_depth[loc])
            count = iSNV_dict[key].count

            if count <= 0:
                continue

            if elems[1] == "deletion" and len(elems[2]) > kmer_dim:
                depth = min(base_depth[loc], base_depth[loc + len(elems[2]) - 1])

            elif key in simple_snp_set:
                depth = snp_depth[loc] + count

            for another_key in filtered_isnv_dict.keys():
                diff = int(another_key.split()[0]) - loc + 1
                if 0 < diff < kmer_dim:
                    depth = max(base_depth[loc], base_depth[loc + diff])

            freq = count / depth
            freq = min(freq, 1)
            freq = np.round(freq, 2)
            if freq < snp_limit:
                continue

            if False:
                string_count = 0
                reads_string = ""
                for seq in iSNV_dict[key].hyplotypes:
                    string_count = max(reads_string.count(seq), string_count)

                if string_count < snp_limit * depth / 4:
                    if verbose:
                        print("fp variant at", loc, string_count, snp_limit * depth / 4)
                    continue

            row = [
                elems[0],
                elems[2],
                elems[3],
                str(freq),
                str(depth * downsample),
            ]  # , str(int(local_abun))]
            output_dict[key] = freq
            new_dict[key] = iSNV_dict[key]

            rows.append(row)
            k += 1

        self.index_dict = index_dict
        self.candidate_dict = iSNV_dict

        self.variant_dict = new_dict
        self.output_dict = output_dict

        self.connected_kmer_counts = cal_connected_kmer_counts(
            high_frequency_kmer_dict, new_dict, kmer_dim, kernel, TRAN_TABLE
        )

        self.base_depth = base_depth
        self.snp_depth = snp_depth

        write_csv(output_file_name, rows)

        return k


class ReadKmer:
    def __init__(self, bases, loc):
        self._kmer = 0
        self._loc = loc
        self._bases = bases

    def __hash__(self):
        return hash(self._bases)

    def __eq__(self, obj):
        if isinstance(obj, ReadKmer):
            if self._bases == obj._bases:
                return True
            else:
                return False

    def set_kmer_code(self, kmer_code):
        self._kmer = kmer_code

    def get_kmer_code(self):
        return self._kmer

    def get_kmer(self):
        return self._bases

    def get_position(self):
        return self._loc


class SNVArea:
    def __init__(
        self, index, loc, ref_seq, var_seq, var_kmers, kmer_dim, pre_dim, sur_dim
    ):
        self.loc = int(loc)
        self.ref_seq = ref_seq
        self.var_seq = var_seq
        self.kmer_dim = kmer_dim

        self.count = 0
        self.index = index

        if not hasattr(var_kmers, "__len__"):
            var_seq = decode(var_kmers, kmer_dim)
            self.var_kmers = [var_kmers for _ in range(len(var_seq))]
        else:
            var_seq = decode(var_kmers[0], kmer_dim)[: kmer_dim - 1]
            for kmer in var_kmers:
                var_seq += decode(kmer, kmer_dim)[-1]
            self.var_kmers = list(var_kmers)

        self.var_seq = var_seq

        self.key = "{}-{}".format(self.ref_seq, self.var_seq)

        self.ref_length = len(ref_seq)
        self.var_length = len(var_seq)

        if pre_dim == 0 and sur_dim == 0:
            pre_dim, sur_dim = get_border_dim(self.ref_seq, self.var_seq)

        self.pre_dim = pre_dim
        self.sur_dim = sur_dim

        aligned_length = min(self.ref_length, self.var_length)
        length_diff = abs(self.var_length - self.ref_length)
        aligned_ref_seq = ["" for _ in range(aligned_length)]
        aligned_var_seq = ["" for _ in range(aligned_length)]

        for i in range(aligned_length):
            if ref_seq[i] == var_seq[i]:
                aligned_ref_seq[i] = ref_seq[i]
                aligned_var_seq[i] = var_seq[i]
            else:
                break

        for j in range(aligned_length):
            if ref_seq[-j] == var_seq[-j]:
                aligned_ref_seq[-j] = ref_seq[-j]
                aligned_var_seq[-j] = var_seq[-j]
            else:
                break

        if self.ref_length > self.var_length:
            aligned_ref_seq[-j] = ref_seq[-length_diff - j : 1 - j]
        elif self.ref_length < self.var_length:
            aligned_var_seq[-j] = var_seq[-length_diff - j : 1 - j]
        else:
            aligned_var_seq[-j] = var_seq[-length_diff - j : 1 - j]
            aligned_ref_seq[-j] = ref_seq[-length_diff - j : 1 - j]

        if "" in aligned_ref_seq or "" in aligned_var_seq:
            aligned_var_seq[i : 1 - j] = [var_seq[i : 1 - j]]
            aligned_ref_seq[i : 1 - j] = [ref_seq[i : 1 - j]]

        if "" in aligned_ref_seq or "" in aligned_var_seq:
            print("warning: might cause depth error")

        # self.match_mask = [True if ref_base == var_base else False for ref_base in aligned_ref_seq \
        #                for var_base in aligned_var_seq]

        self.aligned_ref_seq = aligned_ref_seq
        self.aligned_var_seq = aligned_var_seq

        if min(self.ref_length, self.var_length) - (pre_dim + sur_dim) > 1:
            self.is_simple_area = False
        else:
            self.is_simple_area = True

    def __eq__(self, obj):
        if isinstance(obj, SNVArea):
            if self.key == obj.key:
                return True
            else:
                return False

    def __hash__(self):
        return hash(self.key)

    def _get_ref_seq(self):
        var_kmers = self.var_kmers
        kmer_dim = self.kmer_dim
        if self.var_seq != "":
            return self.var_seq

        var_seq = decode(var_kmers[0], kmer_dim)[: kmer_dim - 1]
        for kmer in var_kmers:
            var_seq += decode(kmer, kmer_dim)[-1]
        self.var_seq = var_seq
        self.var_length = len(var_seq)

        if min(self.ref_length, self.var_length) - (self.pre_dim + self.sur_dim) > 1:
            self.is_simple_area = False
        else:
            self.is_simple_area = True

        return self.var_seq

    def split_var_seq(self, threshold):

        ref_seq = self.aligned_ref_seq
        var_seq = self.aligned_var_seq

        i = -1
        nb_changes = 0
        while True:
            i += 1
            if i >= len(ref_seq):
                break
            if min(len(ref_seq[i]), len(var_seq[i])) <= 1 or len(ref_seq[i]) == len(
                var_seq[i]
            ):
                continue

            lcs_seq, nb_lcstr = longestCommonSubstr(ref_seq[i], var_seq[i])
            ref_lcs_loc = ref_seq[i].find(lcs_seq)

            if (
                self.kmer_dim > len(lcs_seq) >= threshold
                and ref_seq[i][ref_lcs_loc + 1 :].find(lcs_seq) == -1
            ):

                var_lcs_loc = var_seq[i].find(lcs_seq)

                common_seq = list(lcs_seq)

                ref_seq_1 = ref_seq[i][:ref_lcs_loc]
                ref_seq_2 = ref_seq[i][ref_lcs_loc + len(lcs_seq) :]

                var_seq_1 = var_seq[i][:var_lcs_loc]
                var_seq_2 = var_seq[i][var_lcs_loc + len(lcs_seq) :]

                ref_seq[i : i + 1] = [ref_seq_1] + common_seq + [ref_seq_2]
                var_seq[i : i + 1] = [var_seq_1] + common_seq + [var_seq_2]

                seq_length = len(ref_seq)
                for k in range(seq_length):
                    if k + 1 > len(ref_seq):
                        break
                    if ref_seq[k] == "":
                        ref_seq[k:] = ref_seq[k + 1 :]
                        var_seq[k - 1] = var_seq[k - 1] + var_seq[k]
                        var_seq[k:] = var_seq[k + 1 :]

                    elif var_seq[k] == "":
                        var_seq[k:] = var_seq[k + 1 :]
                        ref_seq[k - 1] = ref_seq[k - 1] + ref_seq[k]
                        ref_seq[k:] = ref_seq[k + 1 :]

                nb_changes += 1

        if len(ref_seq) != len(var_seq):
            print("error spliting snv area at", self.loc)
            return 0
        self.aligned_ref_seq = ref_seq
        self.aligned_var_seq = var_seq

        self.local_alignment()

        return nb_changes

    def explain_by_snv(self, snv_set, indel_only=True):

        loc = self.loc
        ref_seq = self.aligned_ref_seq.copy()
        var_seq = self.aligned_var_seq.copy()

        i = -1
        nb_changes = 0
        rel_loc = 0
        while True:
            i += 1
            if i >= len(ref_seq) - 1:
                break
            rel_loc += len(ref_seq[i])
            if (
                len(ref_seq[i]) == len(var_seq[i])
                or min(len(ref_seq[i]), len(var_seq[i])) > 5
            ):
                continue

            coverage_range = [loc + rel_loc - len(ref_seq[i]), loc + rel_loc]

            for snv_case in snv_set:
                if (
                    loc <= snv_case.loc <= coverage_range[1]
                    and snv_case.loc + len(snv_case.ref_base) - 1 <= coverage_range[1]
                ):

                    if snv_case.type == "deletion":
                        pos_diff = snv_case.loc - (loc + rel_loc - len(ref_seq[i])) - 1
                        temp_pos = i + pos_diff
                        if pos_diff < 0:
                            temp_seq = "".join(ref_seq[temp_pos : i + 1])
                        elif pos_diff == 0:
                            temp_seq = ref_seq[i]
                        else:
                            # print(2)
                            temp_seq = ref_seq[i][pos_diff:]

                        if not temp_seq[: len(snv_case.ref_base)] == snv_case.ref_base:
                            continue

                        # length_diff = len(temp_seq)-len(snv_case.ref_base)

                        if 2 * len(snv_case.ref_base) < len(
                            "".join(ref_seq[temp_pos : i + 1])
                        ):
                            continue

                        if pos_diff <= 0:
                            ref_seq[temp_pos : i + 1] = [snv_case.ref_base] + list(
                                temp_seq[len(snv_case.ref_base) :]
                            )

                            if len(ref_seq[temp_pos : i + 1]) != len(
                                [snv_case.ref_base]
                                + list(temp_seq[len(snv_case.ref_base) :])
                            ):
                                length_diff = len(
                                    [snv_case.ref_base]
                                    + list(temp_seq[len(snv_case.ref_base) :])
                                ) - len(ref_seq[temp_pos : i + 1])
                                ref_seq[i : i + length_diff + 1] = [
                                    "".join(ref_seq[i : i + length_diff + 1])
                                ]
                        else:
                            ref_seq[i : i + pos_diff] = list(ref_seq[i][:pos_diff]) + [
                                snv_case.ref_base
                            ]
                            var_seq[i : i + pos_diff + 1] = list(
                                "".join(var_seq[i : i + pos_diff + 1])
                            )

                        nb_changes += 1
                        if len(ref_seq) != len(var_seq):
                            # print('error at explaining snv area', loc)
                            return 0
                        self.aligned_ref_seq = ref_seq
                        self.aligned_var_seq = var_seq
                        return nb_changes

                    elif snv_case.type == "insertion":
                        pos_diff = snv_case.loc - (loc + rel_loc - len(ref_seq[i])) - 1
                        temp_pos = i + pos_diff
                        if pos_diff < 0:
                            temp_seq = "".join(var_seq[temp_pos : i + 1])
                        elif pos_diff == 0:
                            temp_seq = var_seq[i]
                        else:
                            temp_seq = var_seq[i][pos_diff:]
                            # print(2)

                        if not temp_seq[: len(snv_case.var_base)] == snv_case.var_base:
                            continue

                        # length_diff = len(temp_seq)-len(snv_case.var_base)

                        if 2 * len(snv_case.var_base) < len(
                            "".join(var_seq[temp_pos : i + 1])
                        ):
                            continue

                        if pos_diff <= 0:
                            var_seq[temp_pos : i + 1] = [snv_case.var_base] + list(
                                temp_seq[len(snv_case.var_base) :]
                            )

                            if len(var_seq[temp_pos : i + 1]) != len(
                                [snv_case.var_base]
                                + list(temp_seq[len(snv_case.var_base) :])
                            ):
                                length_diff = len(
                                    [snv_case.var_base]
                                    + list(temp_seq[len(snv_case.var_base) :])
                                ) - len(var_seq[temp_pos : i + 1])
                                var_seq[i : i + length_diff + 1] = [
                                    "".join(var_seq[i : i + length_diff + 1])
                                ]
                        else:
                            var_seq[i : i + pos_diff] = list(var_seq[i][:pos_diff]) + [
                                snv_case.var_base
                            ]
                            ref_seq[i : i + pos_diff + 1] = list(
                                "".join(ref_seq[i : i + pos_diff + 1])
                            )

                        nb_changes += 1
                        if len(ref_seq) != len(var_seq):
                            # print('error at explaining snv area', loc)
                            return 0
                        self.aligned_ref_seq = ref_seq
                        self.aligned_var_seq = var_seq
                        return nb_changes

        return 0

    def recognise_dense_isnv(self):

        ref_seq = self.aligned_ref_seq.copy()
        var_seq = self.aligned_var_seq.copy()

        i = -1
        nb_changes = 0

        diff_seq = []
        for i in range(len(ref_seq)):
            if ref_seq[i] != var_seq[i]:
                diff_seq.append(i)
        nb_snv = len(diff_seq)
        if 1 < nb_snv < 3 and min(diff_seq) + nb_snv - 1 == max(diff_seq):
            temp_ref_seq_list = list(
                "".join(ref_seq[min(diff_seq) : max(diff_seq) + 1])
            )
            temp_var_seq_list = list(
                "".join(var_seq[min(diff_seq) : max(diff_seq) + 1])
            )

            if len(temp_ref_seq_list) > len(temp_var_seq_list):
                s_dim = len(temp_var_seq_list) // 2
                p_dim = len(temp_var_seq_list) - s_dim
                temp_ref_seq = temp_ref_seq_list[:p_dim] + temp_ref_seq_list[-s_dim:]
                temp_ref_seq[p_dim - 1] = "".join(temp_ref_seq_list[p_dim - 1 : -s_dim])

                ref_seq[min(diff_seq) : max(diff_seq) + 1] = temp_ref_seq
                var_seq[min(diff_seq) : max(diff_seq) + 1] = temp_var_seq_list
            else:
                s_dim = len(temp_ref_seq_list) // 2
                p_dim = len(temp_ref_seq_list) - s_dim
                temp_var_seq = temp_var_seq_list[:p_dim] + temp_var_seq_list[-s_dim:]
                temp_var_seq[p_dim - 1] = "".join(temp_var_seq_list[p_dim - 1 : -s_dim])

                ref_seq[min(diff_seq) : max(diff_seq) + 1] = temp_ref_seq_list
                var_seq[min(diff_seq) : max(diff_seq) + 1] = temp_var_seq

            nb_changes = nb_snv
            self.aligned_ref_seq = ref_seq
            self.aligned_var_seq = var_seq
            return nb_changes

        while True:
            i += 1
            if i >= len(ref_seq):
                break

            var_length = len(var_seq[i])
            ref_length = len(ref_seq[i])
            nb_snv = min(var_length, ref_length)

            if 1 < nb_snv < 3 and var_length != ref_length:

                temp_ref_seq = list(ref_seq[i])[: nb_snv - 1] + [
                    ref_seq[i][nb_snv - 1 :]
                ]
                temp_var_seq = list(var_seq[i])[: nb_snv - 1] + [
                    var_seq[i][nb_snv - 1 :]
                ]

                nb_changes = nb_snv
                ref_seq[i : i + 1] = temp_ref_seq
                var_seq[i : i + 1] = temp_var_seq

        self.aligned_ref_seq = ref_seq
        self.aligned_var_seq = var_seq
        return nb_changes

    def recognise_adjacent_snp(self):

        ref_seq = self.aligned_ref_seq.copy()
        var_seq = self.aligned_var_seq.copy()

        i = -1
        nb_changes = 0

        while True:
            i += 1
            if i >= len(ref_seq):
                break

            var_length = len(var_seq[i])
            ref_length = len(ref_seq[i])

            if min(var_length, ref_length) > 1 and var_length == ref_length:
                temp_ref_seq = list(ref_seq[i])
                temp_var_seq = list(var_seq[i])
                diff_seq = [
                    ref_base
                    for ref_base, var_base in zip(temp_ref_seq, temp_var_seq)
                    if ref_base != var_base
                ]
                if len(diff_seq) > 5:
                    continue

                nb_changes += 1
                ref_seq[i : i + 1] = temp_ref_seq
                var_seq[i : i + 1] = temp_var_seq

        self.aligned_ref_seq = ref_seq
        self.aligned_var_seq = var_seq
        return nb_changes

    def local_alignment(self):

        ref_seq = self.aligned_ref_seq.copy()
        var_seq = self.aligned_var_seq.copy()

        i = -1
        nb_changes = 0
        rel_los = 0
        while True:
            i += 1
            if i >= len(ref_seq):
                break
            rel_los += len(ref_seq[i])
            if len(ref_seq[i]) > len(var_seq[i]):
                k = 0
                if ref_seq[i][0] == ref_seq[i][-1]:
                    k = 1
                    for j in range(1, len(ref_seq[i])):
                        # print(ref_seq[i][-1-j],ref_seq[i-j])
                        if ref_seq[i][-1 - j] == ref_seq[i - j]:
                            k += 1
                        else:
                            break
                if k > 0:
                    temp_seq = "".join(ref_seq[i - k : i + 1])
                    ref_seq[i - k : i + 1] = [temp_seq[: len(ref_seq[i])]] + list(
                        temp_seq[len(ref_seq[i]) :]
                    )
                    nb_changes += 1

            elif len(ref_seq[i]) < len(var_seq[i]):
                k = 0
                if var_seq[i][0] == var_seq[i][-1]:
                    k = 1
                    for j in range(1, len(var_seq[i])):
                        # print(var_seq[i][-1-j],var_seq[i-j])
                        if var_seq[i][-1 - j] == var_seq[i - j]:
                            k += 1
                        else:
                            break
                if k > 0:
                    temp_seq = "".join(var_seq[i - k : i + 1])
                    var_seq[i - k : i + 1] = [temp_seq[: len(var_seq[i])]] + list(
                        temp_seq[len(var_seq[i]) :]
                    )
                    nb_changes += 1

        self.aligned_ref_seq = ref_seq
        self.aligned_var_seq = var_seq
        return nb_changes

    def recognise_isnv(self):

        if self.ref_seq == self.var_seq:
            return None, None

        ref_seq = self.aligned_ref_seq
        var_seq = self.aligned_var_seq
        loc = self.loc + 1
        shift_size = self.kmer_dim // 2

        iSNV_keys = []
        kmer_tuples = []
        i = -1
        index = 0
        var_index = 0

        while True:
            i += 1
            if i >= len(ref_seq):
                break

            var_length = len(var_seq[i])
            ref_length = len(ref_seq[i])

            if var_seq[i] != ref_seq[i] and min(var_length, ref_length) <= 1 and i >= 1:
                if ref_length > var_length:
                    snv_type = " deletion "
                    iSNV_key = (
                        str(loc + index) + snv_type + ref_seq[i] + " " + var_seq[i]
                    )
                elif ref_length == var_length:
                    snv_type = " snp "
                    iSNV_key = (
                        str(loc + index) + snv_type + ref_seq[i] + " " + var_seq[i]
                    )
                else:
                    snv_type = " insertion "
                    iSNV_key = (
                        str(loc + index) + snv_type + ref_seq[i] + " " + var_seq[i]
                    )

                pos1 = var_index
                pos1 = min(max(shift_size + 1, pos1), self.var_length - shift_size - 1)
                kmer_seq1 = self.var_seq[pos1 - shift_size : pos1 + shift_size + 1]
                kmer1 = ReadKmer(kmer_seq1, loc + pos1)

                pos2 = var_index + len(var_seq[i])
                pos2 = min(max(shift_size + 1, pos2), self.var_length - shift_size - 1)
                kmer_seq2 = self.var_seq[pos2 - shift_size : pos2 + shift_size + 1]
                kmer2 = ReadKmer(kmer_seq2, loc + pos2)

                # if len(kmer_seq1) != 21 or len(kmer_seq2) != 21:
                #    print(1)

                iSNV_keys.append(iSNV_key)
                kmer_tuples.append((kmer1, kmer2))
            elif var_seq[i] != ref_seq[i] and len(var_seq[i]) == len(ref_seq[i]):
                k = 0
                for j in range(len(var_seq[i])):
                    if var_seq[i][j] != ref_seq[i][j]:
                        k += 1
                if k > 5:
                    continue
                for j in range(len(var_seq[i])):
                    if var_seq[i][j] != ref_seq[i][j]:
                        snv_type = " snp "
                        iSNV_key = (
                            str(loc + index + j)
                            + snv_type
                            + ref_seq[i][j]
                            + " "
                            + var_seq[i][j]
                        )

                        pos1 = var_index + j
                        pos1 = min(
                            max(shift_size, pos1), self.var_length - shift_size - 1
                        )
                        kmer_seq1 = self.var_seq[
                            pos1 - shift_size : pos1 + shift_size + 1
                        ]
                        kmer1 = ReadKmer(kmer_seq1, loc + pos1)

                        pos2 = var_index + j
                        pos2 = min(
                            max(shift_size, pos2), self.var_length - shift_size - 1
                        )
                        kmer_seq2 = self.var_seq[
                            pos2 - shift_size : pos2 + shift_size + 1
                        ]
                        kmer2 = ReadKmer(kmer_seq2, loc + pos2)

                        # if len(kmer_seq1) != 21 or len(kmer_seq2) != 21:
                        #    print(1)

                        iSNV_keys.append(iSNV_key)
                        kmer_tuples.append((kmer1, kmer2))

            index += len(ref_seq[i])
            var_index += len(var_seq[i])

        return iSNV_keys, kmer_tuples


class SNVCase:
    def __init__(self, key):
        self._key = key
        self.kmer_tuples = set()
        self.hyplotypes = set()

        loc, snv_type, ref_base, var_base = key.split()
        self.loc = int(loc)
        self.type = snv_type
        self.ref_base = ref_base
        self.var_base = var_base

        self.count = 0
        self.depth = 0.001

    def __eq__(self, obj):
        if isinstance(obj, SNVCase):
            if self._key == obj._key:
                return True
            else:
                return False

    def __hash__(self):
        return hash(self._key)

    def get_count(self, count_dict):

        added_kmers = set()
        count = 0
        counts = np.array([0, 0])
        for temp_tuple in self.kmer_tuples:

            temp_count = [
                count_dict.get(kmer.get_kmer_code(), 0) for kmer in temp_tuple
            ]
            adding_kmer = temp_tuple[np.argmin(temp_count)]
            if adding_kmer not in added_kmers:
                counts += np.array(temp_count)
                count += int(np.min(temp_count))
                added_kmers.add(adding_kmer)

        if count == 0:
            for temp_tuple in self.kmer_tuples:

                temp_count = [
                    count_dict.get(kmer.get_kmer_code(), 0) for kmer in temp_tuple
                ]
                count += int(np.mean(temp_count))

        if np.sum(counts) == 0:
            for temp_tuple in self.kmer_tuples:

                temp_count = [
                    count_dict.get(kmer.get_kmer_code(), 0) for kmer in temp_tuple
                ]
                counts += np.array(temp_count)

        self.count = count
        self.counts = counts
        return self.count

    def get_count_bk(self, count_dict):
        count = 0
        for temp_tuple in self.kmer_tuples:

            temp_count = [
                count_dict.get(kmer.get_kmer_code(), 0) for kmer in temp_tuple
            ]
            count += int(np.min(temp_count))

        self.count = count
        return self.count

    def add_kmer_tuple(self, kmer_tuple, kernel, tran_table):
        kmer_code1 = encode(kmer_tuple[0].get_kmer(), kernel, tran_table)
        kmer_tuple[0].set_kmer_code(kmer_code1)
        kmer_code2 = encode(kmer_tuple[1].get_kmer(), kernel, tran_table)
        kmer_tuple[1].set_kmer_code(kmer_code2)
        self.kmer_tuples.add(kmer_tuple)

    def add_var_sequence(self, var_seq):
        self.hyplotypes.add(var_seq)

    def merge(self, obj):
        assert isinstance(obj, SNVCase)

        self.kmer_tuples.update(obj.kmer_tuples)

    def get_depth(self, reads_abun, snv_dict):
        self.depth = reads_abun[self.loc - 1]

    def find_lcs(self, snv_set, kmer_dim):
        for snv_case in snv_set:
            if abs(snv_case.loc - self.loc) < kmer_dim and snv_case.type == self.type:
                if snv_case.type == "insertion":
                    lcs, _ = longestCommonSubstr(self.var_base, snv_case.var_base)
                    if (
                        len(lcs) / max(len(self.var_base), len(snv_case.var_base))
                        > 0.67
                    ):
                        return True
                elif snv_case.type == "deletion":
                    lcs, _ = longestCommonSubstr(self.ref_base, snv_case.ref_base)
                    if (
                        len(lcs) / max(len(self.ref_base), len(snv_case.ref_base))
                        > 0.67
                    ):
                        return True

        return False


def checkAndMerge(s1, s2):
    m = min(len(s1), len(s2))
    for i in range(m, 0, -1):
        if s1[-i:] == s2[:i]:
            return s1 + s2[i:]
    return False


def max_filter(a, n):
    temp_array = np.zeros(len(a) + n)
    temp_array[n - 1 : -1] = a
    pooled = np.zeros(len(a))
    for i in range(len(a)):
        window = temp_array[i : i + n]
        pooled[i] = np.max(window)
    return pooled


def median_filter(a, n):
    shift_size = n // 2
    temp_array = a.copy()
    pooled = a.copy()
    for i in range(shift_size, len(a) - shift_size):
        window = temp_array[i - shift_size : i + shift_size]
        pooled[i] = np.median(window)
    return pooled.astype(int)


def correct_depth(base_depth, snv_depth, kmer_dim):
    temp_array = np.zeros(len(snv_depth) + kmer_dim)
    temp_array[: len(snv_depth)] = snv_depth
    pooled = np.zeros(len(base_depth))
    for i in range(len(base_depth)):
        window = temp_array[i : i + kmer_dim]
        pooled[i] = max(base_depth[i], np.max(window))
    return pooled


def correct_depth_old(a, b, threshold=0.9):
    corrected_array = a.copy()
    diff_rate = abs(b - a) / (a + 0.001)
    corrected_array[diff_rate < threshold] = b[diff_rate < threshold]
    return corrected_array


def get_border_dim(word1: str, word2: str) -> int:
    p_dim = 0
    s_dim = 0
    for i in range(min(len(word1), len(word2))):
        if word1[i] == word2[i]:
            p_dim += 1
        else:
            break
    for i in range(min(len(word1), len(word2))):
        if word1[-1 - i] == word2[-1 - i]:
            s_dim += 1
        else:
            break
    return p_dim, s_dim


def cal_median_freq(var_freq_list):
    if len(var_freq_list) == 1:
        return np.median(var_freq_list), 0

    var_freq_array = np.array(var_freq_list)
    # var_freq = var_freq_array[:-1].mean()/4350.3/0.818/5
    if len(var_freq_array) % 2 == 0:
        var_freq = np.median(var_freq_array)
    else:
        var_freq = np.median(var_freq_array[:-1])

    # var_var = abs(np.mean(var_freq_array)-np.median(var_freq_array))/4350.3/0.818/5

    if var_freq != 0:
        var_var = (var_freq_array.max() - var_freq_array.min()) / var_freq
        # var_var = var_freq_array[:-1].std()/var_freq_array[:-1].mean()
        var_var = np.around(var_var, 2)
    else:
        var_var = 0

    return var_freq, var_var


def sort_by_count(source_list, column):
    temp_list = []
    for item in source_list:
        temp_count = cal_median_freq(item[column])
        temp_item = [temp_count]
        temp_item.extend(item)
        temp_list.append(temp_item)
    temp_list.sort()
    sorted_list = [item[1:] for item in temp_list[::-1]]

    return sorted_list


def powerSetsBinary(items):
    combo_list = []
    N = len(items)
    for i in range(2**N):
        combo = []
        for j in range(N):
            if (i >> j) % 2 == 1:
                combo.append(items[j])
        # print(combo)
        combo_list.append(combo)
    return combo_list


def longestCommonSubstr(word1: str, word2: str) -> int:

    m = len(word1)
    n = len(word2)
    # dp = [[0] * (n + 1) for _ in range(m + 1)]
    dp = np.zeros([m + 1, n + 1], int)

    max_len = 0
    row = 0
    col = 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if word1[i - 1] == word2[j - 1]:
                dp[i, j] = dp[i - 1, j - 1] + 1
                if max_len < dp[i, j]:
                    max_len = dp[i, j]
                    row = i
                    col = j

    max_str = ""
    i = row
    j = col
    while i > 0 and j > 0:
        if dp[i, j] == 0:
            break
        i -= 1
        j -= 1
        max_str += word1[i]

    lcstr = max_str[::-1]

    nb_lcstr = len(np.where(dp == max_len)[0])

    return lcstr, nb_lcstr


def recognise_isnv(iSNV_area_list, iSNV_dict, kmer_dim, snv_length_limt):

    kernel = 4 ** np.array(range(kmer_dim), dtype=np.int64)

    snv_area_list = []
    snv_area_case_list = []
    for i in range(len(iSNV_area_list)):

        loc = iSNV_area_list[i][0]

        # if loc == 1:
        #    print(1)

        ref_seq = iSNV_area_list[i][1]
        var_kmers = iSNV_area_list[i][3]
        var_count = iSNV_area_list[i][4]

        snv_area_case_list.append(set())

        if len(ref_seq) > snv_length_limt:
            print("ignoring: ref too long:", len(ref_seq))
            continue

        snv_area = SNVArea(i, loc, ref_seq, "", var_kmers, kmer_dim, kmer_dim, kmer_dim)
        snv_area_list.append(snv_area)

    processed_set = set()
    while len(snv_area_list) > 0:

        snv_area = snv_area_list.pop(0)

        # if snv_area.loc == 1:
        #    print(snv_area.loc)
        if snv_area.is_simple_area:
            nb_changes = 0
            # snv_area.local_alignment()
            iSNV_keys, kmer_tuples = snv_area.recognise_isnv()

            for iSNV_key in iSNV_keys:
                if not iSNV_key or len(iSNV_key.split()) != 4:
                    continue

                # print(iSNV_key)
                new_case = SNVCase(iSNV_key)

                if snv_area not in processed_set and new_case.find_lcs(
                    iSNV_dict.values(), kmer_dim
                ):
                    # new_case.find_lcs(iSNV_dict.values(), kmer_dim)
                    nb_changes = snv_area.explain_by_snv(iSNV_dict.values(), True)

                    if nb_changes != 0 and snv_area not in processed_set:
                        snv_area_list.append(snv_area)
                        processed_set.add(snv_area)
                        break

            if nb_changes != 0:
                continue

            for iSNV_key, kmer_tuple in zip(iSNV_keys, kmer_tuples):
                if not iSNV_key or len(iSNV_key.split()) != 4:
                    continue

                if not iSNV_dict.get(iSNV_key):
                    iSNV_dict[iSNV_key] = SNVCase(iSNV_key)
                iSNV_dict[iSNV_key].add_kmer_tuple(kmer_tuple, kernel, TRAN_TABLE)
                iSNV_dict[iSNV_key].add_var_sequence(snv_area.var_seq)

        elif snv_area.ref_length == snv_area.var_length:
            nb_changes = snv_area.recognise_adjacent_snp()
            if nb_changes != 0:
                snv_area.is_simple_area = True
                snv_area_list.append(snv_area)
        else:
            nb_changes = snv_area.explain_by_snv(iSNV_dict.values(), True)
            if nb_changes == 0:
                nb_changes = snv_area.split_var_seq(5)
            if nb_changes != 0:
                snv_area.is_simple_area = True
                snv_area_list.append(snv_area)
            elif snv_area not in processed_set:
                nb_changes = snv_area.recognise_dense_isnv()
                if nb_changes != 0:
                    snv_area.is_simple_area = True
                snv_area_list.append(snv_area)
                processed_set.add(snv_area)

    return iSNV_dict


def recognise_simple_snv(snp_list, iSNV_dict, kmer_dim):
    temp_snv_dict = {}
    snv_area_case_list = []

    kernel = 4 ** np.array(range(kmer_dim), dtype=np.int64)

    for i in range(len(snp_list)):

        item = snp_list[i]
        snv_area_case_list.append(set())
        loc = item[0]

        ref_seq = item[1]
        var_seq = item[2]
        var_kmers = item[3]

        p_dim, s_dim = get_border_dim(ref_seq, var_seq)

        """
        if not hasattr(var_kmers, "__len__"):
            read_kmer_dict[loc] = [loc, var_kmers, decode(var_kmers)]
        else:
            read_kmer_dict[loc] = [loc, var_kmers[0], decode(var_kmers[0])]
        """

        snv_area = SNVArea(i, loc, ref_seq, "", var_kmers, kmer_dim, p_dim, s_dim)

        if snv_area.is_simple_area:
            pass
        elif snv_area.ref_length == snv_area.var_length:
            nb_changes = snv_area.recognise_adjacent_snp()
            if nb_changes == 0:
                continue
        else:
            nb_changes = snv_area.explain_by_snv(iSNV_dict.values(), True)
            if nb_changes == 0:
                nb_changes = snv_area.split_var_seq(5)

        iSNV_keys, kmer_tuples = snv_area.recognise_isnv()
        for iSNV_key, kmer_tuple in zip(iSNV_keys, kmer_tuples):
            if not iSNV_key or len(iSNV_key.split()) != 4:
                continue

            if not iSNV_dict.get(iSNV_key):
                if not temp_snv_dict.get(iSNV_key):
                    temp_snv_dict[iSNV_key] = SNVCase(iSNV_key)
                temp_snv_dict[iSNV_key].add_kmer_tuple(kmer_tuple, kernel, TRAN_TABLE)
                temp_snv_dict[iSNV_key].add_var_sequence(snv_area.var_seq)

    return temp_snv_dict


def write_csv(output_file_name, rows):
    f = open(output_file_name, "w")
    f.write("Position	Reference	Variation	Frequency	Depth\n")

    for row in rows:

        f.write("\t".join(row) + "\n")

    f.close()
    return True


def hyplotype2kmers(hyplotype, kmer_dim, kernel, tran_table):
    nb_kmer = len(hyplotype) - kmer_dim + 1
    kmers = []
    for i in range(nb_kmer):
        kmer_string = hyplotype[i : i + kmer_dim]
        kmers.append(encode(kmer_string, kernel, tran_table))
    # print(hyplotype, kmer_string, len(kmers))
    return kmers


def cal_connected_kmer_counts(
    high_frequency_kmer_dict, iSNV_dict, kmer_dim, kernel, tran_table
):
    kmer_set = set()
    for key in iSNV_dict.keys():
        hyplotypes = iSNV_dict[key].hyplotypes
        for hyplotype in hyplotypes:
            kmers = hyplotype2kmers(hyplotype, kmer_dim, kernel, tran_table)
            kmer_set.update(kmers)

    k = 0
    for kmer in kmer_set:
        k += high_frequency_kmer_dict.get(kmer, 0)
    return k


def output_isnv_bk(
    seq,
    iSNV_area_list,
    snp_list,
    border_list,
    qc_dir,
    high_frequency_kmer_dict,
    output_file_name,
    kmer_dim,
    reads_abun,
    snp_limit,
    downsample,
    count_threshold=5,
    average_reads_length=150,
    error_rate=0.01,
    vcf_threshold=0.9999,
    correct_factor=0.8,
    dfs_limt=1000,
):

    iSNV_dict = {}

    iSNV_area_list = sort_by_count(iSNV_area_list, -1)

    kernel = 4 ** np.array(range(kmer_dim), dtype=np.int64)

    iSNV_dict = recognise_isnv(iSNV_area_list, iSNV_dict, qc_dir, kmer_dim, dfs_limt)

    border_SNV_dict = recognise_simple_snv(border_list, iSNV_dict, kmer_dim)

    simple_SNV_dict = recognise_simple_snv(snp_list, iSNV_dict, kmer_dim)

    for key in border_SNV_dict.keys():
        if not iSNV_dict.get(key):
            iSNV_dict[key] = border_SNV_dict[key]

    for key in simple_SNV_dict.keys():
        if not iSNV_dict.get(key):
            iSNV_dict[key] = simple_SNV_dict[key]

    index_dict = {}
    pos_snv_dict = {}
    for key in iSNV_dict.keys():
        if not key:
            continue
        loc = int(iSNV_dict[key].loc) - 1
        ref_length = len(iSNV_dict[key].ref_base)

        if iSNV_dict[key].type == "snp":
            pos1 = loc
            pos2 = loc + 1
        elif (
            iSNV_dict[key].type == "deletion"
            and len(iSNV_dict[key].ref_base) <= average_reads_length
        ):
            pos1 = loc
            pos2 = loc + ref_length
        else:
            pos1 = loc
            pos2 = loc + 1

        for j in range(pos1, pos2):

            if not pos_snv_dict.get(j):
                pos_snv_dict[j] = set()
            pos_snv_dict[j].add(key)

        if not index_dict.get(loc):
            index_dict[loc] = []
        index_dict[loc].append(key)

        if loc == 23054:
            print(1)
        count = iSNV_dict[key].get_count(high_frequency_kmer_dict)

    snv_depth = np.zeros(reads_abun.shape)
    for pos in pos_snv_dict.keys():
        iSNV_keys = pos_snv_dict.get(pos)
        for key in iSNV_keys:
            # if iSNV_dict[key].type == 'snp' or iSNV_dict[key].type == 'deletion':
            snv_depth[pos] += iSNV_dict[key].count

    base_depth = max_filter(snv_depth + reads_abun, kmer_dim)

    # filter isnv cases
    filtered_isnv_dict = dict()
    filtered_pos_snv_dict = dict()
    key_list = list(index_dict.keys())
    key_list.sort()
    k = 0
    for key in key_list:

        for item in index_dict[key]:
            elems = item.split()
            loc = int(elems[0]) - 1

            if iSNV_dict[item].count <= 0:
                continue

            if base_depth[loc] * snp_limit <= count_threshold:
                if (
                    binom.cdf(iSNV_dict[item].count, base_depth[loc], error_rate)
                    < vcf_threshold
                ):
                    # if binom.pmf(iSNV_dict[item].count, base_depth[loc], 0.01)>0.001:
                    print("ignored low cdf", loc)
                    continue
            else:

                if iSNV_dict[item].count < count_threshold:
                    print("ignored low count isnv at", loc)
                    continue

            if elems[1] == "deletion":
                freq1 = iSNV_dict[item].counts[0] / base_depth[loc]
                freq2 = iSNV_dict[item].counts[1] / base_depth[loc + len(elems[2]) - 1]
                freq_limit = snp_limit * 0.5
            else:
                freq1 = iSNV_dict[item].counts[0] / base_depth[loc]
                freq2 = iSNV_dict[item].counts[1] / base_depth[loc]
                freq_limit = snp_limit * 0.9

            if freq1 < freq_limit or freq2 < freq_limit:

                if freq1 > 0.009 and freq2 > 0.009:
                    print("ignored low freq isnv at", loc, freq1, freq2)
                continue

            if elems[1] == "deletion":
                freq = iSNV_dict[item].count / max(
                    base_depth[loc], base_depth[loc + len(elems[2]) - 1]
                )
                # freq_limit = snp_limit*0.5
            else:
                # freq_limit = snp_limit*0.9
                freq = iSNV_dict[item].count / base_depth[loc]

            # print(freq)
            if freq < freq_limit:

                if freq > 0.009:
                    print("ignored low freq isnv at", loc, freq)
                continue
            freq = np.round(freq, 2)
            if freq > 1:
                print("warning: freq>1, freq=", freq, "at", loc)
                freq = 1
            filtered_isnv_dict[item] = iSNV_dict[item].count
            if not filtered_pos_snv_dict.get(loc):
                filtered_pos_snv_dict[loc] = set()
            filtered_pos_snv_dict[loc].add(item)

    # compensite reads depth
    kmer_depth = reads_abun.copy()
    for key in filtered_isnv_dict.keys():

        elems = key.split()
        loc = int(elems[0]) - 1
        ref_length = len(elems[2])

        neighbor_snv = []
        for i in range(loc + 1, loc + kmer_dim):
            if filtered_pos_snv_dict.get(i) and i != loc:
                keys = filtered_pos_snv_dict[i]
                for pos_key in keys:
                    if pos_key.split()[1] == "snp":
                        neighbor_snv.append(pos_key)

        if len(neighbor_snv) > 0:
            ref_kmer = seq[loc : loc + kmer_dim]
            neighbor_combo_list = powerSetsBinary(neighbor_snv)
            for neighbor_combo in neighbor_combo_list:
                if neighbor_combo == []:
                    continue
                var_kmer = list(ref_kmer)
                # test_kmer = list(ref_kmer)
                for neighbor_key in neighbor_combo:
                    var_loc, _, _, var_base = neighbor_key.split()
                    loc_diff = int(var_loc) - loc - 1
                    var_kmer[loc_diff] = var_base
                    # test_kmer[loc_diff] = var_base.lower()

                var_kmer_code = encode("".join(var_kmer), kernel, TRAN_TABLE)
                count = high_frequency_kmer_dict.get(var_kmer_code, 0)

                if count > kmer_depth[loc] * 2:
                    kmer_depth[loc] += count

    updated_base_depth = max_filter(snv_depth + kmer_depth, kmer_dim)

    # avoid false compensation
    for i in range(len(base_depth)):
        if base_depth[i] > updated_base_depth[i] * correct_factor:
            base_depth[i] = updated_base_depth[i]

    # calculate freq for simple snp
    snp_depth = reads_abun.copy()
    simple_snp_set = set()
    for key in filtered_isnv_dict.keys():

        elems = key.split()
        loc = int(elems[0]) - 1
        if elems[1] != "snp" or len(filtered_pos_snv_dict[loc]) != 1:
            continue

        mark = False
        for key in pos_snv_dict.get(loc):
            if filtered_isnv_dict.get(key):
                mark = True
                break
        if mark:
            continue

        shift_size = kmer_dim // 2

        neighbor_snv = []
        type_set = set()
        for i in range(loc + 1, loc + kmer_dim):
            if filtered_pos_snv_dict.get(i) and i != loc:
                keys = filtered_pos_snv_dict[i]
                for pos_key in keys:
                    neighbor_snv.append(pos_key)
                    type_set.add(pos_key.split()[1])

        if "insertion" in type_set or "deletion" in type_set:
            continue

        simple_snp_set.add(key)
        if len(neighbor_snv) == 0:
            continue

        ref_kmer = seq[loc : loc + kmer_dim]
        neighbor_combo_list = powerSetsBinary(neighbor_snv)
        for neighbor_combo in neighbor_combo_list:
            if neighbor_combo == []:
                continue
            var_kmer = list(ref_kmer)
            test_kmer = list(ref_kmer)
            for neighbor_key in neighbor_combo:
                var_loc, _, _, var_base = neighbor_key.split()
                loc_diff = int(var_loc) - loc - 1
                var_kmer[loc_diff] = var_base
                test_kmer[loc_diff] = var_base.lower()

            var_kmer_code = encode("".join(var_kmer), kernel, TRAN_TABLE)
            count = high_frequency_kmer_dict.get(var_kmer_code, 0)

            snp_depth[loc] += count

    # ouput
    k = 0
    # t0 = time.time()
    f = open(output_file_name, "w")
    f.write("Position	Reference	Variation	Frequency	Depth\n")

    for key in filtered_isnv_dict.keys():

        elems = key.split()
        loc = int(elems[0]) - 1
        depth = int(base_depth[loc])
        count = iSNV_dict[key].count

        if count <= 0:
            continue

        if elems[1] == "deletion":
            freq = count / min(base_depth[loc], base_depth[loc + len(elems[2]) - 1])
            freq = min(freq, 1)
        elif key in simple_snp_set:
            depth = snp_depth[loc] + count
            freq = count / depth
        else:
            freq = count / base_depth[loc]

        freq = np.round(freq, 2)
        if freq < snp_limit:
            continue

        if False:
            string_count = 0
            reads_string = ""
            for seq in iSNV_dict[key].hyplotypes:
                string_count = max(reads_string.count(seq), string_count)

            if string_count < snp_limit * depth / 4:
                print("fp variant at", loc, string_count, snp_limit * depth / 4)
                continue

        row = [
            elems[0],
            elems[2],
            elems[3],
            str(freq),
            str(depth * downsample),
        ]  # , str(int(local_abun))]

        f.write("\t".join(row) + "\n")
        k += 1

    f.close()
    return k, index_dict, iSNV_dict
