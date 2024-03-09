import numpy as np
import os
import time
from pyiSNV.utils import (
    Configuration,
    get_file_name,
    get_file_type,
    is_second,
    find_paired,
    ACCEPT_TYPE_LIST,
)
from pyiSNV.building_ref_db import GenomeMapBuilder
from pyiSNV.building_reads_db import SampleMapBuilder
from pyiSNV.building_connection_table import ConnectMapBuilder
from pyiSNV.output_isnv_area import VariantLocator
from pyiSNV.recognise_isnv import VariantCaller
from multiprocessing import Pool


def detect_variant(
    R1_sequence_file: str,
    R2_sequence_file: str,
    ref_file: str,
    config: Configuration,
    output_path: str,
):

    assert get_file_type(R1_sequence_file) in ACCEPT_TYPE_LIST
    if output_path == "":
        output_path = get_file_name(R1_sequence_file)

    T0 = time.time()

    ref_builder = GenomeMapBuilder(config)

    ref_builder.build_ref_db(ref_file)

    print("1/5 ref db building complete")
    T1 = time.time()
    print("time using: ", np.round(T1 - T0, 2), "seconds")
    print("-" * 30)

    sample_mapper = SampleMapBuilder(config)

    high_frequency_kmer_dict, reads_kmers, reads_depth = sample_mapper.build_reads_db(
        R1_sequence_file, R2_sequence_file, ref_builder
    )

    print("2/5 reads db building complete")
    T1 = time.time()
    print("time using: ", np.round(T1 - T0, 2), "seconds")
    print("-" * 30)

    connection_builder = ConnectMapBuilder(config)

    connection_mat, unique_long_kmers, long_kmer_counts = (
        connection_builder.build_connection_table(
            high_frequency_kmer_dict, reads_kmers, ref_builder
        )
    )

    print("3/5 connection table building complete")
    T1 = time.time()
    print("time using: ", np.round(T1 - T0, 2), "seconds")
    print("-" * 30)

    variant_locator = VariantLocator(config)

    iSNV_area_list = variant_locator.output_iSNV_area(
        high_frequency_kmer_dict,
        connection_mat,
        unique_long_kmers,
        long_kmer_counts,
        reads_depth,
        ref_builder,
    )

    print("4/5 iSNV position calculated")
    T1 = time.time()
    print("time using: ", np.round(T1 - T0, 2), "seconds")
    print("-" * 30)

    variant_caller = VariantCaller(config)

    nb_iSNV = variant_caller.output_isnv(
        iSNV_area_list,
        high_frequency_kmer_dict,
        reads_depth,
        output_path,
        ref_builder,
        sample_mapper,
    )

    print("total cases of iSNV:", nb_iSNV)

    print("5/5 iSNV table output")
    T1 = time.time()
    print("time using:", np.round(T1 - T0, 2), "seconds")
    print("-" * 30)

    reads_kmer_count = len(reads_kmers[reads_kmers != -1]) // 2
    genome_kmer_count = sample_mapper.nb_genome_kmer // 2

    print("total kmer:", reads_kmer_count)
    print("counts of genome kmer:", genome_kmer_count)
    print("filtered-in kmer:", variant_locator.hf_kmer_count)
    print("connected high freq kmer:", variant_locator.connected_kmer_count)
    print("variant kmer:", variant_caller.connected_kmer_counts)
    print("-" * 30)

    genome_kmer_rate = genome_kmer_count / reads_kmer_count
    excluded_kmer_count = (
        reads_kmer_count - genome_kmer_count - variant_caller.connected_kmer_counts
    )
    excluded_kmer_rate = excluded_kmer_count / reads_kmer_count
    hf_kmer_rate = variant_locator.hf_kmer_count / reads_kmer_count
    variant_kmer_rate = variant_caller.connected_kmer_counts / reads_kmer_count
    connection_rate = min(
        1, variant_locator.connected_kmer_count / variant_locator.hf_kmer_count
    )
    utilization_rate = min(
        1, variant_caller.connected_kmer_counts / variant_locator.connected_kmer_count
    )

    print("genome kmer proportion:", f"{genome_kmer_rate * 100:.2f}%")
    print("excluded kmer proportion:", f"{excluded_kmer_rate * 100:.2f}%")
    print("filtered-in kmer proportion:", f"{hf_kmer_rate * 100:.2f}%")
    print("filtered-in kmer connection rate:", f"{connection_rate * 100:.2f}%")
    print("connected kmer utlization rate:", f"{utilization_rate * 100:.2f}%")
    print("variant kmer proportion:", f"{variant_kmer_rate * 100:.2f}%")

    print("-" * 30)
    return (
        reads_kmer_count,
        genome_kmer_count,
        genome_kmer_rate,
        excluded_kmer_count,
        excluded_kmer_rate,
        hf_kmer_rate,
        variant_kmer_rate,
        connection_rate,
        utilization_rate,
    )


def process_folder(data_dir, ref_file, output_dir, config):
    files = os.listdir(data_dir)
    stat_dict = dict()
    for file in files:
        if is_second(file):
            continue
        if get_file_type(file) not in ACCEPT_TYPE_LIST:
            continue
        R1_file = data_dir + file
        R2_file = find_paired(R1_file)
        output_file = output_dir + get_file_name(R1_file) + ".csv"

        (
            reads_kmer_count,
            genome_kmer_count,
            genome_kmer_rate,
            excluded_kmer_count,
            excluded_kmer_rate,
            hf_kmer_rate,
            variant_kmer_rate,
            connection_rate,
            utilization_rate,
        ) = detect_variant(R1_file, R2_file, ref_file, config, output_file)

        stat_dict[file] = [
            reads_kmer_count,
            genome_kmer_count,
            genome_kmer_rate,
            excluded_kmer_count,
            excluded_kmer_rate,
            hf_kmer_rate,
            variant_kmer_rate,
            connection_rate,
            utilization_rate,
        ]

    return stat_dict


def parallel_process_folder(
    data_dir, ref_file, output_dir, config, max_process, interval=1
):

    pool = Pool(max_process)

    files = os.listdir(data_dir)
    for file in files:
        if is_second(file):
            continue
        if get_file_type(file) not in ACCEPT_TYPE_LIST:
            continue
        R1_file = data_dir + file
        R2_file = find_paired(R1_file)
        output_file = os.path.join(output_dir, get_file_name(R1_file) + ".csv")
        pool.apply_async(
            func=detect_variant, args=(R1_file, R2_file, ref_file, config, output_file)
        )
        time.sleep(interval)

    pool.close()
    pool.join()
    return True
