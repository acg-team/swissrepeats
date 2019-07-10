import argparse
import csv
import logging
import os
import pickle
from Bio import SeqIO
from collections import defaultdict

LOG = logging.getLogger(__name__)

pickle_filename = 'nucleotide_counts.pickle'
AMINOACIDS = 'ABCDEFGHIKLMNOPQRSTUVWXYZ'

import homorepeat_expected_distribution
import nucleotide_counts

def split_coords(coords):
    return [tuple([int(y) for y in x.split(':')]) \
            for x in coords.split(';') if x]


def count_homorepeats_in_seq(result_dict, sequence, extreme=10):
    extremes = []
    current = sequence[0]
    count = 1
    for char in sequence[1:]:
        if char == current:
            count += 1
        else:
            result_dict[current][count] += 1
            if count >= extreme:
                extremes.append((current, count))
            current = char
            count = 1
    if count >= extreme:
        extremes.append((current, count))
    result_dict[current][count] += 1
    return extremes


def count_homorepeats_in_fasta(file):
    """
    In all sequence of a fasta file, count the frequency of homorepeat stretches of all lengths.
    """
    result_dict = {aa: defaultdict(int) for aa in AMINOACIDS}
    with open(file, 'r') as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seq = str(record.seq)
            extremes = count_homorepeats_in_seq(result_dict, seq)
    return result_dict


def count_homorepeats_in_annotated_fasta(file, regions, kingdoms, ignore_id,
                                         extreme_result_file="extreme.csv"):
    """
    In all annotated regions of sequences of a fasta file, count the frequency of homorepeat stretches of all lengths.
    """
    result_dict = {kingdom: {aa: defaultdict(int) for aa in AMINOACIDS} for kingdom in set(kingdoms.values())}
    with open(file, 'r') as fh, open(extreme_result_file, "w") as rfh:
        rfh.write("ID,aa,extreme_count\n")
        for record in SeqIO.parse(fh, 'fasta'):
            id = str(record.id.split('|')[1])
            if id in ignore_id:
                continue
            if id in regions:
                seq = str(record.seq)
                for region in regions[id]:
                    if seq[region[0]-1:region[1]] == "":
                        continue
                    try:
                        extremes = count_homorepeats_in_seq(result_dict[kingdoms[id]], seq[region[0]-1:region[1]])
                        for extreme in extremes:
                            rfh.write("{},{},{}\n".format(id, extreme[0], str(extreme[1])))
                    except:
                        print(id)
                        import pdb; pdb.set_trace()
    return result_dict


def write_dict(dict):
    return "{" + ", ".join(["{}: {}".format(str(key), str(value)) for key, value in dict.items()]) + "}"


def extract_and_write_homorepeat_counts(min_disorder_lengh, path_trunk, regions_file, ignore_id,
                                        fasta_file="data/seq/uniprot_sprot.fasta",
                                        results_trunk="results/empirical_and_expected_homorepeat_counts/swiss_prot_homorepeat_frequencies/empirical"):

    """
    ignore_id = ['C9JRZ8', 'Q02E40', 'M3XYN3', 'A6NES4', 'B8H0N0', 'D4AY02', 'B8GYG4', 'E2ASG4', 'Q4PE27', 'Q5LRI7', 'O82377', 'Q0DK78', 'Q4P821', 'Q66619', 'P9WFN3', 'Q7FAD5', 'P22105', 'Q4PIC4', 'F1MYR9', 'Q4P4U7', 'Q8NH95', 'Q4PEX7', 'Q4PD88', 'E1BM58', 'A6NFD8', 'Q4PIP8', 'Q03835', 'F1LW30', 'C1CW06', 'P84718', 'Q7Z3U7', 'P54099', 'P68929', 'Q9SUC3', 'C0SC25', 'Q9H0D2', 'C0S902', 'Q68772', 'Q4P553', 'Q9DCZ4', 'P9WJ45', 'Q5LUP9', 'Q9SGX9', 'P68928', 'V5IPE4', 'K7LFJ0', 'Q4PC85', 'D4B327', 'D3Z3C6', 'Q66653', 'P07059', 'P67225', 'Q9BXW9', 'Q94BT9', 'Q4P4G2', 'Q4PA86', 'Q4P6D5', 'P0CF96', 'P9WQM1', 'C1FZ15', 'Q016E7', 'Q66618', 'Q8X6R3', 'Q9BKZ9', 'L0T5V6', 'P9WQM0', 'Q4PB95', 'Q9FG29', 'G0S7R3', 'Q5LRZ7', 'G5BQH4', 'Q4P3I9', 'Q9AAH2', 'P54763', 'C1GX39', 'C1GVA2', 'F4KAF2', 'G0SC29', 'P11553', 'P30598', 'Q9FLH0', 'P54307', 'Q4P0N6', 'H7BZ55', 'A4D1Z8', 'O60936', 'Q9L268', 'Q5N6V8', 'C0S4K1', 'A5F5X1', 'D4AVJ0', 'G0S2X1', 'Q4P4C1', 'P0CAT6', 'A0QZ13', 'Q4QEB3', 'P23759', 'P20866', 'Q8RWR2']
    """

    results_trunk = os.path.join(path_trunk, results_trunk)
    if not os.path.exists(results_trunk):
        os.makedirs(results_trunk)

    results_file = os.path.join(results_trunk, "homorepeat_counter.txt")
    print("extract_and_write_homorepeat_counts results file: {}".format(results_file))
    result_pickle = os.path.join(results_trunk, "homorepeat_counter.pickle")
    extremes_region = os.path.join(results_trunk, "homorepeat_extremes_region.csv")
    extremes_rest = os.path.join(results_trunk, "homorepeat_extremes_rest.csv")

    regions, kingdoms = nucleotide_counts.read_regions(os.path.join(path_trunk, regions_file), min_disorder_lengh)
    result_dict_region = count_homorepeats_in_annotated_fasta(os.path.join(path_trunk, fasta_file), regions, kingdoms, ignore_id, extremes_region)

    # WARNING! This inverse includes only those sequences present in the other set - this does not necessarily include all sequences.
    regions_inverse = nucleotide_counts.get_inverse_annotation(regions)

    result_dict_rest = count_homorepeats_in_annotated_fasta(os.path.join(path_trunk, fasta_file), regions_inverse, kingdoms, ignore_id, extremes_rest)

    results = {"region": result_dict_region, "rest": result_dict_rest}

    with open(os.path.join(path_trunk, result_pickle), "wb") as fh:
        pickle.dump(results, fh)

    with open(os.path.join(path_trunk, results_file), "w") as fh:
        for region, data in results.items():
            fh.write("{}\n".format(region))
            for kingdom, data in data.items():
                fh.write("{}\n".format(kingdom))
                for aa in AMINOACIDS:
                    fh.write("{}: {}\n".format(aa, write_dict(data[aa])))

    return results


def calculate_expected_homorepeat_distribution(amino_acid_counts,
                                               region_annotation_file,
                                               results_file,
                                               ignore_id,
                                               k_max=10):

    if os.path.isfile(results_file):
        with open(results_file, "rb") as fh:
            return pickle.load(fh)

    ## Load empirical region lengths
    regions, kingdoms = nucleotide_counts.read_regions(region_annotation_file, min_disorder_lengh)

    # Assemble the region lengths for all kingdoms.
    ns = {kingdom: [] for kingdom in set(kingdoms.values())}
    for id, region in regions.items():
        if id not in ignore_id:
            ns[kingdoms[id]] += [i[1]-i[0]+1 for i in region]

    LOG.debug(ns.keys())

    ## Calculate empiricial amino acid frequencies from empirical_aa_counts
    amino_acid_total_counts = {kingdom: sum(counts.values()) for kingdom, counts in amino_acid_counts.items()}
    amino_acid_frequencies = {kingdom: {aa: count/amino_acid_total_counts[kingdom] for aa, count in counts.items()} for kingdom, counts in amino_acid_counts.items()}

    results = {}
    for kingdom, regions in ns.items():
        print(kingdom)
        print(amino_acid_frequencies[kingdom])
        results[kingdom] = homorepeat_expected_distribution.expected_number_runs_cumulated_over_ns(ns=regions, ps=amino_acid_frequencies[kingdom], k_max=k_max)

    with open(results_file, "wb") as fh:
        pickle.dump(results, fh)

    return results


def calculate_expected_homorepeat_distribution_unbound(amino_acid_counts,
                                                       results_pickle,
                                                       results_txt,
                                                       delimiter=",",
                                                       k_max=10):


    ## Calculate empiricial amino acid frequencies from amino_acid_counts
    amino_acid_total_counts = {kingdom: sum(counts.values()) for kingdom, counts in amino_acid_counts.items()}
    amino_acid_frequencies = {kingdom: {aa: count/amino_acid_total_counts[kingdom] for aa, count in counts.items()} for kingdom, counts in amino_acid_counts.items()}

    results = {kingdom: {aa: {} for aa in AMINOACIDS} for kingdom in amino_acid_frequencies.keys()}
    for kingdom, aa_data in amino_acid_frequencies.items():
        for aa, p in aa_data.items():
            for k in range(1, k_max+1):
                results[kingdom][aa][k] = homorepeat_expected_distribution.expected_frequency_runs(k, p)

    with open(results_pickle, "wb") as fh:
        pickle.dump(results, fh)

    with open(results_txt, "w") as fh:
        fh.write(delimiter.join(["kingdom", "aa", "n", "expected_frequency"]) + "\n")
        for kingdom, data_k in results.items():
            for aa, data_aa in data_k.items():
                for k, f in data_aa.items():
                    fh.write(delimiter.join([kingdom, aa, str(k), str(f)]) + "\n")

    return results



def write_expected_and_empirical_homorepeat_distributions(empirical_homorepeat_counts,
                                                          expected_homorepeat_counts,
                                                          results_file,
                                                          delimiter=","):

    with open(results_file, "w") as fh:
        fh.write(delimiter.join(["Kingdom", "aa", "n", "count", "type"]) + "\n")

        for kingdom, expected_counts in expected_homorepeat_counts.items():
            empirical_counts = empirical_homorepeat_counts[kingdom]
            for aa in AMINOACIDS:
                for n, count in expected_counts[aa].items():
                    fh.write(delimiter.join([kingdom, aa, str(n), str(count), "expected"]) + "\n")
                    if n in empirical_counts[aa]:
                        fh.write(delimiter.join([kingdom, aa, str(n), str(empirical_counts[aa][n]), "empirical"]) + "\n")
                    else:
                        fh.write(delimiter.join([kingdom, aa, str(n), "0", "empirical"]) + "\n")




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--regions', required=True,
                        help='The file containing the regions definitions.')
    parser.add_argument('-t', '--tag', required=True,
                        help='e.g. swissprot_disorder_inverse_kingdomwise.')
    parser.add_argument('-p', '--path_trunk', required=True,
                        help='The trunk path prepended to all other paths.')
    args = parser.parse_args()
    tag = args.tag
    path_trunk = args.path_trunk
    region_annotation_file = os.path.join(path_trunk, args.regions)

    #path_trunk = "/Users/elkeschaper/Documents/SIB/Tandem_repeats/stockholm_university/swissrepeat"
    #path_trunk = "/proj/bioinfo/users/x_oxasa/sp"

    #region_annotation_file = os.path.join(path_trunk, "results/disorder_annotations/mobidb_regions.csv")
    #region_annotation_file = os.path.join(path_trunk, "data/uniprot/swissprot_annotations_for_aa_counts.csv")

    min_disorder_lengh = 1
    k_max = 50

    # Some sequences are ignored, as disorder annotations were build on outdated Swiss-Prot sequences. To be fixed later:
    #ignore_id = ['C9JRZ8', 'Q02E40', 'M3XYN3', 'A6NES4', 'B8H0N0', 'D4AY02', 'B8GYG4', 'E2ASG4', 'Q4PE27', 'Q5LRI7', 'O82377', 'Q0DK78', 'Q4P821', 'Q66619', 'P9WFN3', 'Q7FAD5', 'P22105', 'Q4PIC4', 'F1MYR9', 'Q4P4U7', 'Q8NH95', 'Q4PEX7', 'Q4PD88', 'E1BM58', 'A6NFD8', 'Q4PIP8', 'Q03835', 'F1LW30', 'C1CW06', 'P84718', 'Q7Z3U7', 'P54099', 'P68929', 'Q9SUC3', 'C0SC25', 'Q9H0D2', 'C0S902', 'Q68772', 'Q4P553', 'Q9DCZ4', 'P9WJ45', 'Q5LUP9', 'Q9SGX9', 'P68928', 'V5IPE4', 'K7LFJ0', 'Q4PC85', 'D4B327', 'D3Z3C6', 'Q66653', 'P07059', 'P67225', 'Q9BXW9', 'Q94BT9', 'Q4P4G2', 'Q4PA86', 'Q4P6D5', 'P0CF96', 'P9WQM1', 'C1FZ15', 'Q016E7', 'Q66618', 'Q8X6R3', 'Q9BKZ9', 'L0T5V6', 'P9WQM0', 'Q4PB95', 'Q9FG29', 'G0S7R3', 'Q5LRZ7', 'G5BQH4', 'Q4P3I9', 'Q9AAH2', 'P54763', 'C1GX39', 'C1GVA2', 'F4KAF2', 'G0SC29', 'P11553', 'P30598', 'Q9FLH0', 'P54307', 'Q4P0N6', 'H7BZ55', 'A4D1Z8', 'O60936', 'Q9L268', 'Q5N6V8', 'C0S4K1', 'A5F5X1', 'D4AVJ0', 'G0S2X1', 'Q4P4C1', 'P0CAT6', 'A0QZ13', 'Q4QEB3', 'P23759', 'P20866', 'Q8RWR2']
    ignore_id = []

    # The following file was created with `extract_and_write_homorepeat_counts(annotation_file, min_disorder_lengh)`
    empirical_homorepeat_count_trunk = os.path.join(path_trunk, "results", "empirical_and_expected_homorepeat_counts", "swiss_prot_homorepeat_frequencies", "empirical", tag)
    empirical_homorepeat_count_pickle = os.path.join(empirical_homorepeat_count_trunk, "homorepeat_counter.pickle")
    LOG.info("empirical_homorepeat_count_pickle: {}".format(empirical_homorepeat_count_pickle))

    if not os.path.isfile(empirical_homorepeat_count_pickle):
        extract_and_write_homorepeat_counts(path_trunk=path_trunk, regions_file=region_annotation_file,
                                            min_disorder_lengh=min_disorder_lengh, ignore_id=ignore_id,
                                            results_trunk=empirical_homorepeat_count_trunk)

    with open(empirical_homorepeat_count_pickle, "rb") as fh:
        empirical_homorepeat_counts = pickle.load(fh)
    empirical_homorepeat_counts = empirical_homorepeat_counts["region"] # default: region

    empirical_aa_counts_pickle = \
        os.path.join(path_trunk, "results/amino_acid_frequencies", tag, "nucleotide_counts.pickle")
    LOG.info("empirical_aa_counts_pickle: {}".format(empirical_aa_counts_pickle))
    amino_acid_counts = nucleotide_counts.read(empirical_aa_counts_pickle)

    expected_homorepeat_distribution_file = \
        os.path.join(path_trunk, "results/empirical_and_expected_homorepeat_counts/swiss_prot_homorepeat_frequencies/expected", "{}.pickle".format(tag))

    expected_homorepeat_counts = \
        calculate_expected_homorepeat_distribution(region_annotation_file=region_annotation_file,
                                                   amino_acid_counts=amino_acid_counts,
                                                   results_file=expected_homorepeat_distribution_file,
                                                   k_max=k_max,
                                                   ignore_id=ignore_id)

    expected_and_empirical_homorepeat_distribution_file = \
        os.path.join(path_trunk, "results/empirical_and_expected_homorepeat_counts/swiss_prot_homorepeat_frequencies/empirical_and_expected", "{}.csv".format(tag))
    write_expected_and_empirical_homorepeat_distributions(empirical_homorepeat_counts=empirical_homorepeat_counts,
                                                          expected_homorepeat_counts=expected_homorepeat_counts,
                                                          results_file=expected_and_empirical_homorepeat_distribution_file)