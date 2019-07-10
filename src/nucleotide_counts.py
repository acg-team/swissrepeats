import argparse
import csv
import os
import pickle
from Bio import SeqIO
from collections import defaultdict

pickle_filename = 'nucleotide_counts.pickle'
aminoacids = 'ABCDEFGHIKLMNOPQRSTUVWXYZ'

def dd():
    return defaultdict(int)

def split_coords(coords):
    return [tuple([int(y) for y in x.split(':')]) \
            for x in coords.split(';') if x]


def read_regions(region_annot_file, min_region_length):
    """
        Region annotation file format is expected as follows:
        uniprotID, disorder, taxon
        Q54XS2,1:9;19:19;21:21;450:451;505:505;508:509;511:538;652:660;765:771;,Eukaryota
        Q8FHK7,1:4;15:21;307:320;,Bacteria
    """
    regions = defaultdict(list)
    kingdoms = defaultdict(list)
    with open(region_annot_file, 'r') as f:
        reader = csv.reader(f)
        reader.__next__()
        for line in reader:
            regions[line[0]] = [i for i in split_coords(line[1]) if i[1]-i[0] + 1 >= min_region_length]
            kingdoms[line[0]] = line[2]
    assert len(regions) == len(kingdoms)
    return regions, kingdoms


def create_regions_file_from_annotations_file(annotation_input_file,
                                              regions_result_file,
                                              delimiter_input="\t",
                                              delimiter_result=",",
                                              inverse_regions=None):
    """
    Convert a Swiss-Prot annotation file to an input file as expected by `read_regions()`.

    If `inverse_regions` is not defined, the region covers the entire protein, such that amino acid frequencies for
    the entire protein can be calculated using this `nucleotide_counts` script.

    If `inverse_regions` is defined as a region file, the resulting annotation file will contain a) all proteins not
    listed in `inverse_regions` and b) the inverse regions to those defined.
    """

    if inverse_regions:
        regions = get_inverse_annotation(inverse_regions)

    with open(annotation_input_file, "r") as ifh, open(regions_result_file, "w") as rfh:
        input_reader = csv.reader(ifh, delimiter=delimiter_input, quotechar='|')
        header = next(input_reader)
        uniprotID_id = header.index("Entry")
        length_id = header.index("Length")
        superkingdom_id = header.index("Taxonomic lineage (SUPERKINGDOM)")
        rfh.write(delimiter_result.join(["uniprotID", "disorder", "taxon"]))
        for line in input_reader:
            if inverse_regions and line[uniprotID_id] in inverse_regions:
                region = regions[line[uniprotID_id]]
                # When the inverse region extended to the end of the protein, the last region element is a hoax.
                last_region = region.pop()
                if last_region[0] <= int(line[length_id]):
                    region.append((last_region[0], line[length_id]))
                region = ";".join("{}:{}".format(i[0], i[1]) for i in region)
            else:
                region = "1:{}".format(line[length_id])
            rfh.write("\n" + delimiter_result.join([line[uniprotID_id], region, line[superkingdom_id]]))



def get_fasta_data_handle(fasta_file):
    handle = SeqIO.parse(open(fasta_file, 'r'), 'fasta')
    return handle


def get_counts(fasta, regions, kingdoms):
    counts = defaultdict(dd)
    for record in fasta:
        id = str(record.id.split('|')[1])
        seq = str(record.seq)
        for region in regions[id]:
            for aa in aminoacids:
                counts[kingdoms[id]][aa] += seq.count(aa, region[0] - 1,
                                                      region[1])
    return counts


def get_inverse_annotation(regions):
    """
    Get a data structure similar to `regions`, however with the inverse regions defined.
    Assume that counting starts on 1.
    """
    regions_inverse = {}
    for id, region in regions.items():
        regions_inverse[id] = []
        start = 1
        for i in region:
            end = i[0] - 1
            if end >= start:
                regions_inverse[id].append((start, end))
            start = i[1] + 1
        regions_inverse[id].append((start, None))
    return regions_inverse


def create_folder(relative_path):
    absolute_path = os.path.abspath(relative_path)
    if not os.path.exists(absolute_path):
        os.makedirs(absolute_path)
    return absolute_path


def create_output(counts, output, delimiter=","):
    counts_simplified = dict(counts)
    for key in counts_simplified.keys():
        counts_simplified[key] = dict(counts_simplified[key])
    with open(os.path.join(output, pickle_filename), 'wb') as handle:
        pickle.dump(counts_simplified, handle)
    for key, value in counts_simplified.items():
        with open(os.path.join(output, key + '.csv'), 'w') as kingdom_file:
            kingdom_file.write(delimiter.join(aminoacids) + "\n")
            kingdom_file.write(delimiter.join([str(value[char]) for char in aminoacids]) + "\n")


def read(file):
    with open(file, 'rb') as fh:
        return pickle.load(fh)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', required=True,
                        help='path to the input fasta file')
    parser.add_argument('-a', '--annotations', required=True,
                        help='path to the annotations file')
    parser.add_argument('-o', '--output', required=True,
                        help='path to the output folder')
    args = parser.parse_args()
    fasta_file = args.fasta
    annotations_file = args.annotations
    output = create_folder(args.output)
    regions, kingdoms = read_regions(annotations_file, min_region_length=1)
    fasta = get_fasta_data_handle(fasta_file)
    create_output(get_counts(fasta, regions, kingdoms), output)
