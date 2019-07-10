# From a .csv with PFAM annotations, filter those annotations where the PFAM model is not part of <list>
# Create <list> e.g. by seeing which PFAM model files are available.

import csv
import os
import sys

def filter_annotations(file, result_file, possible_annotations, delimiter_csv="\t", delimiter_pfam=";"):

    """
    Currently assuming that the annotation of interest is in the last column of `file`.
    """

    filtered_annotations = []
    with open(file, "r") as fh:
         reader = csv.reader(fh, delimiter=delimiter_csv)
         filtered_annotations.append(next(reader))
         for row in reader:
             annotations = [i for i in row[-1].split(delimiter_pfam) if i in possible_annotations]
             filtered_row = row[:]
             filtered_row[-1] = delimiter_pfam.join(annotations)
             filtered_annotations.append(filtered_row)

    with open(result_file, 'w') as fh:
        for row in filtered_annotations:
            fh.write(delimiter_csv.join(row)+ "\n")


def get_possible_annotations_from_dir(directory):

    files = os.listdir(directory)
    possible_annotations = [i.split(".")[0] for i in files]
    return possible_annotations


if __name__ == "__main__":
    # Update

    pfam_annotations_file = sys.argv[1]
    filtered_pfam_annotations_result_file = sys.argv[2]
    pfam_files_directory = sys.argv[3]
    pfam_files_summary = sys.argv[4]

    possible_annotations = get_possible_annotations_from_dir(pfam_files_directory)
    with open(pfam_files_summary, "w") as fh:
        fh.write("\n".join(possible_annotations))

    filter_annotations(pfam_annotations_file, filtered_pfam_annotations_result_file, possible_annotations)
