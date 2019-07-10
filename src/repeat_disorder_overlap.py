import argparse
import csv
import os


from nucleotide_counts import read_regions
from nucleotide_counts import create_folder

pickle_filename = 'nucleotide_counts.pickle'
aminoacids = 'ARNDCEQGHILKMFPSTWYVUO'

def get_overlap(region, begin, end):
    overlap_start = begin if begin > region[0] else region[0]
    overlap_end = end if end < region[1] else region[1]
    return overlap_end - overlap_start

def get_overlap_coordinates(region, begin, end):
    overlap_start = begin if begin > region[0] else region[0]
    overlap_end = end if end < region[1] else region[1]
    return overlap_start, overlap_end


def find_overlap(regions, id, begin, end):
    counter = 0
    try:
        disordered_regions = regions[id]
    except KeyError as e:
        return counter
    for reg in disordered_regions:
        if reg[1] < begin:
            continue
        if reg[0] > end:
            break
        counter += get_overlap(reg, begin, end)
    return counter

def find_overlap_coordinates(regions, id, begin, end):
    overlap_data = []
    try:
        disordered_regions = regions[id]
    except KeyError as e:
        return overlap_data
    for reg in disordered_regions:
        if reg[1] < begin:
            continue
        if reg[0] > end:
            break
        overlap_start, overlap_end = get_overlap_coordinates(reg, begin, end)
        overlap_data = [overlap_start, overlap_end, end, reg[0], reg[1]]
    return overlap_data

def compute_overlaps(repeats_file, regions, output_folder):
    output_file = os.path.join(output, 'repeat_disorder_overlap.csv')
    with open(repeats_file, 'r') as in_handle,\
         open(output_file, 'w') as out_handle:
        reader = csv.reader(in_handle)
        writer = csv.writer(out_handle)
        header = reader.__next__()
        writer.writerow(header + ['disordered_overlap'])
        for row in reader:
            id = row[0]
            begin = int(row[2])
            length = int(float(row[4]) * float(row[6]))
            end = begin + length - 1
            disordered_count = find_overlap(regions, id, begin, end)
            writer.writerow(row + [disordered_count])

def compute_overlap_regions(repeats_file, disorder_file, output_folder):
    output_file = os.path.join(output_folder, 'repeat_disorder_overlap_regions.csv')
    with open(repeats_file, 'r') as in_handle,\
         open(output_file, 'w') as out_handle:
        reader = csv.reader(in_handle)
        writer = csv.writer(out_handle)
        header = reader.__next__()
        writer.writerow(header + ['overlap_start', 'overlap_end', 'end', 'idr_start', 'idr_end'])
        for row in reader:
            id = row[0]
            begin = int(row[2])
            length = int(float(row[4]) * float(row[6]))
            end = begin + length - 1
            overlap_data= find_overlap_coordinates(regions, id, begin, end)
            if len(overlap_data) !=0:
                writer.writerow(row + overlap_data)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--repeats', required=True,
                        help='path to the repeat description file')
    parser.add_argument('-a', '--annotations', required=True,
                        help='path to the annotations file')
    parser.add_argument('-o', '--output', required=True,
                        help='path to the output folder')
    parser.add_argument('-c', '--coors', required=False,
                        help='compute coordinates')
    args = parser.parse_args()
    repeats_file = args.repeats
    annotations_file = args.annotations
    output = create_folder(args.output)
    regions, kingdoms = read_regions(annotations_file, 30)
    if args.coors:
        compute_overlap_regions(repeats_file, regions, output)
    else:
        compute_overlaps(repeats_file, regions, output)
