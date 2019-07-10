"""
The purpose of this class is to model and handle overlapping sequence annotations.
"""
from collections import defaultdict
import pandas as pd
import pickle

class Sequence(object):

    def __init__(self, ID, length, **meta):
        self.ID = ID
        self.length = length

        self.meta = meta
        # Save all elements in kwargs as attributes.
        #for key, value in kwargs.items():
        #    if not hasattr(self, key):
        #        setattr(self, key, value)

        self.annotations = {}

        self.is_clean = False


    def add_single_annotation(self, start, end, tag):
        """

        Args:
            start int: Coordinates are human style (starting to count on 1)
            end int: Coordinates are human style (starting to count on 1)
            tag str:

        Returns:
        """
        if not tag in self.annotations:
            self.annotations[tag] = []
        self.annotations[tag].append((start, end))

        self.is_clean = False


    def add_annotations(self, annotations, tag):
        """
        Args:
            annotations list of (start, end): Coordinates are human style (starting to count on 1)
            tag str:

        Returns:
        """
        if not tag in self.annotations:
            self.annotations[tag] = []

        self.annotations[tag] += annotations

        self.is_clean = False


    def clean_sort_annotations(self):
        """
            Remove overlapping annotations.
            Sort annotations.
        """

        if self.is_clean:
            return

        annotations = {}
        for tag, anno in self.annotations.items():
            # anno = [(3,4), (4,16), (3,7), (20,24)]
            if len(anno) == 0:
                annotations[tag] = []
                continue
            anno = sorted(anno, key=lambda x: (x[0], -x[1]), reverse=True)
            result = [anno.pop()]
            while(anno):
                n = anno.pop()
                if n[0] > result[-1][1] + 1:
                    result.append(n)
                else:
                    result[-1] = (result[-1][0], n[1])
            annotations[tag] = result

        self.annotations = annotations
        self.is_clean = True



    def compute_overlap_intersection(self, tag_new, tag_1, tag_2):
        pass


    def filter_annotations_with_minimum_overlap(self, tag_new, tag_main, tag_overlapping, overlap_ratio):

        self.clean_sort_annotations()

        if (not tag_main in self.annotations) or (not tag_overlapping in self.annotations):
            return

        a_new = []
        a_main = self.annotations[tag_main]
        a_overlapping = self.annotations[tag_overlapping]


        for am in a_main:
            # 1. Calculate overlap
            overlap = 0
            for ao in a_overlapping:
                if ao[1] < am[0]:
                    continue
                if ao[0] > am[1]:
                    break
                overlap += min(ao[1], am[1]) - max(ao[0], am[0]) + 1

            # 2. Filter annotations with too little overlap.
            if overlap >= overlap_ratio * (am[1] - am[0] + 1):
                a_new.append(am)

        self.add_annotations(a_new, tag_new)



    def get_annotations_coverage(self, tag):
        """
        Computes the number of chars covered by the annotations tagged with tag, and the number of chars divided by the sequence length
        """

        if not tag in self.annotations:
            return 0, 0
        else:
            coverage = 0
            for annotation in self.annotations[tag]:
                coverage += annotation[1] - annotation[0] + 1
            return coverage, coverage/self.length



class Sequences(object):

    def __init__(self, sequences):
        self.sequences = sequences

    @classmethod
    def create_sequences_from_csv(cls, swissprot_filename, id_name="Entry", length_name="Length",
                                  superkingdom_name="Taxonomic lineage (SUPERKINGDOM)",
                                  kingdom_name="Taxonomic lineage (KINGDOM)",
                                  order_name="Taxonomic lineage (ORDER)",
                                  class_name="Taxonomic lineage (CLASS)",
                                  species_name="Taxonomic lineage (SPECIES)",
                                  delimiter="\t"):
        """
        Creates a dictionary of sequences.
        """

        df = pd.read_table(swissprot_filename, sep=delimiter)
        sequences = {entry: Sequence(ID=entry, length=length, species=species, species_class=species_class, order=order, kingdom=kingdom, superkingdom=superkingdom)
                     for entry, length, species, species_class, order, kingdom, superkingdom
                     in zip(df[id_name], df[length_name], df[species_name], df[class_name], df[order_name], df[kingdom_name], df[superkingdom_name])}

        return Sequences(sequences=sequences)

    @classmethod
    def create(cls, file):
        """
        Creates a dictionary of sequences.
        """

        with open(file, "rb") as fh:
            return pickle.load(fh)


    def write(self, file):

        with open(file, "wb") as fh:
            pickle.dump(self, fh)


    def add_tr_annotations_from_csv(self, tr_filename, tag="TR", id_name="ID", delimiter="\t"):

        tr_df = pd.read_table(tr_filename, sep=delimiter)
        tr_df['end'] = tr_df.apply (lambda row: int(row["begin"] + row["l_effective"]*row["n_effective"] - 1), axis=1)

        for entry, start, end, l_effective in zip(tr_df[id_name], tr_df["begin"], tr_df["end"], tr_df["l_effective"]):
            self.sequences[entry].add_single_annotation(start=start, end=end, tag=tag)
            if l_effective == 1:
                self.sequences[entry].add_single_annotation(start=start, end=end, tag=tag+"_homo")
            elif l_effective < 4:
                self.sequences[entry].add_single_annotation(start=start, end=end, tag=tag+"_micro")
            elif l_effective < 15:
                self.sequences[entry].add_single_annotation(start=start, end=end, tag=tag+"_mini")
            else:
                self.sequences[entry].add_single_annotation(start=start, end=end, tag=tag+"_long")


    def add_disorder_annotations_from_csv(self, disorder_filename, tag="IDR", id_name="uniprotID", coordinate_name=" disorder", delimiter="\t"):

        idr_df = pd.read_table(disorder_filename, sep=delimiter)

        for entry, coordinates in zip(idr_df[id_name], idr_df[coordinate_name]):
            if type(coordinates) != str:
                continue
            for coord in coordinates.split(";"):
                if not ":" in coord:
                    continue
                start, end = coord.split(":")
                start = int(start)
                end = int(end)
                self.sequences[entry].add_single_annotation(start=start, end=end, tag=tag)
                if end-start < 5:
                    self.sequences[entry].add_single_annotation(start=start, end=end, tag=tag+"_short")
                elif end-start <= 30:
                    self.sequences[entry].add_single_annotation(start=start, end=end, tag=tag+"_medium")
                else:
                    self.sequences[entry].add_single_annotation(start=start, end=end, tag=tag+"_long")


    def count_annotations(self, tags, groupby_meta="superkingdom"):

        results = {tag: defaultdict(list) for tag in tags}
        for sequence in self.sequences.values():
            for tag in tags:
                if tag in sequence.annotations:
                    results[tag][sequence.meta[groupby_meta]].append(len(sequence.annotations[tag]))
                else:
                    results[tag][sequence.meta[groupby_meta]].append(0)

        return results


    def filter_annotations_with_minimum_overlap(self, tag_new, tag_main, tag_overlapping, overlap_ratio=0.9):
        for sequence in self.sequences.values():
            sequence.filter_annotations_with_minimum_overlap(tag_new, tag_main, tag_overlapping, overlap_ratio=overlap_ratio)


    def get_sequence_length(self, groupby_meta="superkingdom"):

        result = defaultdict(int)
        for sequence in self.sequences.values():
            result[sequence.meta[groupby_meta]] += sequence.length
        return result


    def length_annotations(self, tags, groupby_meta="superkingdom"):

        results = {tag: {"coverage": defaultdict(list), "relative_coverage": defaultdict(list)} for tag in tags}
        for sequence in self.sequences.values():
            for tag in tags:
                if tag in sequence.annotations:
                    coverage, relative_coverage = sequence.get_annotations_coverage(tag=tag)
                    results[tag]["coverage"][sequence.meta[groupby_meta]].append(coverage)
                    results[tag]["relative_coverage"][sequence.meta[groupby_meta]].append(relative_coverage)
                else:
                    results[tag]["coverage"][sequence.meta[groupby_meta]].append(0)
                    results[tag]["relative_coverage"][sequence.meta[groupby_meta]].append(0)

        return results
