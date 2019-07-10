#!python

import argparse
import csv
import datetime
import logging
import logging.config
import os
import pickle
import sys

from pyfaidx import Fasta

from tral import configuration
from tral.paths import config_file
from tral.repeat_list import repeat_list
from tral.sequence import sequence
from tral.hmm import hmm
from tral.hmm import hmm_io

logging.config.fileConfig(config_file("logging.ini"))
LOG = logging.getLogger('root')


# Shift REPEAT_LIST_TAG to configuration file and rename?
REPEAT_LIST_TAG = "all"
DE_NOVO_ALL_TAG = "denovo_all"
PFAM_ALL_TAG = "pfam_all"
DE_NOVO_TAG = "denovo"
PFAM_TAG = "pfam"
DE_NOVO_REFINED_TAG = "denovo_refined"
DE_NOVO_FINAL_TAG = "denovo_final"
FINAL_TAG = "final"

POST_FILTER_TAG = "post_denovo_filtered"
POST_FINAL_TAG = "post_final"


def refilter_all_tandem_repeats(
                            annotation_pickle_dir,
                            result_file_serialized,
                            format="tsv",
                            **kwargs):
    """ Perform refilter_tandem_repeats an all TR pickles in a dir.
    """
    first = True
    for file in os.listdir(annotation_pickle_dir):
        if file.endswith(".pickle"):
            print(file)
            file = os.path.join(annotation_pickle_dir, file)
            if first:
                refilter_tandem_repeats(annotation_pickle=file, result_file_serialized=result_file_serialized, format=format, result_file_serialized_open="w")
                first = False
            else:
                refilter_tandem_repeats(annotation_pickle=file, result_file_serialized=result_file_serialized, format=format, result_file_serialized_open="a")


def refilter_tandem_repeats(
        annotation_pickle,
        result_file_serialized,
        result_pickle=None,
        format="tsv",
        result_file_serialized_open="w",
        **kwargs):
    ''' Refilter TRs from an annotation pickle.

    Save results as .csv

     Args:
         annotation_pickle (str): Path to the pickle file containing ``Repeat`` annotations


     Raises:
        Exception: If the pickle ``annotation_pickle`` cannot be loaded

    '''

    try:
        with open(annotation_pickle, 'rb') as fh:
            tr_annotations = pickle.load(fh)
    except:
        raise Exception(
            "Cannot load annotation_pickle: {}".format(annotation_pickle))

    try:
        if not os.path.isdir(os.path.dirname(result_pickle)):
            os.makedirs(os.path.dirname(result_pickle))
    except:
        LOG.info("Could not create path to result_pickle directory: {}".format(result_pickle))

    try:
        if not os.path.isdir(os.path.dirname(result_file_serialized)):
            os.makedirs(os.path.dirname(result_file_serialized))
    except:
        raise Exception(
            "Could not create path to result_file_serialized directory: {}".format(
                os.path.dirname(result_file_serialized)))

    # According to what criterion to we post filter?
    criterion_filter_order = {
            "func_name": "none_overlapping", "overlap": (
                "common_ancestry", None), "l_criterion": [
                ("pvalue", "phylo_gap01"), ("divergence", "phylo_gap01")]}


    for protein_id, iS in tr_annotations.items():

        iS.d_repeatlist[POST_FILTER_TAG] = iS.get_repeatlist(DE_NOVO_FINAL_TAG).filter(**criterion_filter_order)
        if len(iS.d_repeatlist[POST_FILTER_TAG].repeats) != len(iS.d_repeatlist[DE_NOVO_FINAL_TAG].repeats):
            print(- len(iS.d_repeatlist[POST_FILTER_TAG].repeats) + len(iS.d_repeatlist[DE_NOVO_FINAL_TAG].repeats))

        iS.set_repeatlist(
            iS.get_repeatlist(POST_FILTER_TAG) +
            iS.get_repeatlist(PFAM_TAG),
            POST_FINAL_TAG)

    dResults = tr_annotations

    # 6.a Save results as pickle
    if result_pickle:
        with open(result_pickle, 'wb') as fh:
            pickle.dump(dResults, fh)

    # 6.b Save serialized results
    with open(result_file_serialized, result_file_serialized_open) as fh_o:

        if result_file_serialized_open == "w":
            if format == 'tsv':
                header = [
                    "ID",
                    "begin",
                    "pvalue",
                    "l_effective",
                    "n",
                    "n_effective",
                    "TRD",
                    "model",
                    "MSA",
                    ]
            fh_o.write("\t".join(header))

        for iS in dResults.values():
            for iTR in iS.get_repeatlist(POST_FINAL_TAG).repeats:
                if format == 'tsv':
                    try:
                        data = [
                            str(i) for i in [
                                iS.name,
                                iTR.begin,
                                iTR.pvalue("phylo_gap01"),
                                iTR.l_effective,
                                iTR.n,
                                iTR.n_effective,
                                iTR.TRD,
                                iTR.model,
                                " ".join(iTR.msa),
                                ]
                            ]
                    except:
                        print(iTR)
                        raise Exception(
                            "(Could not save data for the above TR.)")

                fh_o.write("\n" + "\t".join(data))

    print("DONE")


def main():

    pars = read_commandline_arguments()


    if pars["method"] == "refilter_tandem_repeats":
        refilter_tandem_repeats(
            pars["input"],
            pars["output_serialized"],
            pars["output"],
        )
    elif pars["method"] == "refilter_all_tandem_repeats":
        refilter_all_tandem_repeats(
            pars["input"],
            pars["output_serialized"],
        )
    else:
        raise Exception("method {} not implemented".format(pars["method"]))


def read_commandline_arguments():
    parser = argparse.ArgumentParser(
        description='Process tandem repeat detection options')
    parser.add_argument('method', metavar='method_name', type=str,
                        help='The name of the method to be executed.')
    parser.add_argument('-i', '--input', type=str,
                        help='The path to the input file.')
    parser.add_argument('-o', '--output', type=str, required=False,
                        help='The path to the output file.')
    parser.add_argument('-os', '--output_serialized',
                        help='The path to the serialized output file.')


    pars = vars(parser.parse_args())
    pars = {key: value for key, value in pars.items() if value is not None}
    return pars


if __name__ == "__main__":
    main()

    """
Example execution:
python tandem_repeat_post_annotation_filtering_scripts.py refilter_tandem_repeats -i ../../results/tr_annotations/split/frag782.pickle -o ~/Downloads/frag782.pickle -os ~/Downloads/frag782.csv
python tandem_repeat_post_annotation_filtering_scripts.py refilter_all_tandem_repeats -i ../../results/tr_annotations/split -os ~/Downloads/tr_annotations_post_filter.csv

    """