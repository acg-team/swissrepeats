import mobidb
import argparse


def read_commandline_arguments():
    parser = argparse.ArgumentParser(
        description='Download and save mobidb annotations for a given set of proteins')
    parser.add_argument('-i', '--input', type=str,
                        help='The path to the swissprot file.')
    parser.add_argument('-o', '--outdir', type=str,
                        help='The path to the output directory.')

    pars = vars(parser.parse_args())
    pars = {key: value for key, value in pars.iteritems() if value is not None}
    return pars


if __name__ == '__main__':
    pars = read_commandline_arguments()
    mobidb.get_swissprot(pars['input'], pars['outdir'])
