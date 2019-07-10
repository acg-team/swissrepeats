import mobidb
import argparse


def read_commandline_arguments():
    parser = argparse.ArgumentParser(
        description='Process disorder')
    parser.add_argument('-i', '--input', type=str,
                        help='The path to the mobidb folder.')
    parser.add_argument('-t', '--type', type=str,
                        help='Type of disorder annotation to perform')
    parser.add_argument('-m', '--method', type=str,
                        help='Disorder annotation method')
    parser.add_argument('-s', '--swissprot', type=str,
                        help='Path to the swissprot file')
    pars = vars(parser.parse_args())
    pars = {key: value for key, value in pars.iteritems() if value is not None}
    return pars


if __name__ == '__main__':
    pars = read_commandline_arguments()
    if pars['type'] == 'r':
        mobidb.process_disorder_regions(pars['input'], pars['swissprot'], pars['method'])
    elif pars['type'] == 'c':
        mobidb.process_disorder_coordinates(pars['input']) 
    else:
        mobidb.process_disorder(pars['input'])
