__author__ = 'Alexendar Perez'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""CRISPR Specificity Correction

National Cancer Institute, National Institutes of Health, United States of America
Developer: Alexendar R. Perez M.D., Ph.D
Primary Investigator: Joana A. Vidigal Ph.D
Laboratory: Vidigal Laboratory, 2020

"""

#################
#               #
#   Libraries   #
#               #
#################

import sys
import argparse
import numpy as np
import pickle as pl
import pandas as pd
from math import sqrt

from pyearth import Earth
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

from pkg_resources import resource_exists, resource_filename

#########################
#                       #
#   Auxillary Function  #
#                       #
#########################

def arg_parser():
    parser = argparse.ArgumentParser()
    hamming_data = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('-i','--infile',help='absolute filepath to input file',required=True)
    hamming_data.add_argument('-l',dest='library',help='avana, brunello, geckov1, geckov2, tkov: default is None',default=None)
    hamming_data.add_argument('-g',dest='Hamming',help='absolute filepath to CSC Hamming pickle file: default is None',default=None)

    args = parser.parse_args()
    in_file = args.infile
    library = args.library
    genome = args.Hamming

    return in_file,library,genome

def writeout(df, hamming_string_dict, outfile):
    """write out

    :param df: input dataframe that is made from input file
    :param hamming_string_dict: dictionary object created from hamming string pickle
    :param outfile: opened outfile object
    :return: output file

    """
    df_v = np.asarray(df)
    for i in range(len(df_v)):
        grna = df_v[i][0]
        for j in df_v[i]:
            outfile.write('%s,' % j)
        try:
            for jj in hamming_string_dict[grna]:
                outfile.write('%s,' % jj)
            specific, h0 = float(hamming_string_dict[grna][0]), float(hamming_string_dict[grna][1])
            if specific >= 0.16 and h0 == 1:
                c = 'above GuideScan specificity threshold'
            else:
                c = 'below GuideScan specificity threshold'
            outfile.write('%s\n' % c)

        except KeyError:
            sys.stderr.write('\n%s not found in selected library: passing\n' % grna)
            outfile.write('\n%s not present in library\n' % grna)

def specificity_metrics(outdir, filename, df, hamming_string_dict):
    """

    :param outdir: absolute filepath to output directory
    :param filename: name of input file to be used as part of output filename
    :param df: pandas dataframe with first column as gRNA
    :param hamming_string_dict: CSC onboard dictionary object with key as gRNA and value as Hamming metrics
    :return: file with gRNA and specificity metrics

    """
    with open('%s/%s_CSC_gRNA_Hamming_neighborhood.csv' % (outdir, filename), 'w') as outfile:
        outfile.write('%s,%s,%s,%s,%s,%s,%s\n' % (
            'gRNA', 'specificity', 'h0', 'h1', 'h2', 'h3', 'classification'))
        writeout(df, hamming_string_dict, outfile)
    sys.stdout.write('write out complete\n%s/%s_CSC_gRNA_Hamming_neighborhood.csv' % (outdir, filename))

def csc(df, hamming_string_dict, outdir, filename):
    """CRISPR Specificity Correction

    :param df: pandas dataframe with first column as gRNA and second column as logFC/metric
    :param hamming_string_dict: CSC onboard dictionary object with key as gRNA and value as Hamming metrics
    :param outdir: absolute filepath to output directory
    :param filename: name of input file to be used as part of output filename
    :return: CSC adjustment

    """
    # MARS compatible file
    df_mars_lst = []
    df_v = np.asarray(df)
    for i in range(len(df_v)):
        row_lst = []
        grna, metric = df_v[i][0], df_v[i][1]
        try:
            metric = float(metric)
        except ValueError:
            sys.stdout.write('WARNING: encountered %s which is not float compatible, skipping\n' % metric)
            continue
        row_lst.append(grna)
        try:
            for jj in hamming_string_dict[grna]:
                row_lst.append(jj)
            row_lst.append(metric)
            df_mars_lst.append(row_lst)
        except KeyError:
            sys.stdout.write('\n%s not found in selected library: passing\n' % grna)
            continue

    df = pd.DataFrame(df_mars_lst, columns=['gRNA', 'specificity', 'h0', 'h1', 'h2', 'h3', 'original_value'])

    # exclude infinte specificity non-target gRNAs
    df = df[df['h0'] != 0]

    # isolate pertinent confounder variables
    df_confounders = df[['specificity', 'h0', 'h1', 'h2', 'h3']]

    # knots
    knots = df['original_value'].quantile([0.25, 0.5, 0.75, 1])

    # training and testing data
    train_x, test_x, train_y, test_y = train_test_split(df_confounders, df['original_value'], test_size=0.10,
                                                        random_state=1)

    # Fit an Earth model
    model = Earth(feature_importance_type='gcv')
    try:
        model.fit(train_x, train_y)
    except ValueError:
        sys.stdout.write('\nValue Error encountered. Model unable to be trained. Exiting CSC Novo\n')
        model_processed = 'F'
        sys.stdout.write('training input x data\n %s\ntraining input y data\n %s\n' % (train_x,train_y))
        return model_processed

    # Print the model
    print(model.trace())
    print(model.summary())
    print(model.summary_feature_importances())

    # Plot the model
    y_hat = model.predict(test_x)

    # calculating RMSE values
    rms1 = sqrt(mean_squared_error(test_y, y_hat))
    print('\n\nRMSE on Predictions\n\n')
    print(rms1)

    # calculating R^2 for training
    print('\n\nR^2 on Training Data\n\n')
    print(model.score(train_x, train_y))

    # calculating R^2 for testing
    print('\n\nR^2 on Testing Data\n\n')
    print(model.score(test_x, test_y))

    # write out model metrics
    with open('%s/csc_model_metrics_%s.txt' % (outdir, filename), 'w') as outfile:
        outfile.write('%s\n%s\n%s\nRMSE on Predictions\n%s' % (
            model.trace(), model.summary(), model.summary_feature_importances(), rms1))

    if rms1 <= 1.0:

        #model processed
        model_processed = 'T'

        # full data prediction
        df['earth_adjustment'] = model.predict(df_confounders)

        # CSC correction
        df['earth_corrected'] = df['original_value'] - df['earth_adjustment']

        # main write out
        df.to_csv('%s/csc_output_%s_earth_patched.csv' % (outdir, filename))

        # pickle write out
        model_file = open('%s/csc_output_%s_earth_model.pl' % (outdir, filename), 'wb')
        pl.dump(model, model_file)
        model_file.close()

        sys.stdout.write('\nCSC adjustment complete\n')
        sys.stdout.write('\nCSC output files written to %s\n' % outdir)
        return model_processed

    else:
        sys.stdout.write('\nCSC adjustment not computed as model residual mean squared error exceeds 1.0\n')
        model_processed = 'F'
        return model_processed

def read_in(in_file):
    """multiple attempt read in for generic file

    :param in_file: absolute filepath to input file
    :return: opened file, classification of opening method

    """
    classification = '.csv'
    if '\t' in open(in_file).readline():
        classification = '.txt'

    try:
        infile = pd.read_excel(in_file)
        sys.stdout.write('file read in as Excel\n')

    except:

        try:
            if classification == '.csv':
                infile = pd.read_csv(in_file)
                sys.stdout.write('file read in as csv\n')
            else:
                infile = pd.read_csv(in_file, sep='\t')
                sys.stdout.write('file read in as txt\n')

        except:
            infile = pd.DataFrame(open(in_file, 'r'))
            sys.stdout.write('file read in with python open function and cast as pandas DataFrame\n')

    return infile

def csc_processing(in_file, hamming_string_dict):
    """control function that assessed if CSC adjustment/model deployed or if specificity metrics only are given

    :param in_file: absolute filepath to input file
    :param hamming_string_dict: dictionary object with gRNA as key and hamming string as value
    :return: CSC adjustment or specificity metric output

    """
    # read in file
    df = read_in(in_file)
    filename, outdir = in_file.split('/')[-1].split('.')[0], '/'.join(in_file.split('/')[:-1])
    columns, rows = len(df.columns), df.shape[0]

    # ensure columns named correctly
    if columns > 1:
        sys.stdout.write(
            '\n%s columns detected\nfirst two columns will be used\n---column one = gRNA---\n---column two = value---\n' % columns)
        df = df.iloc[:, 0:2]
        df.columns = ['gRNA', 'original_value']

        model_processed = csc(df, hamming_string_dict, outdir, filename)
        if model_processed == 'T':
            pass
        else:
            specificity_metrics(outdir, filename, df, hamming_string_dict)

    elif columns == 1:
        sys.stdout.write('\nfile determined to have only one column\n---column one = gRNA---\n')
        specificity_metrics(outdir, filename, df, hamming_string_dict)

    else:
        sys.stdout.write('\nfile determined to have no columns. Unable to process\n')
        sys.exit(1)

def load_pickle(f):
    """load pickle file and generate dictionary

    :param f: absolute filepath to CSC library pickle files
    :return: dictionary object (Pandas)

    """
    with open(f, 'rb') as infile:
        pickle_dataframe = pl.load(infile,encoding='latin1')

        try:
            pickle_dictionary = pickle_dataframe.set_index('gRNA').to_dict()
            return pickle_dictionary

        except AttributeError:
            if type(pickle_dataframe) == dict:
                sys.stdout.write('\n%s is a dictionary object\n' % f)
                pickle_dictionary = pickle_dataframe
                return pickle_dictionary

            else:
                sys.stderr.write('\n%s is incompatible pickle file\nHave pickle file be dictionary with gRNA as key and specificity string as value\n' % f)
                sys.exit(1)

def file_load(infile):
    """input parameter selections

    :param infile: name of screen

    :return: filepath for Hamming and correction factor pickles for library

    """
    if infile == 'avana':
        infile_h = 'screen_models/Hamming/avana_patched_Hamming_string.pl'
        h = resource_filename(__name__, infile_h)
        return h

    elif infile == 'brunello':
        infile_h = 'screen_models/Hamming/brunello_patch_format_screen_Hamming_string.pl'
        h = resource_filename(__name__, infile_h)
        return h

    elif infile == 'geckov1':
        infile_h = 'screen_models/Hamming/geckov1_patch_format_screen_Hamming_string.pl'
        h = resource_filename(__name__, infile_h)
        return h

    elif infile == 'geckov2':
        infile_h = 'screen_models/Hamming/geckov2_patch_format_screen_Hamming_string.pl'
        h = resource_filename(__name__, infile_h)
        return h

    elif infile == 'tkov':
        infile_h = 'screen_models/Hamming/tkov_patch_format_screen_Hamming_string.pl'
        h = resource_filename(__name__, infile_h)
        return h

    elif infile == 'example_grna_logfc':
        infile_h = 'screen_models/examples/avana_patched_sample_gRNA_lognorm_lnfc.csv'
        h = resource_filename(__name__, infile_h)
        return h

    elif infile == 'example_grna':
        infile_h = 'screen_models/examples/avana_patched_sample_gRNA.csv'
        h = resource_filename(__name__, infile_h)
        return h

    else:
        sys.stderr.write('%s not a recognized screen\n' % infile)

def processing(in_file,screen,classification):
    """core processing function

    :param in_file: absolute filepath to input file
    :param screen: string value corresponding to screen name
    :param classification: deploy lite or novo
    :return:
    """

    # supported screens
    screen = screen.lower()
    support_screens = ['avana', 'brunello', 'geckov1', 'geckov2']

    if classification == 'l':
        # ensure strings all lowercase
        sys.stdout.write('\nCSC Lite deployed\n')
    elif classification == 'g':
        sys.stdout.write('\nCSC Novo deployed\n')

    # check if support screen queried
    if screen in support_screens:
        sys.stdout.write('loading %s screen data\n' % screen)
        h = file_load(screen)

        # load pickle and generate dictionaries
        hamming_dict = load_pickle(h)

        # translate hamming string
        hamming_string_dict = {}
        for key in hamming_dict['Hamming_string'].keys():
            float_casted = [float(i) for i in hamming_dict['Hamming_string'][key].split('_')]
            hamming_string_dict[key] = float_casted

        csc_processing(in_file, hamming_string_dict)

    elif screen == 'example':
        if in_file == 'example_grna_logfc':
            in_file = file_load('example_grna_logfc')
        elif in_file == 'example_grna':
            in_file = file_load('example_grna')
        else:
            sys.stderr.write('ENTER\n"csc_process -i example_grna_logfc -l example"\nOR\n"csc_process -i example_grna -l example"\n')
            sys.exit(1)

        sys.stdout.write('Example\n')
        h = file_load('avana')

        # load pickle and generate dictionaries
        hamming_dict = load_pickle(h)

        # translate hamming string
        hamming_string_dict = {}
        for key in hamming_dict['Hamming_string'].keys():
            hamming_string_dict[key] = hamming_dict['Hamming_string'][key].split('_')

        csc_processing(in_file, hamming_string_dict)

    else:
        if screen == 'tkov':
            sys.stdout.write('tkov screen\n')
            h = file_load('tkov')

            # load pickle and generate dictionaries
            hamming_dict = load_pickle(h)

        else:
            sys.stdout.write('\nscreen selection of %s is novel; will attempt load into memory\n' % screen)

            # load pickle and generate dictionaries
            hamming_dict = load_pickle(screen)

        # translate hamming string
        hamming_string_dict = {}
        try:
            for key in hamming_dict['Hamming_string'].keys():
                float_casted = [float(i) for i in hamming_dict['Hamming_string'][key].split('_')]
                hamming_string_dict[key] = float_casted
        except KeyError:
            for key in hamming_dict.keys():
                float_casted = [float(i) for i in hamming_dict[key].split('_')]
                hamming_string_dict[key] = float_casted

        csc_processing(in_file, hamming_string_dict)

#####################
#                   #
#   Main Function   #
#                   #
#####################

def main():

    # user inputs
    in_file,library,genome = arg_parser()

    if library:
        screen = library
        classification = 'l'
    else:
        screen = genome
        classification = 'g'

    # processing
    processing(in_file,screen,classification)

    # user end message
    sys.stdout.write('\nprocessing complete\n')

if __name__ == '__main__':
    main()

