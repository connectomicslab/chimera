#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2019/3/14 10:00
# @Author  : Yasser Alemán-Gómez
# @Site    : Centre Hospitalier Universitaire Vaudois, Lausanne, Switzerland (CHUV)
# @File    : chimera.py
# @Software: VSCODE

# Importing libraries
import os
import sys
import re
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from bids import BIDSLayout
from glob import glob
from operator import itemgetter
from datetime import datetime
import argparse
import csv
import json
import subprocess
import scipy.ndimage as sc
import concurrent.futures
import time
import shutil

import clabtoolkit.parcellationtools as clabparc #_parc_tsv_table, _tissue_seg_table
import clabtoolkit.misctools as clabmisc  #_rgb2hex, _hex2rgb, _printprogressbar
import clabtoolkit.freesurfertools as clabfs #_read_fscolorlut

import utils.utils as ut # SmartFormatter, _load_parctype_json, _search_in_atlas, _print_availab_parcels

def _build_args_parser():

    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=52)

    from argparse import ArgumentParser

    p = argparse.ArgumentParser(formatter_class=ut.SmartFormatter, description='\n Help \n')

    requiredNamed = p.add_argument_group('Required arguments')
    requiredNamed.add_argument('--regions', '-r', action='store_true', required=False,
                                help="R| List of available parcellations for each supra-region. \n"
                                "\n")

    requiredNamed.add_argument('--bidsdir', '-b', action='store', required=False, metavar='BIDSDIR', type=str, nargs=1,
                                help="R| BIDs dataset folder. \n"
                                "\n")
    requiredNamed.add_argument('--derivdir', '-d', action='store', required=False, metavar='DERIVDIR', type=str, nargs=1,
                                help="R| BIDs derivative folder containing the derivatives folder. \n"
                                "\n",
                                default='None')
    requiredNamed.add_argument('--parcodes', '-p', action='store', required=False,
                                metavar='CODE', type=str, nargs=1,
                                help="R| Sequence of nine one-character identifiers (one per each supra-region).\n"
                                    " The defined supra-regions are: 1) Cortex, 2) Basal ganglia, 3) Thalamus, \n"
                                    " 4) Amygdala, 5) Hippocampus, 6) Hypothalamus, 7) Cerebellum, 8) Brainstem. \n"
                                    "\n"
                                    "Example: \n"
                                    "Parcellation code: HFMIIIFIF.\n"
                                    "   1. Cortical parcellation (H): HCP-MMP1 cortical parcellation (Glasser et al, 2016).\n"
                                    "   2. Basal ganglia parcellation (F): FreeSurfer subcortical parcellation (Fischl et al, 2002).\n"
                                    "   3. Thalamic parcellation (M): Atlas-based thalamic parcellation (Najdenovska et al, 2018).\n"
                                    "   4. Amygdala parcellation (I): Amygdala nuclei parcellation (Saygin et al, 2017).\n"
                                    "   5. Hippocampus parcellation (I): Hippocampus subfield parcellation (Iglesias et al, 2015).\n"
                                    "   6. Hypothalamus parcellation (I): Hypothalamus parcellation (Billot et al, 2020).\n"
                                    "   7. Cerebellum parcellation (F): Default FreeSurfer cerebellum segmentation.\n"
                                    "   8. Brainstem parcellation (I): Brainstem parcellation (Iglesias et al, 2015).\n"
                                    "   9. Gyral White Matter parcellation (F): WM parcellation according to the selected cortical parcellation.\n"
                                    "\n"
                                    "Use the --regions or -r options to show all the available parcellations for eact supra-region.\n"
                                    "\n")
    requiredNamed.add_argument('--nthreads', '-n', action='store', required=False, metavar='NTHREADS', type=str, nargs=1,
                                help="R| Number of processes to run in parallel. \n", default=['1'])

    requiredNamed.add_argument('--growwm', '-g', action='store', required=False, metavar='GROWWM', type=str, nargs=1,
                                help="R| Grow of GM labels inside the white matter in mm. \n", default=['2'])

    requiredNamed.add_argument('--t1s', '-t', action='store', required=False, metavar='T1FILE', type=str, nargs=1,
                                help="R| File containing the basename of the NIFTI images that will be ran. \n"
                                    "   This file is useful to tun Chimera, only, on certain T1s in case of multiple T1s \n"
                                    " for the same session.\n"
                                    " Example of this file: \n"
                                    "   sub-00001_ses-0001_run-2 \n"
                                    "   sub-00001_ses-0003_run-1\n"
                                    "   sub-00001_ses-post_acq-mprage\n"
                                    " \n", default='None')
    requiredNamed.add_argument('--force', '-f', action='store_true', required=False,
                                help="R| Overwrite the results. \n"
                                    "\n")
    p.add_argument('--verbose', '-v', action='store', required=False,
                    type=int, nargs=1,
                    help='verbosity level: 1=low; 2=debug')

    args = p.parse_args()

    if args.derivdir is None or args.bidsdir is None or args.parcodes is None :
        print('--bidsdir, --derivdir and --parcodes are REQUIRED arguments')
        sys.exit()

    bids_dir = args.bidsdir[0]
    deriv_dir = args.derivdir[0]

    if args.regions is True:
        print('Available parcellations for each supra-region:')
        ut._print_availab_parcels()
        sys.exit()

    if not os.path.isdir(deriv_dir):
        print("\n")
        print("Please, supply a valid BIDs derivative directory.")
        p.print_help()
        sys.exit()

    if not os.path.isdir(bids_dir):
        print("\n")
        print("Please, supply a valid BIDs directory.")
        p.print_help()
        sys.exit()

    return p

def _subfields2hbt(temp_ip, hipp_codes):
    """
    This function groups hippocampus subfields in head, body and tail.
    This groupping is based on the hippocampus subfields parcellation from Iglesias et al, 2015.
    For more information about the parcellation see: 
    https://surfer.nmr.mgh.harvard.edu/fswiki/HippocampalSubfieldsAndNucleiOfAmygdala

    Parameters
    ----------
    temp_ip : np.array
        Input image containing the hippocampus subfields
    hipp_codes : list
        List of codes for the hippocampus subfields

    Returns
    -------
    hbtimage : np.array
        Image containing the hippocampus subfields grouped in head, body and tail

    """

    # This function groups hippocampus subfields in head, body and tail
    hbtimage = np.zeros(np.shape(temp_ip), dtype='int16')

    # Creating Head
    bool_ind = np.in1d(temp_ip, hipp_codes[0:8])
    bool_ind = np.reshape(bool_ind, np.shape(temp_ip))
    result = np.where(bool_ind == True)
    hbtimage[result[0], result[1], result[2]] = 1

    # Creating Body
    bool_ind = np.in1d(temp_ip, hipp_codes[9:16])
    bool_ind = np.reshape(bool_ind, np.shape(temp_ip))
    result = np.where(bool_ind == True)
    hbtimage[result[0], result[1], result[2]] = 2

    # Creating Tail
    bool_ind = np.in1d(temp_ip, hipp_codes[17])
    bool_ind = np.reshape(bool_ind, np.shape(temp_ip))
    result = np.where(bool_ind == True)
    hbtimage[result[0], result[1], result[2]] = 3

    # Creating Fissure
    bool_ind = np.in1d(temp_ip, hipp_codes[18])
    bool_ind = np.reshape(bool_ind, np.shape(temp_ip))
    result = np.where(bool_ind == True)
    hbtimage[result[0], result[1], result[2]] = 4

    return hbtimage





