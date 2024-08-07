#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import pandas as pd
import tqdm
from pathlib import Path
import bids
from bids import BIDSLayout
import nibabel as nib
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

class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


class MyBIDs:
    def __init__(self, bids_dir: str):
        fold_list = os.listdir(bids_dir)
        subjids = []
        for it in fold_list:
            if 'sub-' in it:
                subjids.append(it)
        subjids.sort()

        bidsdict = {}

        for subjid in subjids:
            subjdir = os.path.join(bids_dir, subjid)
            sesids = []
            if os.path.isdir(subjdir):
                fold_list = os.listdir(subjdir)
                for it in fold_list:
                    if 'ses-' in it:
                        sesids.append(it)
                sesids.sort()
                bidsdict.__setitem__(subjid, sesids)
        self.value = bidsdict

    def add(self, key, value):
        print(key)
        print(value)
        self[key] = value

    def get_subjids(self, bids_dir: str):

        fold_list = os.listdir(bids_dir)
        subjids = []
        for it in fold_list:
            if 'sub-' in it:
                subjids.append(it)
        subjids.sort()
        self.subjids = subjids
        return subjids

    def get_sesids(self, bids_dir: str, subjid):

        sesids = []
        subjdir = os.path.join(bids_dir, subjid)
        if os.path.isdir(subjdir):
            fold_list = os.listdir(subjdir)
            for it in fold_list:
                if 'ses-' in it:
                    sesids.append(it)
            sesids.sort()
            self.sesids = sesids

        return sesids

# Print iterations progress
def _printprogressbar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', printend="\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printend    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledlength = int(length * iteration // total)
    bar = fill * filledlength + '-' * (length - filledlength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=printend)
    # Print New Line on Complete
    if iteration == total:
        print()


# Loading the JSON file containing the available parcellations
def _load_parctype_json():
    cwd = os.getcwd()
    serJSON = os.path.join(cwd, 'parcTypes.json')
    with open(serJSON) as f:
        data = json.load(f)

    return data

def _build_args_parser():

    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=52)

    from argparse import ArgumentParser

    p = argparse.ArgumentParser(formatter_class=SmartFormatter, description='\n Help \n')

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
        _print_availab_parcels()
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


def _print_availab_parcels(reg_name=None):

    data = _load_parctype_json()

    if reg_name is None:
        supra_keys = data.keys()
        parc_help = ''
        for sup in supra_keys:
            parc_opts = data[sup]
            parc_help = '{} "{}:\n"'.format(parc_help, sup)
            print(sup + ':')
            for opts in parc_opts:
                desc = data[sup][opts]["Atlas"]
                cita = data[sup][opts]["Citation"]
                parc_help = '{} "{}: {} {}\n"'.format(parc_help, opts, desc, cita)
                print('     {}: {} {}'.format(opts, desc, cita))
            print('')
    else:
        parc_opts = data[reg_name]
        print(reg_name + ':')
        for opts in parc_opts:
            desc = data[reg_name][opts]["Atlas"]
            cita = data[reg_name][opts]["Citation"]
            print('     {}: {} {}'.format(opts, desc, cita))
        print('')


def _get_region_features(im_parc, st_codes, st_names, st_red, st_green, st_blue):
    # Extract the code, name and RGB triplet from colorLUT file for the structures appearing in the image im_parc
    temp    = np.unique(im_parc)
    ind_non = np.nonzero(temp)
    temp    = temp[ind_non[0]]
    temp    = temp.astype(int)
    pos     = search(st_codes, temp)

    reg_codes = list(itemgetter(*pos)(st_codes))
    reg_names = list(itemgetter(*pos)(st_names))
    reg_red   = list(itemgetter(*pos)(st_red))
    reg_green = list(itemgetter(*pos)(st_green))
    reg_blue  = list(itemgetter(*pos)(st_blue))

    return reg_codes, reg_names, reg_red, reg_green, reg_blue


def _my_ismember(a, b):
    values, indices = np.unique(a, return_inverse=True)
    is_in_list = np.isin(a, b)
    idx = indices[is_in_list].astype(int)

    return values, idx

def _parc_tsv_table(codes, names, colors, tsv_filename):
    # Table for parcellation
    # 1. Converting colors to hexidecimal string
    seg_hexcol = []
    nrows, ncols = colors.shape
    for i in np.arange(0, nrows):
        seg_hexcol.append(rgb2hex(colors[i, 0], colors[i, 1], colors[i, 2]))

    bids_df = pd.DataFrame(
        {
            'index': np.asarray(codes),
            'name': names,
            'color': seg_hexcol
        }
    )
    #     print(bids_df)
    # Save the tsv table
    with open(tsv_filename, 'w+') as tsv_file:
        tsv_file.write(bids_df.to_csv(sep='\t', index=False))


# Find Structures
def _search_in_atlas(in_atlas, st_tolook, out_atlas, labmax):
    for i, v in enumerate(st_tolook):
        result = np.where(in_atlas == v)
        out_atlas[result[0], result[1], result[2]] = i + labmax + 1
    #         print('%u === %u', v, i + labmax + 1)

    labmax = labmax + len(st_tolook)
    return out_atlas, labmax


# Search the value inside a vector
def search(values, st_tolook):
    ret = []
    for v in st_tolook:
        index = values.index(v)
        ret.append(index)
    return ret


def rgb2hex(r, g, b):
    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def hex2rgb(hexcode):
    return tuple(map(ord, hexcode[1:].decode('hex')))


def _find_images_in_path(new_path, str_ext):
    # This function finds images in a folder that contain certain string in its file name
    temp_var = glob(new_path + os.path.sep + str_ext)

    return temp_var


def _subfields2hbt(temp_ip, hipp_codes):
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

def _select_t1s(t1s, t1file):

    with open(t1file) as file:
        t1s2run = [line.rstrip() for line in file]

    out_t1s = [s for s in t1s if any(xs in s for xs in t1s2run)]

    return out_t1s


def tissue_seg_table(tsv_filename):
    # Table for tissue segmentation
    # 1. Default values for tissues segmentation table
    seg_rgbcol = np.array([[172, 0, 0], [0, 153, 76], [0, 102, 204]])
    seg_codes = np.array([1, 2, 3])
    seg_names = ['cerebro_spinal_fluid', 'gray_matter', 'white_matter']
    seg_acron = ['CSF', 'GM', 'WM']

    # 2. Converting colors to hexidecimal string
    seg_hexcol = []
    nrows, ncols = seg_rgbcol.shape
    for i in np.arange(0, nrows):
        seg_hexcol.append(rgb2hex(seg_rgbcol[i, 0], seg_rgbcol[i, 1], seg_rgbcol[i, 2]))

    bids_df = pd.DataFrame(
        {
            'index': seg_codes,
            'name': seg_names,
            'abbreviation': seg_acron,
            'color': seg_hexcol
        }
    )
    # Save the tsv table
    with open(tsv_filename, 'w+') as tsv_file:
        tsv_file.write(bids_df.to_csv(sep='\t', index=False))


def read_fscolorlut(lutFile):
    # Readind a color LUT file
    fid = open(lutFile)
    LUT = fid.readlines()
    fid.close()

    # Make dictionary of labels
    LUT = [row.split() for row in LUT]
    st_names = []
    st_codes = []
    cont = 0
    for row in LUT:
        if len(row) > 1 and row[0][0] != '#' and row[0][0] != '\\\\':  # Get rid of the comments
            st_codes.append(int(row[0]))
            st_names.append(row[1])
            if cont == 0:
                st_colors = np.array([[int(row[2]), int(row[3]), int(row[4])]])
            else:
                ctemp = np.array([[int(row[2]), int(row[3]), int(row[4])]])
                st_colors = np.append(st_colors, ctemp, axis=0)
            cont = cont + 1

    return st_codes, st_names, st_colors


def _launch_annot2ind(fs_annot, ind_annot, hemi, out_dir, fullid, atlas):

    # Creating the hemisphere id
    if hemi == 'lh':
        hemicad = 'L'
    elif hemi == 'rh':
        hemicad = 'R'

    # Moving the Annot to individual space
    subprocess.run(['mri_surf2surf', '--srcsubject', 'fsaverage', '--trgsubject', fullid,
                            '--hemi', hemi, '--sval-annot', fs_annot,
                            '--tval', ind_annot],
                            stdout=subprocess.PIPE, universal_newlines=True)

    # Copying the resulting annot to the output folder
    out_annot = os.path.join(out_dir, fullid + '_hemi-' + hemicad + '_space-orig_' + atlas + '_dparc.annot')
    subprocess.run(['cp', ind_annot, out_annot], stdout=subprocess.PIPE, universal_newlines=True)

    return out_annot

def _launch_gcs2ind(fssubj_dir, fs_gcs, ind_annot, hemi, out_dir, fullid, atlas):

    # Creating the hemisphere id
    if hemi == 'lh':
        hemicad = 'L'
    elif hemi == 'rh':
        hemicad = 'R'

    # Moving the GCS to individual space
    cort_file = os.path.join(fssubj_dir, fullid, 'label', hemi + '.cortex.label')
    sph_file  = os.path.join(fssubj_dir, fullid, 'surf', hemi + '.sphere.reg')
    subprocess.run(['mris_ca_label', '-l', cort_file, fullid, hemi, sph_file,
                    fs_gcs, ind_annot], stdout=subprocess.PIPE, universal_newlines=True)

    # Copying the resulting annot to the output folder
    out_annot = os.path.join(out_dir, fullid + '_hemi-' + hemicad + '_space-orig_' + atlas + '_dparc.annot')
    subprocess.run(['cp', ind_annot, out_annot],
                    stdout=subprocess.PIPE, universal_newlines=True)

    return out_annot

def _launch_freesurfer(t1file:str, fssubj_dir:str, fullid:str):


    os.environ["SUBJECTS_DIR"] = fssubj_dir

    # Computing FreeSurfer
    subprocess.run(['recon-all', '-subjid', '-i', t1file, fullid, '-all'],
                    stdout=subprocess.PIPE, universal_newlines=True)

    return


def _launch_surf2vol(fssubj_dir, out_dir, fullid, atlas, gm_grow):

    if 'desc' not in atlas:
        atlas_str = atlas + '_desc-'
    else:
        atlas_str = atlas

    if atlas == "aparc":
        atlas_str = "atlas-desikan_desc-aparc"
    elif atlas == "aparc.a2009s":
        atlas_str = "atlas-destrieux_desc-a2009s"

    out_parc = []
    for g in gm_grow:
        out_vol = os.path.join(out_dir, fullid + '_space-orig_' + atlas_str + 'grow' + g + 'mm_dseg.nii.gz')

        if g == '0':
            # Creating the volumetric parcellation using the annot files
            subprocess.run(['mri_aparc2aseg', '--s', fullid, '--annot', atlas,
                            '--hypo-as-wm', '--new-ribbon', '--o', out_vol],
                            stdout=subprocess.PIPE, universal_newlines=True)

        else:
            # Creating the volumetric parcellation using the annot files
            subprocess.run(['mri_aparc2aseg', '--s', fullid, '--annot', atlas, '--wmparc-dmax', g, '--labelwm',
                            '--hypo-as-wm', '--new-ribbon', '--o', out_vol],
                            stdout=subprocess.PIPE, universal_newlines=True)


        # Moving the resulting parcellation from conform space to native
        raw_vol = os.path.join(fssubj_dir, fullid, 'mri', 'rawavg.mgz')
        subprocess.run(['mri_vol2vol', '--mov', out_vol, '--targ', raw_vol,
                        '--regheader', '--o', out_vol, '--no-save-reg', '--interp', 'nearest'],
                        stdout=subprocess.PIPE, universal_newlines=True)

        out_parc.append(out_vol)

    return out_parc

def _parc_conform2native(cform_mgz, nat_nii, fssubj_dir, fullid):
    # Moving the resulting parcellation from conform space to native
    raw_vol = os.path.join(fssubj_dir, fullid, 'mri', 'rawavg.mgz')
    subprocess.run(['mri_vol2vol', '--mov', cform_mgz, '--targ', raw_vol,
                    '--regheader', '--o', nat_nii, '--no-save-reg', '--interp', 'nearest'],
                    stdout=subprocess.PIPE, universal_newlines=True)

def _compute_abased_thal_parc(t1, vol_tparc, deriv_dir, pathcad, fullid, aseg_nii, out_str):

    cwd = os.getcwd()
    thal_spam = os.path.join(cwd, 'thalamic_nuclei_MIALatlas', 'Thalamus_Nuclei-HCP-4DSPAMs.nii.gz')
    t1_temp = os.path.join(cwd, 'mni_icbm152_t1_tal_nlin_asym_09c', 'mni_icbm152_t1_tal_nlin_asym_09c.nii.gz')

    # Creating spatial transformation folder
    stransf_dir = os.path.join(deriv_dir, 'ants-transf2mni', pathcad, 'anat')
    if not os.path.isdir(stransf_dir):
        try:
            os.makedirs(stransf_dir)
        except OSError:
            print("Failed to make nested output directory")

    defFile = os.path.join(stransf_dir, fullid + '_space-MNI152NLin2009cAsym_')
    if not os.path.isfile(defFile + 'desc-t12mni_1InverseWarp.nii.gz'):
        # Registration to MNI template
        subprocess.run(['antsRegistrationSyN.sh', '-d', '3', '-f', t1_temp, '-m', t1, '-t', 's',
                        '-o', defFile + 'desc-t12mni_'],
                        stdout=subprocess.PIPE, universal_newlines=True)

    mial_dir = os.path.dirname(vol_tparc)
    # Creating ouput directory
    if not os.path.isdir(mial_dir):
        try:
            os.makedirs(mial_dir)
        except OSError:
            print("Failed to make nested output directory")

    mial_thalparc = os.path.join(mial_dir, fullid + '_space-orig_desc-' + out_str +'_dseg.nii.gz')
    mial_thalspam = os.path.join(mial_dir, fullid + '_space-orig_desc-' + out_str +'_probseg.nii.gz')

    # Applying spatial transform
    subprocess.run(['antsApplyTransforms', '-d', '3', '-e', '3', '-i', thal_spam,
                    '-o', mial_thalspam, '-r', t1, '-t', defFile + 'desc-t12mni_1InverseWarp.nii.gz',
                    '-t','[' + defFile + 'desc-t12mni_0GenericAffine.mat,1]', '-n', 'Linear'],
                    stdout=subprocess.PIPE, universal_newlines=True)

    # Creating MaxProb
    _spams2maxprob(mial_thalspam, 0.05, mial_thalparc, aseg_nii, 10, 49)
    mial_thalparc = [mial_thalparc]

    return mial_thalparc

def _fs_addon_parcellations(vol_tparc, fullid, fssubj_dir, parcid, out_str):

    volatlas_dir = os.path.dirname(vol_tparc)

    # Creating ouput directory
    if not os.path.isdir(volatlas_dir):
        try:
            os.makedirs(volatlas_dir)
        except OSError:
            print("Failed to make nested output directory")

    if parcid == 'thalamus':
    # Running Thalamic parcellation
        process = subprocess.run(
            ['segmentThalamicNuclei.sh', fullid, fssubj_dir],
            stdout=subprocess.PIPE, universal_newlines=True)

        thal_mgz = os.path.join(fssubj_dir, fullid, 'mri', 'ThalamicNuclei.v12.T1.mgz')

        # Moving Thalamic parcellation to native space
        _parc_conform2native(thal_mgz, vol_tparc, fssubj_dir, fullid)

        out_parc = [vol_tparc]

    elif parcid == 'amygdala' or  parcid == 'hippocampus':
        # Running Hippocampal and Amygdala parcellation
        process = subprocess.run(
            ['segmentHA_T1.sh', fullid, fssubj_dir],
            stdout=subprocess.PIPE, universal_newlines=True)

        # Moving Hippocampal and amygdala parcellation to native space
        lh_mgz = os.path.join(fssubj_dir, fullid, 'mri', 'lh.hippoAmygLabels-T1.v21.mgz')
        lh_gz = os.path.join(volatlas_dir, fullid + '_space-orig_hemi-L_desc-' + out_str + '_dseg.nii.gz')
        _parc_conform2native(lh_mgz, lh_gz, fssubj_dir, fullid)

        rh_mgz = os.path.join(fssubj_dir, fullid, 'mri', 'rh.hippoAmygLabels-T1.v21.mgz')
        rh_gz = os.path.join(volatlas_dir, fullid + '_space-orig_hemi-R_desc-' + out_str + '_dseg.nii.gz')
        _parc_conform2native(rh_mgz, rh_gz, fssubj_dir, fullid)
        out_parc = [lh_gz, rh_gz]

    elif parcid == 'hypothalamus':

    # Running Hypothalamus parcellation
        os.system("WRITE_POSTERIORS=1")
        process = subprocess.run(
            ['mri_segment_hypothalamic_subunits', '--s', fullid, '--sd', fssubj_dir, '--write_posteriors'],
            stdout=subprocess.PIPE, universal_newlines=True)

        # Moving Hypothalamus to native space
        hypo_mgz = os.path.join(fssubj_dir, fullid, 'mri', 'hypothalamic_subunits_seg.v1.mgz')
        hypo_gz = os.path.join(volatlas_dir, fullid + '_space-orig_desc-' + out_str + '_dseg.nii.gz')
        _parc_conform2native(hypo_mgz, hypo_gz, fssubj_dir, fullid)
        out_parc = [hypo_gz]

    elif parcid == 'brainstem':

        # Running Brainstem parcellation
        # os.environ["WRITE_POSTERIORS"] = 1
        os.system("WRITE_POSTERIORS=1")
        process = subprocess.run(
            ['segmentBS.sh', fullid, fssubj_dir],
            stdout=subprocess.PIPE, universal_newlines=True)
        
        # Moving Hypothalamus to native space
        bs_mgz = os.path.join(fssubj_dir, fullid, 'mri', 'brainstemSsLabels.v12.mgz')
        bs_gz = os.path.join(volatlas_dir, fullid + '_space-orig_desc-' + out_str + '_dseg.nii.gz')
        _parc_conform2native(bs_mgz, bs_gz, fssubj_dir, fullid)
        out_parc = [bs_gz]

    return out_parc


def _spams2maxprob(spamImage:str, thresh:float=0.05, maxpName:str=None, thalMask:str=None, thl_code:int=10, thr_code:int=49):
    # ---------------- Thalamic nuclei (MIAL) ------------ #
    thalm_codesl       = np.array([1, 2, 3, 4, 5, 6, 7])
    thalm_codesr       = np.array([8, 9, 10, 11, 12, 13, 14])
    thalm_names        =  ['pulvinar', 'ventral-anterior', 'mediodorsal', 'lateral-posterior-ventral-posterior-group', 'pulvinar-medial-centrolateral-group', 'ventrolateral', 'ventral-posterior-ventrolateral-group']
    prefix             = "thal-lh-"
    thalm_namesl       = [prefix + s.lower() for s in thalm_names]
    prefix             = "thal-rh-"
    thalm_namesr       = [prefix + s.lower() for s in thalm_names]
    thalm_colorsl      = np.array([[255,   0,   0], [0, 255,   0], [255, 255, 0], [255, 123, 0], [0, 255, 255], [255, 0, 255], [0, 0, 255]])
    thalm_colorsr      = thalm_colorsl

    # ---------------- Creating output filenames ------------ #
    outDir           = os.path.dirname(spamImage)
    fname            = os.path.basename(spamImage)
    tempList         = fname.split('_')
    tempList[-1]     = 'dseg.nii.gz'
    if not maxpName:
        maxpName         =  os.path.join(outDir, '_'.join(tempList))

    tempList[-1]     = 'dseg.lut'
    lutName          =  os.path.join(outDir, '_'.join(tempList))
    tempList[-1]     = 'dseg.tsv'
    tsvName          =  os.path.join(outDir, '_'.join(tempList))

    maxlist = maxpName.split(os.path.sep)
    tsvlist = tsvName.split(os.path.sep)
    lutlist = lutName.split(os.path.sep)

    # ---------------- Creating Maximum probability Image ------------- #
    # Reading the thalamic parcellation
    spam_Ip          = nib.load(spamImage)
    affine           = spam_Ip.affine
    spam_Ip          = spam_Ip.get_fdata()
    spam_Ip[spam_Ip < thresh] = 0
    spam_Ip[spam_Ip > 1]      = 1

    # 1. Left Hemisphere
    It               = spam_Ip[:, :, :, :7]
    ind              = np.where(np.sum(It, axis=3) == 0)
    maxprob_thl      = spam_Ip[:, :, :, :7].argmax(axis=3) + 1
    maxprob_thl[ind] = 0

    if thalMask:
        Itemp        = nib.load(thalMask)
        Itemp        = Itemp.get_fdata()
        index        = np.where(Itemp != thl_code)
        maxprob_thl[index[0], index[1], index[2]] = 0

    # 2. Right Hemisphere
    It               = spam_Ip[:, :, :, 7:]
    ind              = np.where(np.sum(It, axis=3) == 0)
    maxprob_thr      = spam_Ip[:, :, :, 7:].argmax(axis=3) + 1
    maxprob_thr[ind] = 0

    if thalMask:
        index        = np.where(Itemp != thr_code)
        maxprob_thr[index[0], index[1], index[2]] = 0

    ind              = np.where(maxprob_thr != 0)
    maxprob_thr[ind] = maxprob_thr[ind] + 7

    # Saving the Nifti file
    imgcoll          = nib.Nifti1Image(maxprob_thr.astype('int16') + maxprob_thl.astype('int16'), affine)
    nib.save(imgcoll, maxpName)

    # Creating the corresponding TSV file
    _parc_tsv_table(np.concatenate((thalm_codesl, thalm_codesr)),
                    np.concatenate((thalm_namesl, thalm_namesr)),
                    np.concatenate((thalm_colorsl, thalm_colorsr)),
                    tsvName)

    # Creating and saving the corresponding colorlut table
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    luttable = ['# $Id: <BIDsDirectory>/derivatives/{} {} \n'.format('/'.join(lutlist[-5:]), date_time),
                        '# Corresponding parcellation: ',
                        '# <BIDsDirectory>/derivatives/' + '/'.join(maxlist[-5:]) ,
                        '# <BIDsDirectory>/derivatives/' + '/'.join(tsvlist[-5:]) + '\n']
    luttable.append('{:<4} {:<50} {:>3} {:>3} {:>3} {:>3} \n '.format("#No.", "Label Name:", "R", "G", "B", "A"))

    luttable.append("# Left Hemisphere. Thalamic nuclei parcellation (MIAL, Najdenovska and Alemán-Gómez et al, 2018)")
    for roi_pos, roi_name in enumerate(thalm_namesl):
        luttable.append('{:<4} {:<50} {:>3} {:>3} {:>3} {:>3}'.format(roi_pos + 1, roi_name, thalm_colorsl[roi_pos,0], thalm_colorsl[roi_pos,1], thalm_colorsl[roi_pos,2], 0))
    nright = roi_pos +1

    luttable.append('\n')

    luttable.append("# Right Hemisphere. Thalamic nuclei parcellation")
    for roi_pos, roi_name in enumerate(thalm_namesr):
        luttable.append('{:<4} {:<50} {:>3} {:>3} {:>3} {:>3}'.format(nright + roi_pos + 1, roi_name, thalm_colorsr[roi_pos,0], thalm_colorsr[roi_pos,1], thalm_colorsr[roi_pos,2], 0))

    with open(lutName, 'w') as colorLUT_f:
                colorLUT_f.write('\n'.join(luttable))

# def _build_parcellation(layout, bids_dir, deriv_dir, ent_dict, parccode):
def _build_parcellation(t1, bids_dir, deriv_dir, parccode, growwm):

    # layout = bids.BIDSLayout(bids_dir, validate=False)

    anat_dir = os.path.dirname(t1)
    t1_name = os.path.basename(t1)
    temp_entities = t1_name.split('_')[:-1]
    fullid = "_".join(temp_entities)



    if "ses" in fullid:
        ses_index = [i for i, s in enumerate(temp_entities) if 'ses' in s]
        path_cad = temp_entities[0] + os.path.sep + temp_entities[ses_index[0]]
    else:
        path_cad = temp_entities[0]

    # layout = BIDSLayout(anat_dir, validate=False)
    # ent_dict = layout.parse_file_entities(t1)

    # ######## ------------- Creating the full ID. It is used for a correct image file naming. ------------ #
    # if 'session' in ent_dict.keys():
    #     pattern_fullid = "sub-{subject}_ses-{session}_run-{run}"
    #     path_cad       = "sub-" + ent_dict["subject"] + os.path.sep + "ses-" + ent_dict["session"]
    # else:
    #     pattern_fullid = "sub-{subject}_run-{run}"
    #     path_cad       = "sub-" + ent_dict["subject"]


    # fullid             = os.path.basename(layout.build_path(ent_dict, pattern_fullid, validate=False))

    ######## ------------- Reading the parcellation dictionary ------------ #
    parcdict = _load_parctype_json()

    ######## ------------- Detecting FreeSurfer Subjects Directory  ------------ #
    fshome_dir = os.getenv('FREESURFER_HOME')
    fssubj_dir = os.path.join(deriv_dir, 'freesurfer')
    os.environ["SUBJECTS_DIR"] = fssubj_dir

    if not os.path.isdir(fssubj_dir):
        print("The freesurfer subjects directory is not inside the derivative folder.")
        sys.exit()

    ######## ------------- Reading FreeSurfer color lut table ------------ #
    lutFile = os.path.join(fshome_dir, 'FreeSurferColorLUT.txt')
    st_codes, st_names, st_colors = read_fscolorlut(lutFile)

    ######## ------------- Labelling the structures ------------ #
    # 1. ---------------- Detecting White matter  ------------------------- #
    wm_codesl = np.array([2, 5001])
    idx = search(st_codes, wm_codesl)
    wm_namesl = ['wm-lh-brain-segmented', 'wm-lh-brain-unsegmented']
    wm_colorsl = st_colors[idx]

    wm_codesr = np.array([41, 5002])
    idx = search(st_codes, wm_codesr)
    wm_namesr = ['wm-rh-brain-segmented', 'wm-rh-brain-unsegmented']
    wm_colorsr = st_colors[idx]

    cc_codes = np.array([250, 251, 252, 253, 254, 255])
    idx = search(st_codes, cc_codes)
    cc_names = np.array(st_names)[idx].tolist()
    prefix = 'wm-brain-'
    cc_names = [prefix + s.lower() for s in cc_names]
    cc_names = [s.replace('_', '-').lower() for s in cc_names]
    cc_colors = st_colors[idx]

    wm_codes = np.concatenate((wm_codesl.astype(int), wm_codesr.astype(int), cc_codes.astype(int)))
    wm_names = ['wm-brain-white_matter']
    wm_colors = np.array([[255, 255, 255]])

    # 2. ---------------- Detecting Subcortical structures (Freesurfer) ------------------------- #
    subc_codesl = np.array([11, 12, 13, 26])
    idx = search(st_codes, subc_codesl)
    subc_namesl = np.array(st_names)[idx].tolist()
    subc_namesl = [s.replace('Left-', 'subc-lh-').lower() for s in subc_namesl]
    subc_colorsl = st_colors[idx]

    subc_codesr = np.array([50, 51, 52, 58])
    idx = search(st_codes, subc_codesr)
    subc_namesr = np.array(st_names)[idx].tolist()
    subc_namesr = [s.replace('Right-', 'subc-rh-').lower() for s in subc_namesr]
    subc_colorsr = st_colors[idx]

    # 3. ---------------- Detecting Thalamic structures (Freesurfer) ------------------------- #
    thalf_codesl = np.ones((1,), dtype=int) * 10
    idx = [st_codes.index(thalf_codesl)]
    # idx                = search(st_codes, thalf_codesl)
    thalf_namesl = np.array(st_names)[idx].tolist()
    thalf_namesl = [s.replace('Left-', 'thal-lh-').lower() for s in thalf_namesl]
    thalf_colorsl = st_colors[idx]

    thalf_codesr = np.ones((1,), dtype=int) * 49
    idx = [st_codes.index(thalf_codesr)]
    # idx                = search(st_codes, thalf_codesr)
    thalf_namesr = np.array(st_names)[idx].tolist()
    thalf_namesr = [s.replace('Right-', 'thal-rh-').lower() for s in thalf_namesr]
    thalf_colorsr = st_colors[idx]

    # 3. ---------------- Detection of Thalamic nuclei (Iglesias)------------ #
    thali_codesl = np.array(
        [8103, 8104, 8105, 8106, 8108, 8109, 8110, 8111, 8112, 8113, 8115, 8116, 8117, 8118, 8119, 8120, 8121, 8122,
            8123, 8125, 8126, 8127, 8128, 8129, 8130, 8133, 8134])
    thali_codesr = np.array(
        [8203, 8204, 8205, 8206, 8208, 8209, 8210, 8211, 8212, 8213, 8215, 8216, 8217, 8218, 8219, 8220, 8221, 8222,
            8223, 8225, 8226, 8227, 8228, 8229, 8230, 8233, 8234])

    idx = search(st_codes, thali_codesl)
    thali_namesl = np.array(st_names)[idx].tolist()
    thali_colorsl = st_colors[idx]
    thali_namesl = [s.replace('_', '-') for s in thali_namesl]
    thali_namesl = [s.replace('Left-', 'thal-lh-').lower() for s in thali_namesl]

    idx = search(st_codes, thali_codesr)
    thali_namesr = np.array(st_names)[idx].tolist()
    thali_colorsr = st_colors[idx]
    thali_namesr = [s.replace('_', '-') for s in thali_namesr]
    thali_namesr = [s.replace('Right-', 'thal-rh-').lower() for s in thali_namesr]

    # 3. ---------------- Detection of Thalamic nuclei (MIAL)------------ #
    thalm_codesl = np.array([1, 2, 3, 4, 5, 6, 7])
    thalm_codesr = np.array([8, 9, 10, 11, 12, 13, 14])
    thalm_names = ['pulvinar', 'ventral-anterior', 'mediodorsal', 'lateral-posterior-ventral-posterior-group',
                    'pulvinar-medial-centrolateral-group', 'ventrolateral', 'ventral-posterior-ventrolateral-group']
    prefix = "thal-lh-"
    thalm_namesl = [prefix + s.lower() for s in thalm_names]
    prefix = "thal-rh-"
    thalm_namesr = [prefix + s.lower() for s in thalm_names]
    thalm_colorsl = np.array(
        [[255, 0, 0], [0, 255, 0], [255, 255, 0], [255, 123, 0], [0, 255, 255], [255, 0, 255], [0, 0, 255]])
    thalm_colorsr = thalm_colorsl

    # 4. ---------------- Detecting Amygdala structures (Freesurfer) ------------------------- #
    amygf_codesl = np.ones((1,), dtype=int) * 18
    idx = [st_codes.index(amygf_codesl)]

    # idx                = search(st_codes, amygf_codesl)
    amygf_namesl = np.array(st_names)[idx].tolist()
    amygf_namesl = [s.replace('Left-', 'amygd-lh-').lower() for s in amygf_namesl]
    amygf_colorsl = st_colors[idx]

    amygf_codesr = np.ones((1,), dtype=int) * 54
    idx = [st_codes.index(amygf_codesr)]

    # idx                = search(st_codes, amygf_codesr)
    amygf_namesr = np.array(st_names)[idx].tolist()
    amygf_namesr = [s.replace('Right-', 'amygd-rh-').lower() for s in amygf_namesr]
    amygf_colorsr = st_colors[idx]

    # 4. ---------------- Detecting Amygdala nuclei (Iglesias) ------------------------- #
    amygi_codesl = np.array([7001, 7003, 7005, 7006, 7007, 7008, 7009, 7010, 7015])
    idx = search(st_codes, amygi_codesl)
    amygi_namesl = np.array(st_names)[idx].tolist()
    amygi_colorsl = st_colors[idx]

    amygi_codesr = amygi_codesl
    amygi_namesr = amygi_namesl
    amygi_colorsr = amygi_colorsl

    prefix = "amygd-lh-"
    amygi_namesl = [prefix + s.lower() for s in amygi_namesl]
    amygi_namesl = [s.replace('_', '-') for s in amygi_namesl]

    prefix = "amygd-rh-"
    amygi_namesr = [prefix + s.lower() for s in amygi_namesr]
    amygi_namesr = [s.replace('_', '-') for s in amygi_namesr]

    # 5. ---------------- Detecting Hippocampus structures (Freesurfer) ------------------------- #
    hippf_codesl = np.ones((1,), dtype=int) * 17
    idx = [st_codes.index(hippf_codesl)]

    # idx                = search(st_codes, hippf_codesl)
    hippf_namesl = np.array(st_names)[idx].tolist()
    hippf_namesl = [s.replace('Left-', 'hipp-lh-').lower() for s in hippf_namesl]
    hippf_colorsl = st_colors[idx]

    hippf_codesr = np.ones((1,), dtype=int) * 53
    idx = [st_codes.index(hippf_codesr)]

    # idx                = search(st_codes, hippf_codesr)
    hippf_namesr = np.array(st_names)[idx].tolist()
    hippf_namesr = [s.replace('Right-', 'hipp-rh-').lower() for s in hippf_namesr]
    hippf_colorsr = st_colors[idx]

    # 5. ---------------- Detecting Hippocampus nuclei (Iglesias) ------------------------- #
    hippi_codesl = np.array(
        [203, 233, 235, 237, 239, 241, 243, 245, 211, 234, 236, 238, 240, 242, 244, 246, 212, 226, 215])
    idx = search(st_codes, hippi_codesl)
    hippi_namesl = np.array(st_names)[idx].tolist()
    hippi_colorsl = st_colors[idx]

    hippi_codesr = hippi_codesl
    hippi_namesr = hippi_namesl
    hippi_colorsr = hippi_colorsl

    prefix = "hipp-lh-"
    hippi_namesl = [prefix + s.lower() for s in hippi_namesl]
    hippi_namesl = [s.replace('_', '-') for s in hippi_namesl]

    prefix = "hipp-rh-"
    hippi_namesr = [prefix + s.lower() for s in hippi_namesr]
    hippi_namesr = [s.replace('_', '-') for s in hippi_namesr]

    # 5. ---------------- Detecting Hippocampus nuclei and grouping in Head, Body and Tail (Iglesias) ------------------------- #
    hipph_codesl = np.array([1, 2, 3, 4])
    hipph_namesl = ['hipp-lh-hippocampus-head', 'hipp-lh-hippocampus-body', 'hipp-lh-hippocampus-tail',
                    'hipp-lh-hippocampus-fissure']
    hipph_colorsl = np.array([[255, 0, 0], [0, 255, 0], [255, 255, 0], [255, 123, 0]])

    hipph_codesr = np.array([1, 2, 3, 4])
    hipph_namesr = ['hipp-rh-hippocampus-head', 'hipp-rh-hippocampus-body', 'hipp-rh-hippocampus-tail',
                    'hipp-rh-hippocampus-fissure']
    hipph_colorsr = np.array([[255, 0, 0], [0, 255, 0], [255, 255, 0], [255, 123, 0]])

    # 6. ---------------- Segmenting the hypothalamus using the VentralDC (Connectomics Lab) ------------ #
    # Ventral DC Left
    vdcf_codel = 28
    vdcf_namesl = ["vdc-lh-ventraldc"]
    vdcf_colorsl = np.array([[165, 42, 42]])

    # Ventral DC Right
    vdc_coder = 60
    vdcf_namesr = ["vdc-rh-ventraldc"]
    vdcf_colorsr = np.array([[165, 42, 42]])

    # Third Ventricle
    vent3_code = 14

    # Hypothalamus
    hypf_namesl  = ["hypo-lh-hypothalamus"]
    hypf_colorsl = np.array([[204, 182, 142]])

    hypf_namesr  = ["hypo-rh-hypothalamus"]
    hypf_colorsr = np.array([[204, 182, 142]])

    # 6. ---------------- Detection of hypothalamic nuclei (Iglesias) ------------ #
    hypi_codesl = np.array([801, 802, 803, 804, 805])
    hypi_codesr = np.array([806, 807, 808, 809, 810])

    idx = search(st_codes, hypi_codesl)
    hypi_namesl = np.array(st_names)[idx].tolist()
    hypi_colorsl = st_colors[idx]
    hypi_namesl = [s.replace('_', '-') for s in hypi_namesl]
    hypi_namesl = [s.replace('L-hypothalamus-', 'hypo-lh-').lower() for s in hypi_namesl]

    idx = search(st_codes, hypi_codesr)
    hypi_namesr = np.array(st_names)[idx].tolist()
    hypi_colorsr = st_colors[idx]
    hypi_namesr = [s.replace('_', '-') for s in hypi_namesr]
    hypi_namesr = [s.replace('R-hypothalamus-', 'hypo-rh-').lower() for s in hypi_namesr]

    # 7. ---------------- Detecting Cerbellum structures (Freesurfer) ------------------------- #
    cerebf_codesl = np.ones((1,), dtype=int) * 8
    idx = [st_codes.index(cerebf_codesl)]
    cerebf_namesl = ['cer-lh-cerebellum']
    cerebf_colorsl = st_colors[idx]

    cerebf_codesr = np.ones((1,), dtype=int) * 47
    idx = [st_codes.index(cerebf_codesr)]
    cerebf_namesr = ['cer-rh-cerebellum']
    cerebf_colorsr = st_colors[idx]

    # 8. ---------------- Detection of Brainstem (FreeSurfer) ------------ #
    bstemf_codes = np.ones((1,), dtype=int) * 16
    idx = search(st_codes, bstemf_codes)
    bstemf_names = ['brain-stem-brainstem']
    bstemf_colors = st_colors[idx]

    # 8. ---------------- Detection of Brainstem (Iglesias) ------------ #
    bstemi_codes = np.array([173, 174, 175, 178])
    idx = search(st_codes, bstemi_codes)
    bstemi_names = np.array(st_names)[idx].tolist()
    bstemi_colors = st_colors[idx]
    prefix = "brain-stem-"
    bstemi_names = [prefix + s.lower() for s in bstemi_names]

    ## ================ Creating the new parcellation
    lut_lines = ['{:<4} {:<40} {:>3} {:>3} {:>3} {:>3} \n \n'.format("#No.", "Label Name:", "R", "G", "B", "A")]
    parc_desc_lines = ["# Parcellation code: " + parccode]

    if len(parccode) != 9:  # Length of the parcellation string
        parccode = parccode.ljust(9, 'N')

    ##### ========== Selecting the cortical parcellation ============== #####
    try:
        atlas_str     = parcdict["Cortical"][parccode[0]]["String"]
        atlas_desc    = parcdict["Cortical"][parccode[0]]["Description"]
        atlas_cita    = parcdict["Cortical"][parccode[0]]["Citation"]
        atlas_type    = parcdict["Cortical"][parccode[0]]["Type"]
        atlas_names   = parcdict["Cortical"][parccode[0]]["Name"]
        atlas_surfloc = parcdict["Cortical"][parccode[0]]["OutSurfLocation"]
        atlas_volloc  = parcdict["Cortical"][parccode[0]]["OutVolLocation"]
        surfatlas_dir = os.path.join(deriv_dir, atlas_surfloc, path_cad, 'anat')
        volatlas_dir = os.path.join(deriv_dir, atlas_volloc, path_cad, 'anat')
        parc_desc_lines.append(atlas_desc + ' ' + atlas_cita)
    except:
        print("Incorrect cortical parcellation. The available parcellations are: ")
        _print_availab_parcels("Cortical")
        sys.exit(1)

    # Finding FreeSurfer folder
    fs_indivdir = os.path.join(fssubj_dir, fullid)

    if not os.path.isfile(os.path.join(fs_indivdir,'mri', 'aparc+aseg.mgz')):
        _launch_freesurfer(t1, fssubj_dir, fullid)

    if os.path.isfile(os.path.join(fs_indivdir,'mri', 'aparc+aseg.mgz')):
        # Finding the cortical parcellations
        out_sparc = glob(surfatlas_dir + os.path.sep + fullid + '*' + atlas_str + '*.annot')

        if len(out_sparc) != len(atlas_names)*2:
            print("The selected cortical parcellation (" +  parcdict["Cortical"][parccode[0]]["Atlas"] + ") is not computed.")
            print("Trying to compute the correponding cortical parcellation.")

            cwd = os.getcwd()
            atlas_dir = os.path.join(cwd, atlas_type.upper() + '_atlases')

            # Creating the link for fsaverage
            fsave_dir = os.path.join(fshome_dir,'subjects','fsaverage')
            if not os.path.isdir(os.path.join(fssubj_dir,'fsaverage')):
                process = subprocess.run(['ln', '-s', fsave_dir, fssubj_dir],
                                        stdout=subprocess.PIPE, universal_newlines=True)

            out_sparc = []
            for atlas in atlas_names:

                # Mapping the annot file to individual space
                # 1. Left Hemisphere
                ind_annot     = os.path.join(fssubj_dir, fullid, 'label', 'lh.' + atlas + '.annot') # Annot in individual space (freesurfer subject's directory)
                out_annot = os.path.join(surfatlas_dir, fullid + '_hemi-L_space-orig_' + atlas + '_dparc.annot')
                if not os.path.isfile(out_annot):
                    if atlas_type == 'annot':
                        fs_annot  = os.path.join(atlas_dir,
                                                'lh.' + atlas + '.annot')  # Annot in fsaverage space (Atlas folder)
                        out_annot = _launch_annot2ind(fs_annot, ind_annot, 'lh', surfatlas_dir, fullid, atlas)

                    elif atlas_type == 'gcs':
                        fs_gcs    = os.path.join(atlas_dir,
                                                'lh.' + atlas + '.gcs')  # GCS in fsaverage space (Atlas folder)
                        out_annot = _launch_gcs2ind(fssubj_dir, fs_gcs, ind_annot, 'lh', surfatlas_dir, fullid, atlas)

                    elif atlas_type == 'freesurfer':

                        # Copying the resulting annot to the output folder
                        if atlas == "aparc":
                            atlas_str = "atlas-desikan_desc-aparc"
                        elif atlas == "aparc.a2009s":
                            atlas_str = "atlas-destrieux_desc-a2009s"

                        out_annot = os.path.join(surfatlas_dir, fullid + '_hemi-L_space-orig_' + atlas_str + '_dparc.annot')
                        subprocess.run(['cp', ind_annot, out_annot], stdout=subprocess.PIPE, universal_newlines=True)


                out_sparc.append(out_annot) # Annot in individual space (Atlases subject's directory)

                # 2. Right Hemisphere
                ind_annot = os.path.join(fssubj_dir, fullid, 'label', 'rh.' + atlas + '.annot') # Annot in individual space (freesurfer subject's directory)
                out_annot = os.path.join(surfatlas_dir, fullid + '_hemi-R_space-orig_' + atlas + '_dparc.annot')
                if not os.path.isfile(out_annot):
                    if atlas_type == 'annot':
                        fs_annot  = os.path.join(atlas_dir,
                                                'rh.' + atlas + '.annot')  # Annot in fsaverage space (Atlas folder)
                        out_annot = _launch_annot2ind(fs_annot, ind_annot, 'rh', surfatlas_dir, fullid, atlas)

                    elif atlas_type == 'gcs':
                        fs_gcs    = os.path.join(atlas_dir,
                                                'rh.' + atlas + '.gcs')  # GCS in fsaverage space (Atlas folder)
                        out_annot = _launch_gcs2ind(fssubj_dir, fs_gcs, ind_annot, 'rh', surfatlas_dir, fullid, atlas)

                    elif atlas_type == 'freesurfer':
                        # Copying the resulting annot to the output folder
                        if atlas == "aparc":
                            atlas_str = "atlas-desikan_desc-aparc"
                        elif atlas == "aparc.a2009s":
                            atlas_str = "atlas-destrieux_desc-a2009s"

                        out_annot = os.path.join(surfatlas_dir, fullid + '_hemi-R_space-orig_' + atlas_str + '_dparc.annot')
                        subprocess.run(['cp', ind_annot, out_annot], stdout=subprocess.PIPE, universal_newlines=True)

                out_sparc.append(out_annot) # Annot in individual space (Atlases subject's directory)


        # Right hemisphere (Surface parcellation)
        rh_cparc = [s for s in out_sparc if "hemi-R" in s]  # Right cortical parcellation
        rh_cparc.sort()

        # Left hemisphere (Surface parcellation)
        lh_cparc = [s for s in out_sparc if "hemi-L" in s]  # Left cortical parcellation
        lh_cparc.sort()

        # Selecting the volumetric parcellation
        # vol_cparc = glob(volatlas_dir + os.path.sep + fullid + '*' + atlas_str + '*.nii.gz') # Cortical surface parcellation (.annot, .gii)
        # vol_cparc.sort()  # Volumetric cortical parcellation

        # vol2look = ['grow' + s + 'mm' for s in growwm]

        vol2look = []
        for s in growwm:
            if s.isnumeric():
                vol2look.append('grow' + s + 'mm')
            else:
                vol2look.append('grow' + s)


        vol_cparc = []
        for g in growwm:
            tmp = glob(
                volatlas_dir + os.path.sep + fullid + '*' + atlas_str + '*grow' + g +'*.nii.gz')  # Cortical surface parcellation (.annot, .gii)
            vol_cparc.extend(tmp)
        vol_cparc.sort()  # Volumetric cortical parcellation

        # If the volumetric parcellation does not exist it will try to create it
        if len(vol_cparc) != len(atlas_names)*len(growwm):
            vol_cparc = []
            for atlas in atlas_names:
                out_vol = _launch_surf2vol(fssubj_dir, volatlas_dir, fullid, atlas, growwm)

        if len(rh_cparc) != len(lh_cparc):  # Verifying the same number of parcellations for both hemispheres
            print(
                "Error: Some surface-based cortical parcellations are missing. Different number of files per hemisphere.\n")
            sys.exit(1)

        for f in rh_cparc:  # Verifying the existence of all surface-based cortical parcellations (Right)
            temp_file = Path(f)
            if not temp_file.is_file():
                print("Error: Some surface-based cortical parcellations are missing.\n")
                sys.exit(1)

        for f in lh_cparc:  # Verifying the existence of all surface-based cortical parcellations (Left)
            temp_file = Path(f)
            if not temp_file.is_file():
                print("Error: Some surface-based cortical parcellations are missing.\n")
                sys.exit(1)

        for f in vol_cparc:  # Verifying the existence of volumetric parcellations
            temp_file = Path(f)
            if not temp_file.is_file():
                print("Error: Volumetric parcellations are missing.\n")
                sys.exit(1)

        # Loading Aparc parcellation
        if 'F' in parccode:
            tempDir = os.path.join(deriv_dir, 'freesurfer', fullid)
            aparc_mgz = os.path.join(tempDir, 'mri', 'aseg.mgz')
            raw_mgz = os.path.join(tempDir, 'mri', 'rawavg.mgz')
            aseg_nii = os.path.join(tempDir, 'tmp', 'aseg.nii.gz')
            process = subprocess.run(['mri_vol2vol', '--mov', aparc_mgz, '--targ', raw_mgz, '--regheader', '--o', aseg_nii, '--no-save-reg', '--interp', 'nearest'],
                                    stdout=subprocess.PIPE, universal_newlines=True)
            temp_file = Path(aseg_nii)
            if temp_file.is_file():
                aseg_parc = nib.load(temp_file)
                aseg_parc = aseg_parc.get_fdata()
            else:
                print("Error: Cannot create the parcellation because there are missing files.\n")
                sys.exit(1)

        if 'aseg_parc' in locals():
            dim = aseg_parc.shape
        else:
            aseg_parc = nib.load(vol_cparc[0])
            dim = aseg_parc.shape

        outparc_lh = np.zeros((dim[0], dim[1], dim[2]), dtype='int16')  # Temporal parcellation for the left hemisphere
        outparc_rh = np.zeros((dim[0], dim[1], dim[2]), dtype='int16')  # Temporal parcellation for the right hemisphere

        # Loading FIRST parcellation
        if parccode[1] == 'R' or parccode[2] == 'R' or parccode[3] == 'R' or parccode[4] == 'R' or parccode[7] == 'R':
            atlas_str = parcdict["Subcortical"]["R"]["String"]
            atlas_desc = parcdict["Subcortical"]["R"]["Description"]
            atlas_volloc = parcdict["Subcortical"]["R"]["OutVolLocation"]
            volatlas_dir = os.path.join(deriv_dir, atlas_volloc, path_cad, 'anat')

            first_parc = os.path.join(volatlas_dir, fullid + '_space-orig_atlas-' + atlas_str + '_dseg.nii.gz')
            if not os.path.isfile(first_parc):

                # Creating ouput directory
                if not os.path.isdir(volatlas_dir):
                    try:
                        os.makedirs(volatlas_dir)
                    except OSError:
                        print("Failed to make nested output directory")

                # Running FIRST subcortical parcellation
                process = subprocess.run(
                    ['run_first_all', '-i', t1, '-o', volatlas_dir + os.path.sep + 'temp'],
                    stdout=subprocess.PIPE, universal_newlines=True)

                # Changing name
                process = subprocess.run(
                    ['mv', 'temp_all_fast_firstseg.nii.gz', first_parc],
                    stdout=subprocess.PIPE, universal_newlines=True)

                # Deleting temporary files
                process = subprocess.run(
                    ['rm', '-rf', 'temp*'],
                    stdout=subprocess.PIPE, universal_newlines=True)

            first_parc = nib.load(first_parc)
            first_parc = aseg_parc.get_fdata()

        ##### ========== Selecting Subcortical parcellation ============== #####
        try:
            atlas_str  = parcdict["Subcortical"][parccode[1]]["String"]
            atlas_desc = parcdict["Subcortical"][parccode[1]]["Description"]
            atlas_cita    = parcdict["Subcortical"][parccode[1]]["Citation"]
            atlas_volloc = parcdict["Subcortical"][parccode[1]]["OutVolLocation"]
            volatlas_dir = os.path.join(deriv_dir, atlas_volloc, path_cad, 'anat')
            parc_desc_lines.append(atlas_desc + ' ' + atlas_cita)
        except:
            _print_availab_parcels("Subcortical")
            sys.exit(1)

        if parccode[1] == 'F':

            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(aseg_parc, subc_codesr, outparc_rh, 0)

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(aseg_parc, subc_codesl, outparc_lh, 0)

        elif parccode[1] == 'R':  # TODO
            # Volumetric subcortical parcellation
            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(first_parc, subc_codesr, outparc_rh, 0)

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(first_parc, subc_codesl, outparc_lh, 0)


        ##### ========== Selecting Thalamic parcellation ============== #####
        try:
            atlas_str = parcdict["Thalamus"][parccode[2]]["String"]
            atlas_desc = parcdict["Thalamus"][parccode[2]]["Description"]
            atlas_cita = parcdict["Thalamus"][parccode[2]]["Citation"]
            atlas_volloc = parcdict["Thalamus"][parccode[2]]["OutVolLocation"]
            volatlas_dir = os.path.join(deriv_dir, atlas_volloc, path_cad, 'anat')
            parc_desc_lines.append(atlas_desc + ' ' + atlas_cita)

        except:
            _print_availab_parcels("Thalamus")
            sys.exit(1)

        if parccode[2] == 'F':
            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(aseg_parc, thalf_codesr, outparc_rh, st_lengtrh)

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(aseg_parc, thalf_codesl, outparc_lh, st_lengtlh)

        elif parccode[2] == 'R':  # TODO
            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(first_parc, thalf_codesr, outparc_rh, st_lengtrh)

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(first_parc, thalf_codesl, outparc_lh, st_lengtlh)

        elif parccode[2] == 'I':

            vol_tparc = glob(os.path.join(volatlas_dir, fullid + '_space-orig_desc-' + atlas_str + '_dseg.nii.gz'))

            if not vol_tparc:
                vol_tparc = os.path.join(volatlas_dir, fullid + '_space-orig_desc-' + atlas_str + '_dseg.nii.gz')
                vol_tparc = _fs_addon_parcellations(vol_tparc, fullid, fssubj_dir, 'thalamus', atlas_str)

            # Reading the thalamic parcellation
            temp_iparc = nib.load(vol_tparc[0])
            temp_iparc = temp_iparc.get_fdata()

            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(temp_iparc, thali_codesr, outparc_rh, st_lengtrh)

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(temp_iparc, thali_codesl, outparc_lh, st_lengtlh)

        elif parccode[2] == 'M':

            # Thalamic parcellation based on Najdenovska et al, 2018
            vol_tparc = glob(os.path.join(volatlas_dir, fullid + '_space-orig_desc-' + atlas_str + '_dseg.nii.gz'))

            if not vol_tparc:
                vol_tparc = os.path.join(volatlas_dir, fullid + '_space-orig_desc-' + atlas_str + '_dseg.nii.gz')

                # Computing thalamic nuclei using atlas-based parcellation
                vol_tparc = _compute_abased_thal_parc(t1, vol_tparc, deriv_dir, path_cad, fullid, aseg_nii, atlas_str)

            # Reading the thalamic parcellation
            temp_iparc = nib.load(vol_tparc[0])
            temp_iparc = temp_iparc.get_fdata()

            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(temp_iparc, thalm_codesr, outparc_rh, st_lengtrh)

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(temp_iparc, thalm_codesl, outparc_lh, st_lengtlh)

        ##### ========== Selecting Amygdala parcellation ============== #####
        try:
            atlas_str = parcdict["Amygdala"][parccode[3]]["String"]
            atlas_desc = parcdict["Amygdala"][parccode[3]]["Description"]
            atlas_cita = parcdict["Amygdala"][parccode[3]]["Citation"]
            atlas_volloc = parcdict["Amygdala"][parccode[3]]["OutVolLocation"]
            volatlas_dir = os.path.join(deriv_dir, atlas_volloc, path_cad, 'anat')
            parc_desc_lines.append(atlas_desc + ' ' + atlas_cita)

        except:
            _print_availab_parcels("Amygdala")
            sys.exit(1)

        if parccode[3] == 'F':
            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(aseg_parc, amygf_codesr, outparc_rh, st_lengtrh)

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(aseg_parc, amygf_codesl, outparc_lh, st_lengtlh)

        elif parccode[3] == 'I':

            # # Volumetric amygdala parcellation (vol_aparc)
            vol_aparc_lh = glob(os.path.join(volatlas_dir, fullid + '*L*' + atlas_str + '*.nii.gz'))
            vol_aparc_rh = glob(os.path.join(volatlas_dir, fullid + '*R*' + atlas_str + '*.nii.gz'))

            if not vol_aparc_lh or not vol_aparc_rh:

                vol_aparc_lh = os.path.join(volatlas_dir, fullid + '_space-orig_hemi-L_desc-' + atlas_str + '_dseg.nii.gz')
                vol_aparc_rh = os.path.join(volatlas_dir, fullid + '_space-orig_hemi-R_desc-' + atlas_str + '_dseg.nii.gz')

                # Computing amygdala nuclei
                vol_tparc = _fs_addon_parcellations(vol_aparc_lh, fullid, fssubj_dir, 'amygdala', atlas_str)
                vol_aparc_lh = [vol_tparc[0]]
                vol_aparc_rh = [vol_tparc[1]]

            # Reading the amygdala parcellation
            temp_iparc = nib.load(vol_aparc_rh[0])
            temp_iparc = temp_iparc.get_fdata()

            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(temp_iparc, amygi_codesr, outparc_rh, st_lengtrh)

            # Reading the amygdala parcellation
            temp_iparc = nib.load(vol_aparc_lh[0])
            temp_iparc = temp_iparc.get_fdata()

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(temp_iparc, amygi_codesl, outparc_lh, st_lengtlh)

        elif parccode[3] == 'R':
            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(first_parc, amygf_codesr, outparc_rh, st_lengtrh)

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(first_parc, amygf_codesl, outparc_lh, st_lengtlh)


        ##### ========== Selecting Hippocampus parcellation ============== #####
        try:
            atlas_str = parcdict["Hippocampus"][parccode[4]]["String"]
            atlas_desc = parcdict["Hippocampus"][parccode[4]]["Description"]
            atlas_cita = parcdict["Hippocampus"][parccode[4]]["Citation"]
            atlas_volloc = parcdict["Hippocampus"][parccode[4]]["OutVolLocation"]
            volatlas_dir = os.path.join(deriv_dir, atlas_volloc, path_cad, 'anat')
            parc_desc_lines.append(atlas_desc + ' ' + atlas_cita)

            if parccode[4] == 'H':
                atlas_str_ig = parcdict["Hippocampus"]["I"]["String"]
                atlas_volloc_ig = parcdict["Hippocampus"]["I"]["OutVolLocation"]
                volatlas_dir_ig = os.path.join(deriv_dir, atlas_volloc_ig, path_cad, 'anat')

        except:
            _print_availab_parcels("Hippocampus")
            sys.exit(1)

        if parccode[4] == 'F':

            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(aseg_parc, hippf_codesr, outparc_rh, st_lengtrh)

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(aseg_parc, hippf_codesl, outparc_lh, st_lengtlh)

        elif parccode[4] == 'I':

            # Hippocampus parcellation based on Iglesias et al, 2015

            vol_hparc_lh = glob(os.path.join(volatlas_dir, fullid + '*L*' + atlas_str + '*.nii.gz'))
            vol_hparc_rh = glob(os.path.join(volatlas_dir, fullid + '*R*' + atlas_str + '*.nii.gz'))

            if not vol_hparc_lh or not vol_hparc_rh:
                vol_hparc_lh = os.path.join(volatlas_dir, fullid + '_space-orig_hemi-L_desc-' + atlas_str + '_dseg.nii.gz')
                vol_hparc_rh = os.path.join(volatlas_dir, fullid + '_space-orig_hemi-R_desc-' + atlas_str + '_dseg.nii.gz')

                # Computing thalamic nuclei using atlas-based parcellation
                vol_tparc = _fs_addon_parcellations(vol_hparc_lh, fullid, fssubj_dir, 'hippocampus', atlas_str)
                vol_hparc_lh = [vol_tparc[0]]
                vol_hparc_rh = [vol_tparc[1]]

            # Reading the hippocampus parcellation
            temp_iparc = nib.load(vol_hparc_rh[0])
            temp_iparc = temp_iparc.get_fdata()

            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(temp_iparc, hippi_codesr, outparc_rh, st_lengtrh)

            # Reading the hippocampus parcellation
            temp_iparc = nib.load(vol_hparc_lh[0])
            temp_iparc = temp_iparc.get_fdata()

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(temp_iparc, hippi_codesl, outparc_lh, st_lengtlh)

        elif parccode[4] == 'H':

            # Detecting the parcellation of hippocampal subfields
            atlas_str_ig = parcdict["Hippocampus"]["I"]["String"]

            # Hippocampus parcellation based on Iglesias et al, 2015
            vol_hparc_lh = glob(volatlas_dir + os.path.sep + fullid + '*L*' + atlas_str_ig + '*.nii.gz')
            vol_hparc_rh = glob(volatlas_dir + os.path.sep + fullid + '*R*' + atlas_str_ig + '*.nii.gz')

            if not vol_hparc_lh or not vol_hparc_rh:
                vol_hparc_lh = os.path.join(volatlas_dir, fullid, '_space-orig_hemi-L_desc-' + atlas_str_ig + '_dseg.nii.gz')
                vol_hparc_rh = os.path.join(volatlas_dir, fullid, '_space-orig_hemi-R_desc-' + atlas_str_ig + '_dseg.nii.gz')

                # Computing thalamic nuclei using atlas-based parcellation
                vol_tparc = _fs_addon_parcellations(vol_hparc_lh, fullid, fssubj_dir, 'hippocampus', atlas_str)
                vol_hparc_lh = [vol_tparc[0]]
                vol_hparc_rh = [vol_tparc[1]]

            # Reading the hippocampus parcellation
            temp_iparc = nib.load(vol_hparc_rh[0])
            temp_iparc = temp_iparc.get_fdata()
            temp_iparc = _subfields2hbt(temp_iparc, hippi_codesr)

            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(temp_iparc, hipph_codesr, outparc_rh, st_lengtrh)

            # Reading the hippocampus parcellation
            temp_iparc = nib.load(vol_hparc_lh[0])
            temp_iparc = temp_iparc.get_fdata()
            temp_iparc = _subfields2hbt(temp_iparc, hippi_codesl)

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(temp_iparc, hipph_codesl, outparc_lh, st_lengtlh)

        elif parccode[4] == 'R':
            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(aseg_parc, hippf_codesr, outparc_rh, st_lengtrh)

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(aseg_parc, hippf_codesl, outparc_lh, st_lengtlh)

        ##### ========== Selecting Hypothalamus parcellation ============== #####
        try:
            atlas_str = parcdict["Hypothalamus"][parccode[5]]["String"]
            atlas_desc = parcdict["Hypothalamus"][parccode[5]]["Description"]
            atlas_cita = parcdict["Hypothalamus"][parccode[5]]["Citation"]
            atlas_volloc = parcdict["Hypothalamus"][parccode[5]]["OutVolLocation"]
            volatlas_dir = os.path.join(deriv_dir, atlas_volloc, path_cad, 'anat')
            parc_desc_lines.append(atlas_desc + ' ' + atlas_cita)

        except:
            _print_availab_parcels("Hypothalamus")
            sys.exit(1)

        if parccode[5] == 'F':
            print("This needs to be implemented")

            im_tmp           = np.zeros(aseg_parc.shape)
            ind_vent         = np.where(aseg_parc == vent3_code)
            im_tmp[ind_vent] = 1
            im_dil = mask_dilation(im_tmp, 5) # Dilating 5mm
            thirdV = op.abspath('{}.nii.gz'.format("ventricle3"))
            hdr = V.get_header()
            hdr2 = hdr.copy()
            hdr2.set_data_dtype(np.int16)
            iflogger.info("    ... Image saved to {}".format(thirdV))
            img = ni.Nifti1Image(tmp, V.get_affine(), hdr2)
            ni.save(img, thirdV)

            iflogger.info("  > Dilate the ventricule image")
            thirdV_dil = op.abspath('{}_dil.nii.gz'.format("ventricle3"))
            fslmaths_cmd = 'fslmaths {} -kernel sphere 5 -dilD {}'.format(thirdV, thirdV_dil)
            iflogger.info("    ... Command: {}".format(fslmaths_cmd))
            process = subprocess.Popen(fslmaths_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            proc_stdout = process.communicate()[0].strip()

            tmp = ni.load(thirdV_dil).get_data()
            indrhypothal = np.where((im_dil == 1) & (aseg_parc == vdcf_codel))
            indlhypothal = np.where((im_dil == 1) & (aseg_parc == vdcf_coder))
            del (tmp)

        elif parccode[5] == 'I':

            # Thalamic parcellation based on Iglesias et al, 2019
            vol_yparc = glob(os.path.join(volatlas_dir, fullid + '*' + atlas_str + '*.nii.gz'))

            if not vol_yparc:
                vol_yparc = os.path.join(volatlas_dir, fullid + '_space-orig_desc-' + atlas_str + '_dseg.nii.gz')
                vol_yparc = _fs_addon_parcellations(vol_yparc, fullid, fssubj_dir, 'hypothalamus', atlas_str)

            # Reading the hypothalamus parcellation
            temp_iparc = nib.load(vol_yparc[0])
            temp_iparc = temp_iparc.get_fdata()

            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(temp_iparc, hypi_codesr, outparc_rh, st_lengtrh)

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(temp_iparc, hypi_codesl, outparc_lh, st_lengtlh)

        ##### ========== Selecting Cerebellum parcellation ============== #####
        try:
            atlas_str = parcdict["Cerebellum"][parccode[6]]["String"]
            atlas_desc = parcdict["Cerebellum"][parccode[6]]["Description"]
            atlas_cita = parcdict["Cerebellum"][parccode[6]]["Citation"]
            atlas_volloc = parcdict["Cerebellum"][parccode[6]]["OutVolLocation"]
            volatlas_dir = os.path.join(deriv_dir, atlas_volloc, path_cad, 'anat')
            parc_desc_lines.append(atlas_desc + ' ' + atlas_cita)

        except:
            _print_availab_parcels("Cerebellum")
            sys.exit(1)

        if parccode[6] == 'F':
            # Right Hemisphere
            outparc_rh, st_lengtrh = _search_in_atlas(aseg_parc, cerebf_codesr, outparc_rh, st_lengtrh)

            # Left Hemisphere
            outparc_lh, st_lengtlh = _search_in_atlas(aseg_parc, cerebf_codesl, outparc_lh, st_lengtlh)

        ##### ========== Selecting Brainstem parcellation ============== #####
        try:
            atlas_str = parcdict["Brainstem"][parccode[7]]["String"]
            atlas_desc = parcdict["Brainstem"][parccode[7]]["Description"]
            atlas_cita = parcdict["Brainstem"][parccode[7]]["Citation"]
            atlas_volloc = parcdict["Brainstem"][parccode[7]]["OutVolLocation"]
            volatlas_dir = os.path.join(deriv_dir, atlas_volloc, path_cad, 'anat')
            parc_desc_lines.append(atlas_desc + ' ' + atlas_cita)

        except:
            _print_availab_parcels("Brainstem")
            sys.exit(1)

        if parccode[7] == 'F':
            # Adding Brainstem parcellation to the left parcellation
            outparc_lh, st_lengtlh = _search_in_atlas(temp_iparc, bstemf_codes, outparc_lh, st_lengtlh)

        elif parccode[7] == 'I':
            # Brainstem parcellation based on Iglesias et al, 2018
            vol_bparc = glob(os.path.join(volatlas_dir, fullid + '*' + atlas_str + '*.nii.gz'))

            if not vol_bparc:
                vol_bparc = os.path.join(volatlas_dir, fullid + '_space-orig_desc-' + atlas_str + '_dseg.nii.gz')
                vol_bparc = _fs_addon_parcellations(vol_bparc, fullid, fssubj_dir, 'brainstem', atlas_str)

            # Reading the Brainstem parcellation
            temp_iparc = nib.load(vol_bparc[0])
            temp_iparc = temp_iparc.get_fdata()

            # Adding Brainstem parcellation to the left parcellation
            outparc_lh, st_lengtlh = _search_in_atlas(temp_iparc, bstemi_codes, outparc_lh, st_lengtlh)

        ##### ========== Selecting white matter parcellation ============== #####
        try:
            atlas_str = parcdict["GyralWM"][parccode[8]]["String"]
            atlas_desc = parcdict["GyralWM"][parccode[8]]["Description"]
            atlas_cita = parcdict["Cortical"][parccode[8]]["Citation"]
            atlas_volloc = parcdict["GyralWM"][parccode[8]]["OutVolLocation"]
            volatlas_dir = os.path.join(deriv_dir, atlas_volloc, path_cad, 'anat')
            parc_desc_lines.append(atlas_desc + ' ' + atlas_cita)

        except:
            _print_availab_parcels("GyralWM")
            sys.exit(1)

        if parccode[8] == 'F':
            sss = 1


        # Removing temporal aseg image
        os.remove(aseg_nii)

        parc_desc_lines.append("\n")

        # Creating ouput directory
        out_dir = os.path.join(deriv_dir, 'chimera-atlases', path_cad, 'anat')
        if not os.path.isdir(out_dir):
            try:
                os.makedirs(out_dir)
            except OSError:
                print("Failed to make nested output directory")


        # Loop around each parcellation
        for i in np.arange(0, len(rh_cparc)):
            right_sdata = nib.freesurfer.io.read_annot(rh_cparc[i], orig_ids=False)
            rh_codes = right_sdata[0]
            rh_colors = right_sdata[1][1:, 0:3]
            rh_stnames = right_sdata[2][1:]

            left_sdata = nib.freesurfer.io.read_annot(lh_cparc[i], orig_ids=False)
            lh_codes = left_sdata[0]
            lh_colors = left_sdata[1][1:, 0:3]
            lh_stnames = left_sdata[2][1:]

            fname = os.path.basename(rh_cparc[i])

            temp = fname.split('_')
            scaleid = [s for s in temp if "desc-" in s]  # Detect if the label key exist

            # Selecting the volumetric parcellations for all the wm grow levels
            if scaleid:
                grow_parcs = [s for s in vol_cparc if scaleid[0] in s]
            else:
                grow_parcs = vol_cparc

            nctx_rh = len(rh_stnames)  # Number of cortical regions in the right hemisphere
            nctx_lh = len(lh_stnames)  # Number of cortical regions in the left hemisphere
            nroi_right = nctx_rh + st_lengtrh  # Number of regions in the right hemisphere

            rh_luttable = ["# Right Hemisphere. Cortical Structures"]
            lh_luttable = ["# Left Hemisphere. Cortical Structures"]

            ##### ========== LUT Cortical Surface (Right Hemisphere)============== #####
            # rh_scode, rh_ctab,
            for roi_pos, roi_name in enumerate(rh_stnames):
                temp_name = 'ctx-rh-{}'.format(roi_name.decode("utf-8"))
                rh_luttable.append(
                    '{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(roi_pos + 1, temp_name, rh_colors[roi_pos, 0],
                                                                rh_colors[roi_pos, 1], rh_colors[roi_pos, 2],
                                                                0))
            maxlab_rh = roi_pos + 1

            for roi_pos, roi_name in enumerate(lh_stnames):
                temp_name = 'ctx-lh-{}'.format(roi_name.decode("utf-8"))
                lh_luttable.append(
                    '{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(nroi_right + roi_pos + 1, temp_name,
                                                                lh_colors[roi_pos, 0],
                                                                lh_colors[roi_pos, 1], lh_colors[roi_pos, 2], 0))
            maxlab_lh = nroi_right + roi_pos + 1

            rh_luttable.append('\n')
            lh_luttable.append('\n')

            ##### ========== Selecting Subcortical parcellation ============== #####
            if parccode[1] == 'F' or parccode[1] == 'R':
                rh_luttable.append("# Right Hemisphere. Subcortical Structures")
                for roi_pos, roi_name in enumerate(subc_namesr):
                    rh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_rh + roi_pos + 1, roi_name,
                                                                                    subc_colorsr[roi_pos, 0],
                                                                                    subc_colorsr[roi_pos, 1],
                                                                                    subc_colorsr[roi_pos, 2], 0))
                maxlab_rh = maxlab_rh + roi_pos + 1

                lh_luttable.append("# Left Hemisphere. Subcortical Structures")
                for roi_pos, roi_name in enumerate(subc_namesl):
                    lh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_lh + roi_pos + 1, roi_name,
                                                                                    subc_colorsl[roi_pos, 0],
                                                                                    subc_colorsl[roi_pos, 1],
                                                                                    subc_colorsl[roi_pos, 2], 0))
                maxlab_lh = maxlab_lh + roi_pos + 1


            elif parccode[1] == 'R':  # TODO
                t = 1

            rh_luttable.append('\n')
            lh_luttable.append('\n')

            ##### ========== Selecting Thalamic parcellation ============== #####
            if parccode[2] == 'F' or parccode[2] == 'R':

                rh_luttable.append("# Right Hemisphere. Thalamic Structures")
                for roi_pos, roi_name in enumerate(thalf_namesr):
                    rh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_rh + roi_pos + 1, roi_name,
                                                                                    thalf_colorsr[roi_pos, 0],
                                                                                    thalf_colorsr[roi_pos, 1],
                                                                                    thalf_colorsr[roi_pos, 2], 0))
                maxlab_rh = maxlab_rh + roi_pos + 1

                lh_luttable.append("# Left Hemisphere. Thalamic Structures")
                for roi_pos, roi_name in enumerate(thalf_namesl):
                    lh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_lh + roi_pos + 1, roi_name,
                                                                                    thalf_colorsl[roi_pos, 0],
                                                                                    thalf_colorsl[roi_pos, 1],
                                                                                    thalf_colorsl[roi_pos, 2], 0))
                maxlab_lh = maxlab_lh + roi_pos + 1

            elif parccode[2] == 'R':  # TODO

                # Volumetric thalamic parcellation
                volatlas_dir = os.path.join(deriv_dir, 'fsl-subcparc', subjId, sesId,
                                    'anat')  # Subcortical parcellation based on Patenaude et al, 2011
                parc_desc_lines.append(
                    "# 3. Thalamic parcellation (R): FIRST thalamic parcellation")

            elif parccode[2] == 'I':

                rh_luttable.append("# Right Hemisphere. Thalamic Structures")
                for roi_pos, roi_name in enumerate(thali_namesr):
                    rh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_rh + roi_pos + 1, roi_name,
                                                                                    thali_colorsr[roi_pos, 0],
                                                                                    thali_colorsr[roi_pos, 1],
                                                                                    thali_colorsr[roi_pos, 2], 0))
                maxlab_rh = maxlab_rh + roi_pos + 1

                lh_luttable.append("# Left Hemisphere. Thalamic Structures")
                for roi_pos, roi_name in enumerate(thali_namesl):
                    lh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_lh + roi_pos + 1, roi_name,
                                                                                    thali_colorsl[roi_pos, 0],
                                                                                    thali_colorsl[roi_pos, 1],
                                                                                    thali_colorsl[roi_pos, 2], 0))
                maxlab_lh = maxlab_lh + roi_pos + 1

            elif parccode[2] == 'M':

                rh_luttable.append("# Right Hemisphere. Thalamic Structures")
                for roi_pos, roi_name in enumerate(thalm_namesr):
                    rh_luttable.append('{:<4} {:<50} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_rh + roi_pos + 1, roi_name,
                                                                                    thalm_colorsr[roi_pos, 0],
                                                                                    thalm_colorsr[roi_pos, 1],
                                                                                    thalm_colorsr[roi_pos, 2], 0))
                maxlab_rh = maxlab_rh + roi_pos + 1

                lh_luttable.append("# Left Hemisphere. Thalamic Structures")
                for roi_pos, roi_name in enumerate(thalm_namesl):
                    lh_luttable.append('{:<4} {:<50} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_lh + roi_pos + 1, roi_name,
                                                                                    thalm_colorsr[roi_pos, 0],
                                                                                    thalm_colorsr[roi_pos, 1],
                                                                                    thalm_colorsr[roi_pos, 2], 0))
                maxlab_lh = maxlab_lh + roi_pos + 1

            rh_luttable.append('\n')
            lh_luttable.append('\n')

            ##### ========== Selecting Amygdala parcellation ============== #####
            if parccode[3] == 'F' or parccode[3] == 'R':

                rh_luttable.append("# Right Hemisphere. Amygdala")
                for roi_pos, roi_name in enumerate(amygf_namesr):
                    rh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_rh + roi_pos + 1, roi_name,
                                                                                    amygf_colorsr[roi_pos, 0],
                                                                                    amygf_colorsr[roi_pos, 1],
                                                                                    amygf_colorsr[roi_pos, 2], 0))
                maxlab_rh = maxlab_rh + roi_pos + 1

                lh_luttable.append("# Left Hemisphere. Amygdala")
                for roi_pos, roi_name in enumerate(amygf_namesl):
                    lh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_lh + roi_pos + 1, roi_name,
                                                                                    amygf_colorsl[roi_pos, 0],
                                                                                    amygf_colorsl[roi_pos, 1],
                                                                                    amygf_colorsl[roi_pos, 2], 0))
                maxlab_lh = maxlab_lh + roi_pos + 1


            elif parccode[3] == 'I':

                rh_luttable.append("# Right Hemisphere. Amygdala Nuclei")
                for roi_pos, roi_name in enumerate(amygi_namesr):
                    rh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_rh + roi_pos + 1, roi_name,
                                                                                    amygi_colorsr[roi_pos, 0],
                                                                                    amygi_colorsr[roi_pos, 1],
                                                                                    amygi_colorsr[roi_pos, 2], 0))
                maxlab_rh = maxlab_rh + roi_pos + 1

                lh_luttable.append("# Left Hemisphere. Amygdala Nuclei")
                for roi_pos, roi_name in enumerate(amygi_namesl):
                    lh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_lh + roi_pos + 1, roi_name,
                                                                                    amygi_colorsl[roi_pos, 0],
                                                                                    amygi_colorsl[roi_pos, 1],
                                                                                    amygi_colorsl[roi_pos, 2], 0))
                maxlab_lh = maxlab_lh + roi_pos + 1

            rh_luttable.append('\n')
            lh_luttable.append('\n')

            ##### ========== Selecting Hippocampus parcellation ============== #####
            if parccode[4] == 'F' or parccode[4] == 'R':

                rh_luttable.append("# Right Hemisphere. Hippocampus")
                for roi_pos, roi_name in enumerate(hippf_namesr):
                    rh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_rh + roi_pos + 1, roi_name,
                                                                                    hippf_colorsr[roi_pos, 0],
                                                                                    hippf_colorsr[roi_pos, 1],
                                                                                    hippf_colorsr[roi_pos, 2], 0))
                maxlab_rh = maxlab_rh + roi_pos + 1

                lh_luttable.append("# Left Hemisphere. Hippocampus")
                for roi_pos, roi_name in enumerate(hippf_namesl):
                    lh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_lh + roi_pos + 1, roi_name,
                                                                                    hippf_colorsl[roi_pos, 0],
                                                                                    hippf_colorsl[roi_pos, 1],
                                                                                    hippf_colorsl[roi_pos, 2], 0))
                maxlab_lh = maxlab_lh + roi_pos + 1

            elif parccode[4] == 'I':

                rh_luttable.append("# Right Hemisphere. Hippocampus subfields")
                for roi_pos, roi_name in enumerate(hippi_namesr):
                    rh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_rh + roi_pos + 1, roi_name,
                                                                                    hippi_colorsr[roi_pos, 0],
                                                                                    hippi_colorsr[roi_pos, 1],
                                                                                    hippi_colorsr[roi_pos, 2], 0))
                maxlab_rh = maxlab_rh + roi_pos + 1

                lh_luttable.append("# Left Hemisphere. Hippocampus subfields")
                for roi_pos, roi_name in enumerate(hippi_namesl):
                    lh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_lh + roi_pos + 1, roi_name,
                                                                                    hippi_colorsl[roi_pos, 0],
                                                                                    hippi_colorsl[roi_pos, 1],
                                                                                    hippi_colorsl[roi_pos, 2], 0))
                maxlab_lh = maxlab_lh + roi_pos + 1

            elif parccode[4] == 'H':

                rh_luttable.append(
                    "# Right Hemisphere. Hippocampus subfields grouped in Head, Body, Tail and Fissure")
                for roi_pos, roi_name in enumerate(hipph_namesr):
                    rh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_rh + roi_pos + 1, roi_name,
                                                                                    hipph_colorsr[roi_pos, 0],
                                                                                    hipph_colorsr[roi_pos, 1],
                                                                                    hipph_colorsr[roi_pos, 2], 0))
                maxlab_rh = maxlab_rh + roi_pos + 1

                lh_luttable.append(
                    "# Left Hemisphere. Hippocampus subfields grouped in Head, Body, Tail and Fissure")
                for roi_pos, roi_name in enumerate(hipph_namesl):
                    lh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_lh + roi_pos + 1, roi_name,
                                                                                    hipph_colorsl[roi_pos, 0],
                                                                                    hipph_colorsl[roi_pos, 1],
                                                                                    hipph_colorsl[roi_pos, 2], 0))
                maxlab_lh = maxlab_lh + roi_pos + 1

            rh_luttable.append('\n')
            lh_luttable.append('\n')
            ##### ========== Selecting Hypothalamus parcellation ============== #####
            if parccode[5] == 'F':

                # # Volumetric Hypothalamus parcellation (vol_yparc)
                volatlas_dir = os.path.join(deriv_dir, 'freesurfer-volparc', subjId, sesId, 'anat')
                vol_yparc = glob(volatlas_dir + os.path.sep + '*desikan*.nii.gz')
                vol_yparc = vol_yparc[0]

            elif parccode[5] == 'I':
                rh_luttable.append("# Right Hemisphere. Hypothalamus nuclei")
                for roi_pos, roi_name in enumerate(hypi_namesr):
                    rh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_rh + roi_pos + 1, roi_name,
                                                                                    hypi_colorsr[roi_pos, 0],
                                                                                    hypi_colorsr[roi_pos, 1],
                                                                                    hypi_colorsr[roi_pos, 2], 0))
                maxlab_rh = maxlab_rh + roi_pos + 1

                lh_luttable.append("# Left Hemisphere. Hypothalamus nuclei")
                for roi_pos, roi_name in enumerate(hypi_namesl):
                    lh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_lh + roi_pos + 1, roi_name,
                                                                                    hypi_colorsl[roi_pos, 0],
                                                                                    hypi_colorsl[roi_pos, 1],
                                                                                    hypi_colorsl[roi_pos, 2], 0))
                maxlab_lh = maxlab_lh + roi_pos + 1

            rh_luttable.append('\n')
            lh_luttable.append('\n')

            ##### ========== Selecting Cerebellum parcellation ============== #####
            if parccode[6] == 'F':
                rh_luttable.append("# Right Hemisphere. Cerebellum")
                for roi_pos, roi_name in enumerate(cerebf_namesr):
                    rh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_rh + roi_pos + 1, roi_name,
                                                                                    cerebf_colorsr[roi_pos, 0],
                                                                                    cerebf_colorsr[roi_pos, 1],
                                                                                    cerebf_colorsr[roi_pos, 2], 0))
                maxlab_rh = maxlab_rh + roi_pos + 1

                lh_luttable.append("# Left Hemisphere. Cerebellum")
                for roi_pos, roi_name in enumerate(cerebf_namesl):
                    lh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_lh + roi_pos + 1, roi_name,
                                                                                    cerebf_colorsl[roi_pos, 0],
                                                                                    cerebf_colorsl[roi_pos, 1],
                                                                                    cerebf_colorsl[roi_pos, 2], 0))
                maxlab_lh = maxlab_lh + roi_pos + 1

            rh_luttable.append('\n')
            lh_luttable.append('\n')

            ##### ========== Selecting Brainstem parcellation ============== #####
            lh_luttable.append("# Left Hemisphere. BrainStem parcellation")
            if parccode[7] == 'F' or parccode[7] == 'R':
                for roi_pos, roi_name in enumerate(bstemf_names):
                    lh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_lh + roi_pos + 1, roi_name,
                                                                                    bstemf_colors[roi_pos, 0],
                                                                                    bstemf_colors[roi_pos, 1],
                                                                                    bstemf_colors[roi_pos, 2], 0))
                maxlab_lh = maxlab_lh + roi_pos + 1

            elif parccode[7] == 'I':
                for roi_pos, roi_name in enumerate(bstemi_names):
                    lh_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(maxlab_lh + roi_pos + 1, roi_name,
                                                                                    bstemi_colors[roi_pos, 0],
                                                                                    bstemi_colors[roi_pos, 1],
                                                                                    bstemi_colors[roi_pos, 2], 0))
                maxlab_lh = maxlab_lh + roi_pos + 1

            lh_luttable.append('\n')
            ##### ========== Selecting White matter parcellation ============== #####
            if parccode[8] == 'F':

                # # Volumetric White matter parcellation (vol_wparc)
                wm_luttable = ["# Global White Matter"]
                wm_luttable.append(
                    '{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(3000, wm_names[0], wm_colors[0, 0], wm_colors[0, 1],
                                                                wm_colors[0, 2], 0))
                wm_luttable.append('\n')
                wm_luttable.append("# Right Hemisphere. Gyral White Matter Structures")

                # rh_scode, rh_ctab,
                for roi_pos, roi_name in enumerate(rh_stnames):
                    temp_name = 'wm-rh-{}'.format(roi_name.decode("utf-8"))
                    wm_luttable.append(
                        '{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(roi_pos + 3001, temp_name,
                                                                    255 - rh_colors[roi_pos, 0],
                                                                    255 - rh_colors[roi_pos, 1],
                                                                    255 - rh_colors[roi_pos, 2],
                                                                    0))
                wm_luttable.append('\n')

                wm_luttable.append("# Left Hemisphere. Gyral White Matter Structures")
                for roi_pos, roi_name in enumerate(lh_stnames):
                    temp_name = 'wm-lh-{}'.format(roi_name.decode("utf-8"))
                    wm_luttable.append('{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format(nroi_right + roi_pos + 3001, temp_name,
                                                                                    255 - lh_colors[roi_pos, 0],
                                                                                    255 - lh_colors[roi_pos, 1],
                                                                                    255 - lh_colors[roi_pos, 2], 0))

            elif parccode[8] == 'N':
                vol_wparc = []
                parc_desc_lines.append("# 9. White matter parcellation (N): No WM segmentation.")

            for gparc in grow_parcs:

                # Reading the cortical parcellation
                temp_iparc = nib.load(gparc)
                affine = temp_iparc.affine
                temp_iparc = temp_iparc.get_fdata()

                out_atlas = np.zeros(np.shape(temp_iparc), dtype='int16')

                # Adding cortical regions (Right Hemisphere)

                ind = np.where(np.logical_and(temp_iparc > 2000, temp_iparc < 3000))
                out_atlas[ind[0], ind[1], ind[2]] = temp_iparc[ind[0], ind[1], ind[2]] - np.ones((len(ind[0]),)) * 2000

                # Adding the rest of the regions (Right Hemisphere)
                ind = np.where(outparc_rh > 0)
                out_atlas[ind[0], ind[1], ind[2]] = outparc_rh[ind[0], ind[1], ind[2]] + np.ones((len(ind[0]),)) * nctx_rh

                # Adding cortical regions (Left Hemisphere)
                ind = np.where(np.logical_and(temp_iparc > 1000, temp_iparc < 2000))
                out_atlas[ind[0], ind[1], ind[2]] = temp_iparc[ind[0], ind[1], ind[2]] - np.ones(
                    (len(ind[0]),)) * 1000 + np.ones(
                    (len(ind[0]),)) * nroi_right

                # Adding the rest of the regions (Left Hemisphere + Brainstem)
                ind = np.where(outparc_lh > 0)
                out_atlas[ind[0], ind[1], ind[2]] = outparc_lh[ind[0], ind[1], ind[2]] + np.ones(
                    (len(ind[0]),)) * nctx_lh + np.ones((len(ind[0]),)) * nroi_right

                # Adding global white matter
                bool_ind = np.in1d(temp_iparc, wm_codes)
                bool_ind = np.reshape(bool_ind, np.shape(temp_iparc))
                ind      = np.where(bool_ind)
                out_atlas[ind[0], ind[1], ind[2]] = 3000

                # Adding right white matter
                ind      = np.where(np.logical_and(temp_iparc > 4000, temp_iparc < 5000))

                if parccode[8] == 'F':
                    out_atlas[ind[0], ind[1], ind[2]] = temp_iparc[ind[0], ind[1], ind[2]] - np.ones((len(ind[0]),)) * 1000
                else:
                    out_atlas[ind[0], ind[1], ind[2]] = 3000

                    # Adding left white matter
                ind = np.where(np.logical_and(temp_iparc > 3000, temp_iparc < 4000))

                if parccode[8] == 'F':
                    out_atlas[ind[0], ind[1], ind[2]] = temp_iparc[ind[0], ind[1], ind[2]] + np.ones((len(ind[0]),)) * nroi_right
                else:
                    out_atlas[ind[0], ind[1], ind[2]] = 3000

                # Creating output filename using pybids
                ######## ------------- Creating the full ID. It is used for a correct image file naming. ------------ #
                # if 'session' in ent_dict.keys():
                #     pattern = "sub-{subject}_ses-{session}_run-{run}_space-{space}[_atlas-{atlas}][_desc-{desc}]_{suffix}.nii.gz"
                #     patternlut = "sub-{subject}_ses-{session}_run-{run}_space-{space}[_atlas-{atlas}][_desc-{desc}]_{suffix}.lut"
                #     patterntsv = "sub-{subject}_ses-{session}_run-{run}_space-{space}[_atlas-{atlas}][_desc-{desc}]_{suffix}.tsv"
                # else:
                #     pattern = "sub-{subject}_run-{run}_space-{space}[_atlas-{atlas}][_desc-{desc}]_{suffix}.nii.gz"
                #     patternlut = "sub-{subject}_run-{run}_space-{space}[_atlas-{atlas}][_desc-{desc}]_{suffix}.lut"
                #     patterntsv = "sub-{subject}_run-{run}_space-{space}[_atlas-{atlas}][_desc-{desc}]_{suffix}.tsv"

                fname            = os.path.basename(gparc)
                templist         = fname.split('_')
                tempVar          = [s for s in templist if "desc-" in s]  # Left cortical parcellation
                descid           = tempVar[0].split('-')[1]

                # laydict          = layout.parse_file_entities(gparc)
                # laydict["atlas"] = 'chimera' + parccode
                # laydict["desc"]  = descid

                base_id = fullid.split('_')
                base_id.append('space-orig')
                base_id.append('atlas-' + 'chimera' + parccode)
                base_id.append('desc-' + descid)

                # Saving the parcellation
                outparcFilename = os.path.join(out_dir, '_'.join(base_id) + '_dseg.nii.gz')
                imgcoll          = nib.Nifti1Image(out_atlas.astype('int16') , affine)
                nib.save(imgcoll, outparcFilename)

                # Saving the colorLUT
                colorlutFilename = os.path.join(out_dir, '_'.join(base_id) + '_dseg.lut')

                now              = datetime.now()
                date_time        = now.strftime("%m/%d/%Y, %H:%M:%S")
                time_lines       = ['# $Id: {} {} \n'.format(colorlutFilename, date_time),
                                    '# Corresponding parcellation: ',
                                    '# ' + outparcFilename + '\n']

                hdr_lines        = ['{:<4} {:<40} {:>3} {:>3} {:>3} {:>3}'.format("#No.", "Label Name:", "R", "G", "B", "A")]
                lut_lines        = time_lines + parc_desc_lines + hdr_lines + rh_luttable + lh_luttable + wm_luttable
                with open(colorlutFilename, 'w') as colorLUT_f:
                    colorLUT_f.write('\n'.join(lut_lines))

                st_codes_lut, st_names_lut, st_colors_lut = read_fscolorlut(colorlutFilename)

                # Saving the TSV
                tsvFilename = os.path.join(out_dir, '_'.join(base_id) + '_dseg.tsv')
                _parc_tsv_table(st_codes_lut, st_names_lut, st_colors_lut, tsvFilename)



def main():
    # 0. Handle inputs
    parser = _build_args_parser()
    args = parser.parse_args()

    print(args)
    if args.verbose is not None:
        v = np.int(args.verbose[0])
    else:
        v = 0
        print('- Verbose set to 0\n')
    if v:
        print('\nInputs\n')
    #

    # Getting the path of the current running python file
    cwd          = os.getcwd()
    bids_dirercts   = args.bidsdir[0].split(sep=',')
    # deriv_dir    = args.derivdir[0]
    deriv_dirercts   = args.derivdir[0].split(sep=',')
    parcodes     = args.parcodes[0].split(sep=',')
    t1s2run_file = args.t1s[0]
    growwm       = args.growwm[0]
    nthreads     = int(args.nthreads[0])
    growwm       = growwm.split(',')
    
    # Avoiding the layout to decrease indexing time
    suffix = "T1w"
    extension = 'nii'
    cad2look = '*_{}.{}*'.format(suffix, extension)
    t1s = []

    for cont, bids_dir in enumerate(bids_dirercts):
        if len(bids_dirercts) == len(deriv_dirercts):
            deriv_dir = deriv_dirercts[cont]
        else:
            deriv_dir = deriv_dirercts[0]

        for path, subdir, files in os.walk(bids_dir):
            if not deriv_dir in path:
                for file in glob(os.path.join(path, cad2look)):
                    t1s.append(file)


        if os.path.isfile(t1s2run_file):
            t1s = _select_t1s(t1s, t1s2run_file)
        else:
            t12run = t1s2run_file.split(',')
            t1s = [s for s in t1s if any(xs in s for xs in t12run)]

        # parcodes = ['KFMFIIFIF', 'SFMFIIFIF', 'CFMFIIFIF']
        # t1s = t1s[0:93]
        Nparc = len(parcodes)
        nsubj = len(t1s)
        
        for p, parccode in enumerate(parcodes):
            failed = []
            print("Parcellation: % d"% (p+1), "of % d"% (Nparc))
            if nthreads == 1:
                for i, t1 in enumerate(t1s):
                    # ent_dict = layout.parse_file_entities(t1)

                    t1_name = os.path.basename(t1)
                    temp = t1_name.split("_")
                    fullid = '_'.join(temp[:-1])
                    _printprogressbar(i + 1, nsubj,
                                    'Processing T1w --> ' + fullid + ': ' + '(' + str(i + 1) + '/' + str(nsubj) + ')')
                    # _build_parcellation(layout, bids_dir, deriv_dir, ent_dict, parccode)
                    _build_parcellation(t1, bids_dir, deriv_dir, parccode, growwm)
            else:
                start_time = time.perf_counter()
                ndwis = len(t1s)
                ncores = os.cpu_count()

                if nthreads > 4:
                    nthreads = nthreads - 4
                with concurrent.futures.ProcessPoolExecutor(ncores) as executor:
                #     results = [executor.submit(do_something, sec) for sec in secs]
                    results = list(executor.map(_build_parcellation, t1s,
                    [bids_dir] * ndwis, [deriv_dir] * ndwis, [parccode] * ndwis, [growwm] * ndwis))

                end_time = time.perf_counter()

                print(f'Finished in {end_time - start_time} seconds (s)...')


if __name__ == "__main__":
    main()
