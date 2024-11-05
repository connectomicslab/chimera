#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from os import PathLike
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import bids
from bids import BIDSLayout
import nibabel as nib
from glob import glob
from operator import itemgetter
from datetime import datetime
import argparse
from typing import Union
import shutil
import uuid
import csv
import json
import subprocess
import scipy.ndimage as sc
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait
from concurrent.futures import as_completed
import time
import copy
from threading import Lock

from templateflow import api as tflow


import clabtoolkit.misctools as cltmisc
import clabtoolkit.freesurfertools as cltfree
import clabtoolkit.parcellationtools as cltparc
import clabtoolkit.bidstools as cltbids
import clabtoolkit.segmentationtools as cltseg
from rich.progress import Progress

class bcolors:
    """
    Class to define the colors for the terminal output.
    
    
    """
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    OKYELLOW = '\033[93m'
    OKRED = '\033[91m'
    OKMAGENTA = '\033[95m'
    PURPLE = '\033[35m'
    OKCYAN = '\033[96m'
    DARKCYAN = "\033[36m"
    ORANGE = "\033[48:5:208m%s\033[m"
    OKWHITE = '\033[97m'
    DARKWHITE = '\033[37m'
    OKBLACK = '\033[30m'
    OKGRAY = '\033[90m'
    OKPURPLE = '\033[35m'
    
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    ITALIC = '\033[3m'
    UNDERLINE = '\033[4m'

# Define the Chimera class. This class will be used to create and work with Chimera objects
class Chimera:
    """
    Class to create and work with Chimera objects.
    
    Parameters
    ----------
    parc_code : str
        Parcellation code.
    parc_dict_file : str
        Parcellation dictionary file.
    supra_folder : str
        Folder containing the supraregions TSV files.
        
    Returns
    -------
    Chimera object
    
    
    """

    def __init__(self, parc_code,
                        scale: Union[str, list]  = None,
                        seg: Union[str, list]  = None,
                        parc_dict_file: str = None, 
                        supra_folder: str = None):
        """
        Initialize the Chimera class
        
        Parameters
        ----------
        parc_code : str
            Parcellation code.
        parc_dict_file : str
            Parcellation dictionary file.
        supra_folder : str
            Folder containing the supraregions TSV files.
            
        Returns
        -------
        Chimera object

        """

        chim_dir = os.path.dirname(os.path.abspath(__file__))

        # Rise an error if the parcellation code is not provided
        if parc_code is None:
            raise ValueError("Please provide a parcellation code")
        
        if parc_dict_file is not None:
            if not os.path.isfile(parc_dict_file):
                raise ValueError("The parcellation dictionary file does not exist")
        else:
            parc_dict_file = os.path.join(chim_dir, 'config', 'supraregions_dictionary.json')
        
        if supra_folder is not None:
            if not os.path.isdir(supra_folder):
                raise ValueError("The the folder containing the supra-regions TSV files is not valid")
            else:
                self.suprafolder = supra_folder
        else:
            self.suprafolder = os.path.join(chim_dir, 'config', 'supraregions')
        
        self.parc_dict, self.supra_dict = _load_parcellations_info(parc_json=parc_dict_file, supra_folder=supra_folder)
        
        ####  Filtering the parcellation dictionary according to the parcellation code ####
        supra_names = list(self.parc_dict.keys())
        
        temp_dict = {}
        for i in range(len(parc_code)):
            if parc_code[i] in self.parc_dict[supra_names[i]].keys():
                
                # Defining the dictionary for the method
                meth_dict = {}
                meth_dict["code"] = parc_code[i] # Add the parcellation code
                
                # Append the information of the parcellation using that method
                meth_dict.update(self.parc_dict[supra_names[i]][parc_code[i]])

                # Filtering the parcellation names by the scale and segmentation
                if i == 0:
                    parcel_names = meth_dict['parcels']
                    
                    # Filtering the parcellation names by the scale 
                    if scale is not None:
                        if isinstance(scale , list):
                            # If the scale is a list do a loop over the elements and
                            # verify if the scale contains the string '_scale-'
                            scale_tmp = []
                            for sc in scale :
                                if '_scale-' not in sc:
                                    scale_tmp.append('_scale-' + sc)
                            
                        elif isinstance(scale , str):
                            if '_scale-' not in scale :
                                scale_tmp = '_scale-' + scale 

                        parcel_names = cltmisc._filter_by_substring(parcel_names, scale_tmp, boolcase=False)

                    # Filtering the parcellation names by the segmentation
                    if seg is not None:
                        if isinstance(seg, list):
                            # If the seg is a list do a loop over the elements and
                            # verify if the seg contains the string '_seg-'
                            seg_tmp = []
                            for sc in seg:
                                if '_seg-' not in sc:
                                    seg_tmp.append('_seg-' + sc)
                            
                        elif isinstance(seg, str):
                            if '_seg-' not in seg:
                                seg_tmp = '_seg-' + seg

                        parcel_names = cltmisc._filter_by_substring(parcel_names, seg_tmp, boolcase=False)
                    
                    # Saving the new parcels names
                    meth_dict['parcels'] = parcel_names
            
                # Adding the dictionary to the temp_dict
                temp_dict[supra_names[i]] = meth_dict
                
            else:
                # Print a message that the parcellation code is not present in the dictionary
                print("The parcellation code {} is not present in the dictionary for the supra-region {}.".format(parc_code[i],supra_names[i]))
                # Error message and exit
                #sys.exit(1)
    
        self.parc_dict = temp_dict
        self.parc_code = parc_code
        self.scale = scale
        self.seg = seg
    
    
    
    def _prepare_templates(self, fssubj_dir:str = None):
        """
        This method prepares the templates for the Chimera parcellation.
        Based on the code of the parcellation, it will download the necessary templates
        from the TemplateFlow repository or set up the templates directory for each supra region.
        
        Parameters
        ----------
        fssubj_dir : str
            FreeSurfer directory.
        
        """
        
        global pipe_dict
        
        # Setting up the FreeSurfer directory
        cltfree.FreeSurferSubject._set_freesurfer_directory(fssubj_dir)
        
        # Create the simlink to the FreeSurfer directory
        if "Cortical" in self.parc_dict.keys():
            cltfree._create_fsaverage_links(fssubj_dir, fsavg_dir=None, refsubj_name=self.parc_dict["Cortical"]["reference"])
        
        
        # Detecting the base directory
        chim_dir = os.path.dirname(os.path.abspath(__file__))

        # Reading the names of the supra-regions
        supra_names = list(self.parc_dict.keys())
        
        for supra in supra_names:

            atlas_src     = self.parc_dict[supra]["source"]
            atlas_str     = self.parc_dict[supra]["atlas"]
            atlas_ref     = self.parc_dict[supra]["reference"]
            
            if supra == 'Cortical':
                
                # Atributes for the cortical parcellation
                atlas_type    = self.parc_dict[supra]["type"]
                atlas_names   = self.parc_dict[supra]["parcels"]
                
                # Selecting the source and downloading the parcellation
                if atlas_src == 'templateflow':
                    atlas_ext = '.gii'
                    method = 'annot2indiv'
                    
                    tflow_home = _set_templateflow_home(pipe_dict["packages"]["templateflow"]["home_dir"])
                    ctx_parc_lh = tflow.get(template=atlas_ref, atlas=atlas_str, hemi='L' , suffix='dseg', extension='.label.gii')
                    ctx_parc_rh = tflow.get(template=atlas_ref, atlas=atlas_str, hemi='R' , suffix='dseg', extension='.label.gii')
                    
                    # Convert the list of PosfixPath to strings
                    if isinstance(ctx_parc_lh, list):
                        ctx_parc_lh = [str(x) for x in ctx_parc_lh]
                        
                    elif isinstance(ctx_parc_lh, PathLike):
                        ctx_parc_lh = [str(ctx_parc_lh)]
                    
                    if isinstance(ctx_parc_rh, list):
                        ctx_parc_rh = [str(x) for x in ctx_parc_rh]
                        
                    elif isinstance(ctx_parc_rh, PathLike):
                        ctx_parc_rh = [str(ctx_parc_rh)]

                    # Select the files that contain the atlas names
                    ctx_parc_lh = cltmisc._filter_by_substring(ctx_parc_lh, atlas_names, boolcase=False)
                    ctx_parc_rh = cltmisc._filter_by_substring(ctx_parc_rh, atlas_names, boolcase=False)
                    ctx_parc_lh.sort()
                    ctx_parc_rh.sort()
                    
                    ctx_parc_lh_annot = []
                    ctx_parc_rh_annot = []
                    for i, parc_file in enumerate(ctx_parc_lh):
                        
                        # Detect which element in atlas_names is in the string ctx_parc_lh
                        at_name = [s for s in atlas_names if s in ctx_parc_lh[i]]
                        
                        if at_name:
                            
                            # Moving the gifti to native space
                            tmp_annot = os.path.join(fssubj_dir, atlas_ref,'label','lh.' + at_name[0] + '.annot')
                            tmp_refsurf = os.path.join(fssubj_dir, atlas_ref,'surf','lh.white')
                            ctx_parc_lh_annot.append(tmp_annot)
                            lh_obj = cltfree.AnnotParcellation.gii2annot(gii_file = parc_file,
                                                                        ref_surf = tmp_refsurf,
                                                                        annot_file = tmp_annot, 
                                                                        cont_tech=pipe_dict["packages"]["freesurfer"]["cont_tech"], 
                                                                        cont_image=pipe_dict["packages"]["freesurfer"]["container"])
                            
                            tmp_annot   = os.path.join(fssubj_dir, atlas_ref,'label','rh.' + at_name[0] + '.annot')
                            tmp_refsurf = os.path.join(fssubj_dir, atlas_ref,'surf','rh.white')
                            ctx_parc_rh_annot.append(tmp_annot)
                            rh_obj = cltfree.AnnotParcellation.gii2annot(gii_file = ctx_parc_rh[i],
                                                                        ref_surf = tmp_refsurf,
                                                                        annot_file = tmp_annot, 
                                                                        cont_tech=pipe_dict["packages"]["freesurfer"]["cont_tech"], 
                                                                        cont_image=pipe_dict["packages"]["freesurfer"]["container"])
                    
                    if not ctx_parc_lh_annot or not ctx_parc_rh_annot:
                        raise ValueError("Cortical parcellations should be supplied for both hemispheres.")
                    else: 
                        meth_dict = {'method': method,'reference': atlas_ref,
                                                                'labels': {'lh': ctx_parc_lh_annot, 'rh': ctx_parc_rh_annot}}

                elif atlas_src == 'local':

                    if atlas_type == 'annot':
                        atlas_dir = os.path.join(chim_dir, 'data','annot_atlases')
                        atlas_ext = '.annot'
                        method = 'annot2indiv'
                        
                    elif atlas_type == 'gcs':
                        atlas_dir = os.path.join(chim_dir, 'data', 'gcs_atlases')
                        atlas_ext = '.gcs'
                        method = 'gcs2indiv'
                        

                    ctx_parc_lh = glob(os.path.join(atlas_dir, '*-L_*' + atlas_ext))
                    ctx_parc_rh = glob(os.path.join(atlas_dir, '*-R_*' + atlas_ext))
                    
                    # Filtering for selecting the correct cortical parcellation
                    ctx_parc_lh = cltmisc._filter_by_substring(ctx_parc_lh, atlas_names, boolcase=False)
                    ctx_parc_rh = cltmisc._filter_by_substring(ctx_parc_rh, atlas_names, boolcase=False)
                    ctx_parc_lh.sort()
                    ctx_parc_rh.sort()
                    
                    if not ctx_parc_lh or not ctx_parc_rh:
                        raise ValueError("Cortical parcellations should be supplied for both hemispheres.")
                    
                    else:
                    
                        meth_dict = {'method': method,'reference': atlas_ref,
                                            'labels': {'lh': ctx_parc_lh, 'rh': ctx_parc_rh}}
                    
                    
            else:
                atlas_cad = self.parc_dict[supra]['atlas']
                type_cad = self.parc_dict[supra]['type']
                atlas_ref = self.parc_dict[supra]["reference"] 
                
                if atlas_src == 'templateflow':
                    
                    # Getting the templates
                    # Reference space
                    tflow_home = _set_templateflow_home(pipe_dict["packages"]["templateflow"]["home_dir"])
                    t1_temp = tflow.get(atlas_ref, desc=None, resolution=1, suffix='T1w', extension='nii.gz')
                    
                    # Getting the thalamic nuclei spams 
                    if type_cad == 'spam':
                        atlas_file = tflow.get(atlas_ref, desc=None, resolution=1,atlas=atlas_cad, suffix='probseg', extension='nii.gz')
                        
                    elif type_cad == 'maxprob':
                        atlas_file = tflow.get(atlas_ref, desc=None, resolution=1,atlas=atlas_cad, suffix='dseg', extension='nii.gz')
                        
                    meth_dict = {'method': 'atlasbased','type':type_cad,'reference': str(t1_temp),
                                            'labels': str(atlas_file)}
                
                elif atlas_src == 'local':
                    atlas_dir = os.path.join(chim_dir, 'data', 'vol_atlases')
                    
                    t1_temp = glob(os.path.join(atlas_dir, '*' + atlas_cad + '*_T1w.nii.gz'))
                    
                    if type_cad == 'spam':
                        atlas_file = glob(os.path.join(atlas_dir, '*' + atlas_cad + '*_probseg.nii.gz'))
                        
                    elif type_cad == 'maxprob':
                        atlas_file = glob(os.path.join(atlas_dir, '*' + atlas_cad + '*_dseg.nii.gz'))
                
                    meth_dict = {'method': 'atlasbased','type':type_cad,'reference': str(t1_temp),
                                            'labels': str(atlas_file)}
                    
                elif atlas_src == 'freesurfer':
                        
                    meth_dict = {'method': 'comform2native','type':None,'reference': 'native',
                                            'labels': None}
                
                elif atlas_src == 'freesurferextra':
                        meth_dict = {'method': 'comform2native','type':None,'reference': 'native',
                                                'labels': atlas_src.lower() }
                else:
                    meth_dict = {'method': None, 'type':None,'reference': 'native', 'labels': None}
            
            self.parc_dict[supra]["processing"] = meth_dict
    
    def _create_table(self, wm_index_offset:int = 3000, 
                    reg2rem:Union[list, str] = ['unknown', 'medialwall', 'corpuscallosum']):
        """
        This method creates a table with the regions that will be created using the Chimera parcellation.
        It allows to verify the region distribution for a specified parcelation code.
        
        Parameters
        ----------
        wm_index_offset : int
            Index offset for the white matter parcellation (default = 3000).
        
        reg2rem : list
            List of regions to remove from the parcellation.

        """
        
        # Detecting the base directory
        cwd = os.path.dirname(os.path.abspath(__file__))
        chim_dir = os.path.dirname(cwd)

        
        # Reading the names of the supra-regions
        supra_names = list(self.parc_dict.keys())

        lh_noctx_codes  = []
        rh_noctx_codes  = []
        lh_noctx_names  = []
        rh_noctx_names  = []
        lh_noctx_colors = []
        rh_noctx_colors = []
        bs_noctx_codes  = []
        bs_noctx_names  = []
        bs_noctx_colors = []

        ctx_parc_lh = []
        ctx_parc_rh = []
        
        desc_noctx = []
        for supra in supra_names:
            
            # Getting the information of the common atributes
            atlas_code    = self.parc_dict[supra]["code"]
            atlas_str     = self.parc_dict[supra]["atlas"]
            atlas_desc    = self.parc_dict[supra]["description"]
            atlas_cita    = self.parc_dict[supra]["citation"]
            atlas_src     = self.parc_dict[supra]["source"]
            atlas_ref     = self.parc_dict[supra]["reference"]
            
            if supra == 'Cortical':
                
                # Atributes for the cortical parcellation
                atlas_type    = self.parc_dict[supra]["type"]
                atlas_names   = self.parc_dict[supra]["parcels"]
                
                # Selecting the source and downloading the parcellation
                if atlas_src == 'templateflow':
                    tflow_home = _set_templateflow_home(pipe_dict["packages"]["templateflow"]["home_dir"])
                    ctx_parc_lh = tflow.get(template=atlas_ref, atlas=atlas_str, hemi='L' , suffix='dseg', extension='.label.gii')
                    ctx_parc_rh = tflow.get(template=atlas_ref, atlas=atlas_str, hemi='R' , suffix='dseg', extension='.label.gii')
                    
                    # Convert the list of PosfixPath to strings
                    if isinstance(ctx_parc_lh, list):
                        ctx_parc_lh = [str(x) for x in ctx_parc_lh]
                        
                    elif isinstance(ctx_parc_lh, PathLike):
                        ctx_parc_lh = [str(ctx_parc_lh)]
                    
                    if isinstance(ctx_parc_rh, list):
                        ctx_parc_rh = [str(x) for x in ctx_parc_rh]
                        
                    elif isinstance(ctx_parc_rh, PathLike):
                        ctx_parc_rh = [str(ctx_parc_rh)]

                    # Select the files that contain the atlas names
                    ctx_parc_lh = cltmisc._filter_by_substring(ctx_parc_lh, atlas_names, boolcase=False)
                    ctx_parc_rh = cltmisc._filter_by_substring(ctx_parc_rh, atlas_names, boolcase=False)
                    
                elif atlas_src == 'local':

                    if atlas_type == 'annot':
                        atlas_dir = os.path.join(chim_dir, 'data','annot_atlases')
                        atlas_ext = '.annot'

                    elif atlas_type == 'gcs':
                        atlas_dir = os.path.join(chim_dir, 'data', 'gcs_atlases')
                        atlas_ext = '.gcs'
                        
                    ctx_parc_lh = glob(os.path.join(atlas_dir, '*-L_*' + atlas_ext))
                    ctx_parc_rh = glob(os.path.join(atlas_dir, '*-R_*' + atlas_ext))
                    
                    # Filtering for selecting the correct cortical parcellation
                    ctx_parc_lh = cltmisc._filter_by_substring(ctx_parc_lh, atlas_names, boolcase=False)
                    ctx_parc_rh = cltmisc._filter_by_substring(ctx_parc_rh, atlas_names, boolcase=False)
                    
            else:
                
                desc_noctx.append(atlas_desc)
                # Selecting the source and downloading the parcellation
                if atlas_src == 'templateflow':

                    # Reference space
                    tflow_home = _set_templateflow_home(pipe_dict["packages"]["templateflow"]["home_dir"])
                    ref_img = tflow.get(atlas_ref, desc=None, resolution=1, suffix='T1w', extension='nii.gz')
                    
                    # Getting the thalamic nuclei spams 
                    parc_img = tflow.get(atlas_ref, desc=None, resolution=1,atlas=atlas_str, suffix=atlas_type, extension='nii.gz')
                    
                if supra in self.supra_dict.keys():
                    meth_dict = self.parc_dict[supra]
                    st_dict = self.supra_dict[supra][supra][meth_dict["code"]]
                    if len(st_dict) == 1:
                            bs_noctx_codes = bs_noctx_codes + st_dict['mid']['index']
                            bs_noctx_names = bs_noctx_names + st_dict['mid']['name']
                            bs_noctx_colors = bs_noctx_colors + st_dict['mid']['color']
                            
                    elif len(st_dict) == 2:
                        lh_noctx_codes = lh_noctx_codes + st_dict['lh']['index']
                        rh_noctx_codes = rh_noctx_codes + st_dict['rh']['index']
                        
                        lh_noctx_names = lh_noctx_names + st_dict['lh']['name']
                        rh_noctx_names = rh_noctx_names + st_dict['rh']['name']
                        
                        lh_noctx_colors = lh_noctx_colors + st_dict['lh']['color']
                        rh_noctx_colors = rh_noctx_colors + st_dict['rh']['color']

        if rh_noctx_names:
            indexes = cltmisc._get_indexes_by_substring(rh_noctx_names, reg2rem).tolist()
            # Remove the elements in all_names and all_colors
            if indexes:
                for i in indexes:
                    rh_noctx_names.pop(i)
                    rh_noctx_codes.pop(i)
                    rh_noctx_colors.pop(i)

        if lh_noctx_names:
            indexes = cltmisc._get_indexes_by_substring(lh_noctx_names, reg2rem).tolist()
            # Remove the elements in all_names and all_colors
            if indexes:
                for i in indexes:
                    lh_noctx_names.pop(i)
                    lh_noctx_codes.pop(i)
                    lh_noctx_colors.pop(i)
        
        if bs_noctx_names:
            indexes = cltmisc._get_indexes_by_substring(bs_noctx_names, reg2rem).tolist()
            # Remove the elements in all_names and all_colors
            if indexes:
                for i in indexes:
                    bs_noctx_names.pop(i)
                    bs_noctx_codes.pop(i)
                    bs_noctx_codes.pop(i)
        

        # Creating the list of dataframes for the different parcellations
        tab_list = []
        desc_list = []
        parc_id = 'atlas-chimera' + self.parc_code
        parc_id_list = []
        
        # If ctx_parc_lh is empty, it means that the parcellation is not available
        if len(ctx_parc_lh) == 0:
            all_names  = rh_noctx_names  + lh_noctx_names  + bs_noctx_names
            all_colors = rh_noctx_colors + lh_noctx_colors + bs_noctx_colors
            index      = np.arange(1, len(all_names)+1).tolist()
            tab_df     = pd.DataFrame({'index': index, 'name': all_names, 'color': all_colors})
            tab_list.append(tab_df)
            
            gen_desc = ["# Parcellation code: " + self.parc_code ]
            gen_desc.append(desc_noctx)
            
            parc_id_list.append(parc_id)
            
        else:
            for i, parc_file in enumerate(ctx_parc_lh):
                
                gen_desc = ["# Parcellation code: " + self.parc_code ]
                
                tmp_name = os.path.basename(ctx_parc_lh[i])
                tmp_ent = tmp_name.split('_')[:-1]
                
                # Get the element that contains the string 'scale' and extract it
                scale_ent = [s for s in tmp_ent if 'scale' in s]
                if scale_ent:
                    scale_ent = scale_ent[0]
                    scale_ent = scale_ent.split('-')[1]
                    parc_id = parc_id + '_scale-' + scale_ent
                
                    # Add a the segmentation to to the string of the general description list
                    gen_desc[0] = gen_desc[0] + ". Scale: " + scale_ent
                
                # Get the element that contains the string 'seg' and extract it 
                seg_ent = [s for s in tmp_ent if 'seg' in s]
                
                if seg_ent:
                    seg_ent = seg_ent[0]
                    seg_ent = seg_ent.split('-')[1]
                    parc_id = parc_id + '_seg-' + seg_ent
                    
                    # Add a the segmentation to to the string of the general description list
                    gen_desc[0] = gen_desc[0] + ". Segmentation: " + seg_ent
                    
                gen_desc.append(self.parc_dict["Cortical"]["description"])
                gen_desc = gen_desc + desc_noctx
                    
                # Reading the cortical parcellations
                lh_obj = cltfree.AnnotParcellation(parc_file = ctx_parc_lh[i])
                rh_obj = cltfree.AnnotParcellation(parc_file = ctx_parc_rh[i])

                df_lh, out_tsv = lh_obj._export_to_tsv(prefix2add='ctx-lh-')
                df_rh, out_tsv = rh_obj._export_to_tsv(prefix2add='ctx-rh-')

                # Convert the column name of the dataframe to a list
                lh_ctx_name = df_lh['name'].tolist()
                rh_ctx_name = df_rh['name'].tolist()

                # Convert the column color of the dataframe to a list
                lh_ctx_color = df_lh['color'].tolist()
                rh_ctx_color = df_rh['color'].tolist()
                
                ## Removing elements from the table according to their name for both
                indexes = cltmisc._get_indexes_by_substring(lh_ctx_name, reg2rem).tolist()
                if indexes:
                    for i in indexes:
                        lh_ctx_name.pop(i)
                        lh_ctx_color.pop(i)
                
                indexes = cltmisc._get_indexes_by_substring(rh_ctx_name, reg2rem).tolist()
                if indexes:
                    for i in indexes:
                        rh_ctx_name.pop(i)
                        rh_ctx_color.pop(i)

                # Concatenating the lists
                if 'GyralWM' in self.parc_dict.keys():
                    gen_desc.append(self.parc_dict['GyralWM']["description"])
                    
                    wm_rh_name = cltmisc._correct_names(rh_ctx_name, replace=['ctx-rh-','wm-rh-'])
                    wm_rh_indexes = np.arange(1, len(wm_rh_name)+1) + wm_index_offset
                    wm_rh_indexes = wm_rh_indexes.tolist()
                    
                    wm_lh_name = cltmisc._correct_names(lh_ctx_name, replace=['ctx-lh-','wm-lh-'])
                    wm_lh_indexes = np.arange(1, len(wm_lh_name)+1) + len(rh_ctx_name) + len(rh_noctx_names) + wm_index_offset
                    wm_lh_indexes = wm_lh_indexes.tolist()
                    
                    wm_rh_color = rh_ctx_color
                    wm_lh_color = lh_ctx_color

                else:
                    wm_lh_indexes = []
                    wm_rh_indexes = []
                    wm_lh_name = []
                    wm_rh_name = []
                    wm_lh_color = []
                    wm_rh_color = []

                # Right hemisphere
                rh_all_names =  rh_ctx_name + rh_noctx_names
                rh_all_indexes = np.arange(1, len(rh_all_names)+1).tolist()
                
                # Left hemisphere
                lh_all_names =  lh_ctx_name +  lh_noctx_names  + bs_noctx_names + wm_rh_name + wm_lh_name
                lh_all_indexes = np.arange(1, len(lh_ctx_name +  lh_noctx_names  + bs_noctx_names)+1) + np.max(rh_all_indexes)
                lh_all_indexes = lh_all_indexes.tolist()
                
                rh_all_colors =  rh_ctx_color +  rh_noctx_colors
                lh_all_colors =  lh_ctx_color +  lh_noctx_colors + bs_noctx_colors + wm_rh_color + wm_lh_color

                # Concatenating the hemispheres
                all_names = rh_all_names + lh_all_names
                all_colors = rh_all_colors + lh_all_colors
                all_indexes = rh_all_indexes + lh_all_indexes + wm_rh_indexes + wm_lh_indexes

                # Generating a dataframe
                tab_df = pd.DataFrame({'index': all_indexes, 'name': all_names, 'color': all_colors})
                tab_list.append(tab_df)
                desc_list.append(gen_desc)
                parc_id_list.append(parc_id)

        # Add the tab_list as an attribute of the class
        self.regtable = {"parc_id": parc_id_list, "desc": desc_list, "table": tab_list}
        
    def _export_table(self, out_basename:str = None, format:Union[list, str] = 'tsv'):
        """
        This method exports the table of the regions to a TSV or a LUT file.
        
        Parameters
        ----------
        out_basename : str
            Output basename for the TSV or LUT file.
        
        format : str or list
            Format of the output file. It can be 'tsv', 'lut' or ['tsv, lut'].
        
        """
        
        if out_basename is None:
            # Exit if the output basename is not provided
            raise ValueError("Please provide an output basename for the TSV or LUT file.")
        
        out_name = os.path.basename(out_basename)
        
        out_dir = os.path.dirname(out_basename)
        out_dir = Path(out_dir)
        
        # Create the output directory if it does not exist
        out_dir.mkdir(parents=True, exist_ok=True)
        
        parc_ids = self.regtable["parc_id"]
        parc_desc = self.regtable["desc"]
        parc_tables = self.regtable["table"]
        
        # Export the table to a TSV file
        for i, tab_df in enumerate(parc_tables):
            
            if isinstance(format, list):
                if 'tsv' in format:
                    out_file_tsv = os.path.join(str(out_dir), out_name + '_' + parc_ids[i] + '.tsv')
                    cltparc.Parcellation.write_tsvtable(tsv_df=tab_df,out_file=out_file_tsv, force=True)

                if 'lut' in format:
                    out_file_lut = os.path.join(str(out_dir), out_name + '_' + parc_ids[i] + '.lut')
                    codes = tab_df['index'].tolist()
                    names = tab_df['name'].tolist()
                    colors = tab_df['color'].tolist()
                    cltparc.Parcellation.write_luttable(codes=codes, names=names, colors=colors, out_file=out_file_lut, headerlines=parc_desc, force=True)
            else:
                if format == 'tsv':
                    out_file_tsv = os.path.join(str(out_dir), out_name + '_' + parc_ids[i] + '.tsv')
                    cltparc.Parcellation.write_tsvtable(tsv_df=tab_df,out_file=out_file_tsv, force=True)

                if format == 'lut':
                    out_file_lut = os.path.join(str(out_dir), out_name + '_' + parc_ids[i] + '.lut')
                    codes = tab_df['index'].tolist()
                    names = tab_df['name'].tolist()
                    colors = tab_df['color'].tolist()
                    cltparc.Parcellation.write_luttable(codes=codes, names=names, colors=colors, out_file=out_file_lut, headerlines=parc_desc, force=True)
    
    def _build_lut_header(self):
        """
        This method builds the header of the LUT file.
        
        """
        
        # Detecting the base directory
        chim_dir = os.path.dirname(os.path.abspath(__file__))
    
        # Get the absolute of this file
        parc_json = os.path.join(chim_dir, 'config', 'supraregions_dictionary.json')

        # Reading the parcellation dictionary
        with open(parc_json) as f:
            parc_dict = json.load(f)
        
        # Reading the names of the supra-regions
        supra_names = list(parc_dict.keys())
        
        # Reading the parcellation code
        chim_code = self.parc_code
        
        # Creating the header lines
        headerlines = [" # Chimera parcellation code: {}".format(self.parc_code)]
        for i, supra in enumerate(supra_names):
            tmp_dict = parc_dict[supra]
            
            # Check if the parcellation code is in the dictionary
            if chim_code[i] in tmp_dict.keys():
                tmp_dict = tmp_dict[chim_code[i]]
                
                if tmp_dict["description"].endswith('.'):
                    tmp_dict["description"] = tmp_dict["description"][:-1]
                
                cite = "{} {}.". format (tmp_dict["atlas"], tmp_dict["citation"])  
                glob_desc = tmp_dict["description"] + ". Name: " + cite              
                headerlines.append("    " + glob_desc)
            else: 
                headerlines.append("    # {}. The parcellation code {} is not present in the dictionary for the supra-region {}.".format(i+1, chim_code[i], supra))
        
        return headerlines
        
    def _build_parcellation(self, t1:str, bids_dir:str, 
                            deriv_dir:str = None,
                            fssubj_dir:str = None,
                            growwm:Union[str, int] = None, 
                            bool_mixwm:bool = False,
                            force:bool = False):
        """
        This method builds the parcellation for the selected parcellation code.
        
        Parameters
        ----------
        bids_dir : str
            BIDs dataset directory.
        
        deriv_dir : str
            BIDs derivative directory.
            
            fssub
        
        growwm : str or int
            Grow of GM labels inside the white matter in mm.
        
        force : bool
            Overwrite the results.
        
        """
        global pipe_dict
        
        # Detecting the base directory
        cwd = os.path.dirname(os.path.abspath(__file__))
        chim_dir = os.path.dirname(cwd)
        
        if not os.path.isfile(t1):
            raise ValueError("Please provide a valid T1 image.")

        # Getting the entities from the name 
        anat_dir = os.path.dirname(t1)
        t1_name = os.path.basename(t1)
        ent_dict = cltbids._str2entity(t1_name)
        
        temp_entities = t1_name.split('_')[:-1]
        fullid = "_".join(temp_entities)
        ent_dict_fullid = cltbids._str2entity(fullid)
        
        if 'ses' in ent_dict.keys():
            path_cad       = "sub-" + ent_dict["subject"] + os.path.sep + "ses-" + ent_dict["ses"]
        else:
            path_cad       = "sub-" + ent_dict["subject"]
            
        # Creating Chimera directories
        if deriv_dir is None:
            chim_dir = os.path.join(bids_dir, 'derivatives', 'chimera', path_cad, 'anat')
        else:
            chim_dir = os.path.join(deriv_dir,'chimera', path_cad, 'anat')
        
        # Create the Chimera directory if it does not exist
        chim_dir = Path(chim_dir)
        chim_dir.mkdir(parents=True, exist_ok=True)
        
        supra_names = list(self.parc_dict.keys())
        # Detecting if Cortical is on the list of supra-regions
        bool_ctx = False
        if 'Cortical' in supra_names:
            # Remove it from the list
            supra_names.remove('Cortical')
            bool_ctx = True
        
        # ----------- Veryfing the existence of the parcellations, otherwise, compute them  --------- #
        if force: 
            bool_chim_exist = False
        else:
            bool_chim_exist = True
            
            if bool_ctx:
                        
                # Atributes for the cortical parcellation
                atlas_names   = self.parc_dict["Cortical"]["parcels"]
                for at_name in atlas_names:
                    ## -------- Cortical parcellation for the left hemisphere ---------------
                    # Creating the name for the output file
                    
                    if growwm is None:
                        growwm = ['0']
                                
                    for ngrow in np.arange(len(growwm)):
                        if growwm[ngrow] == '0':
                            out_vol_name = fullid + '_' + at_name + '_dseg.nii.gz'
                        else:
                            
                            ent_dict = cltbids._str2entity(at_name)
                            if 'desc' in ent_dict.keys():
                                if bool_mixwm:
                                    ent_dict["desc"] = ent_dict["desc"] + 'grow' + str(growwm[ngrow]) + 'mm+mixwm'
                                else:
                                    ent_dict["desc"] = ent_dict["desc"] + 'grow' + str(growwm[ngrow]) + 'mm'
                                tmp_str = cltbids._entity2str(ent_dict)
                                out_vol_name = fullid + '_' + tmp_str + '_dseg.nii.gz'
                            else:
                                if bool_mixwm:
                                    out_vol_name = fullid + '_' + at_name + '_desc-grow' + str(growwm[ngrow]) + 'mm+mixwm_dseg.nii.gz'
                                else:
                                    out_vol_name = fullid + '_' + at_name + '_desc-grow' + str(growwm[ngrow]) + 'mm_dseg.nii.gz'

                        
                        chim_parc_name = cltbids._replace_entity_value(out_vol_name, {"atlas": "chimera" + chim_code} ) 
                        chim_parc_file = os.path.join(str(chim_dir), chim_parc_name)
                        chim_parc_lut  = os.path.join(str(chim_dir), cltbids._replace_entity_value(chim_parc_name, {"extension": "lut"}))
                        chim_parc_tsv  = os.path.join(str(chim_dir), cltbids._replace_entity_value(chim_parc_name, {"extension": "tsv"}))
                        
                        if not os.path.isfile(chim_parc_file) or not os.path.isfile(chim_parc_lut) or not os.path.isfile(chim_parc_tsv)  or force:
                            bool_chim_exist = False
                            
            else:
                out_vol_name = fullid + '_dseg.nii.gz'

                chim_parc_name = cltbids._insert_entity(out_vol_name, {"atlas": "chimera" + chim_code} ) 
                chim_parc_file = os.path.join(str(chim_dir), chim_parc_name)
                chim_parc_lut  = os.path.join(str(chim_dir), cltbids._replace_entity_value(chim_parc_name, {"extension": "lut"}))
                chim_parc_tsv  = os.path.join(str(chim_dir), cltbids._replace_entity_value(chim_parc_name, {"extension": "tsv"}))
                
                if not os.path.isfile(chim_parc_file) or not os.path.isfile(chim_parc_lut) or not os.path.isfile(chim_parc_tsv) or force:
                    bool_chim_exist = False
        
        # ------- End of veryfing the existence of the parcellations, otherwise, compute them  --------- #
        
        # Run chimera if the desired parcellation does not exist
        if not bool_chim_exist:
        
            ######## ----------- Detecting FREESURFER_HOME directory ------------- #
            if pipe_dict["packages"]["freesurfer"]["cont_tech"] != 'local':
                cont_tech = pipe_dict["packages"]["freesurfer"]["cont_tech"]
                cont_image = pipe_dict["packages"]["freesurfer"]["container"]
                
                cmd_bashargs = ['echo', '$FREESURFER_HOME']
                cmd_cont = cltmisc._generate_container_command(cmd_bashargs, cont_tech, cont_image) # Generating container command
                out_cmd = subprocess.run(cmd_cont, stdout=subprocess.PIPE, universal_newlines=True)
                fslut_file_cont = os.path.join(out_cmd.stdout.split('\n')[0], 'FreeSurferColorLUT.txt')
                tmp_name = str(uuid.uuid4())
                cmd_bashargs = ['cp', 'replace_cad', '/tmp/' + tmp_name]
                cmd_cont = cltmisc._generate_container_command(cmd_bashargs, cont_tech, cont_image)
                
                # Replace the element of the list equal to replace_cad by the path of the lut file
                cmd_cont = [w.replace('replace_cad', fslut_file_cont) for w in cmd_cont]
                subprocess.run(cmd_cont, stdout=subprocess.PIPE, universal_newlines=True)
                fslut_file = os.path.join('/tmp', tmp_name)
                
                ######## ------------- Reading FreeSurfer color lut table ------------ #
                lut_dict= cltparc.Parcellation.read_luttable(fslut_file)
                
                os.remove(fslut_file)
                
            else:

                fshome_dir = os.getenv('FREESURFER_HOME')
                fslut_file = os.path.join(fshome_dir, 'FreeSurferColorLUT.txt')
                
                ######## ------------- Reading FreeSurfer color lut table ------------ #
                lut_dict= cltparc.Parcellation.read_luttable(fslut_file)
                
            # Extracting the information from the lut dictionary
            st_codes = lut_dict['index']
            st_names = lut_dict['name']
            st_colors = lut_dict['color']
            
            
            ######## ----- Running FreeSurfer if it was not previously computed ------ #
            sub2proc = cltfree.FreeSurferSubject(fullid, subjs_dir=fssubj_dir)

            cont_tech_freesurfer  = pipe_dict["packages"]["freesurfer"]["cont_tech"]
            cont_image_freesurfer = pipe_dict["packages"]["freesurfer"]["container"]
            cont_tech_ants        = pipe_dict["packages"]["ants"]["cont_tech"]
            cont_image_ants       = pipe_dict["packages"]["ants"]["container"]
            cont_tech_fsl         = pipe_dict["packages"]["fsl"]["cont_tech"]
            cont_image_fsl        = pipe_dict["packages"]["fsl"]["container"]
            
            # Running FreeSurfer if it was not previously computed is mandatory
            sub2proc._launch_freesurfer(force=force, 
                                        t1w_img=t1,
                                            cont_tech=cont_tech_freesurfer, 
                                            cont_image=cont_image_freesurfer)
            

            nii_image = os.path.join(sub2proc.subjs_dir, sub2proc.subj_id, 'tmp', 'aparc+aseg.nii.gz')
            mgz_image = os.path.join(sub2proc.subjs_dir, sub2proc.subj_id, 'mri', 'aparc+aseg.mgz')
                    
            if not os.path.isfile(nii_image):
                sub2proc._conform2native(mgz_conform=mgz_image,
                                        nii_native=nii_image,
                                        cont_image=cont_image_freesurfer,
                                        cont_tech=cont_tech_freesurfer,
                                        force=force)
                
            if "aseg_parc" not in locals():
                aseg_parc = cltparc.Parcellation(parc_file=nii_image)
                aseg_parc.index = st_codes
                aseg_parc.name = st_names
                aseg_parc.color = st_colors
                aseg_parc._adjust_values()
                
            # Creating the parcellation for the extra regions
            extra_parc = _create_extra_regions_parc(aparc=nii_image)
                
            # Remove the nifti file
            os.remove(nii_image)
            
            # Building the main header information for the LUT file
            glob_header_info = self._build_lut_header()

            # Processing that will be perfomed for multiple supra-regions
            gm_sub_names = list(self.parc_dict.keys())
            
            # Remove 'Cortical', 'GyralWM' and 'WhiteMatter' from the list
            if 'Cortical' in gm_sub_names:
                gm_sub_names.remove('Cortical')
            
            if 'GyralWM' in gm_sub_names:
                gm_sub_names.remove('GyralWM')
                
            if 'WhiteMatter' in supra_names:
                gm_sub_names.remove('WhiteMatter')
                
            # Taking the image dimensions and the affine matrix in native space
            t1_image = nib.load(t1)
            affine   =  t1_image.affine
            
            # Create a numpy array with the same dimensions of the T1 image and fill it with zeros. 
            # The array will be used to store the parcellation. The elements should be integers.
            ref_image = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)
            lh2refill = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)
            rh2refill = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)
            
            # Creating the parcellation objects
            lh_parc  = cltparc.Parcellation(parc_file=ref_image, affine=affine) # It will include the parcellation for the left hemisphere
            rh_parc  = copy.deepcopy(lh_parc) # It will include the parcellation for the right hemisphere
            mid_parc = copy.deepcopy(lh_parc) # It will include the parcellation for structures that do not belong to any hemisphere
            
            files2del = [] # Temporal files that will be deleted
            exec_cmds = []
            for supra in gm_sub_names:
                
                # Getting the information of the common atributes
                atlas_code    = self.parc_dict[supra]["code"]
                atlas_str     = self.parc_dict[supra]["atlas"]
                atlas_desc    = self.parc_dict[supra]["description"]
                atlas_cita    = self.parc_dict[supra]["citation"]
                atlas_src     = self.parc_dict[supra]["source"]
                atlas_ref     = self.parc_dict[supra]["reference"]
                atlas_parcs   = self.parc_dict[supra]["parcels"]
                atlas_mask    = self.parc_dict[supra]["mask"]
                atlas_type    = self.parc_dict[supra]["type"]
                deriv_fold    = self.parc_dict[supra]["deriv_volfold"]
                proc_dict     = self.parc_dict[supra]["processing"]
                
                
                if proc_dict["method"] == 'comform2native':
                    
                    if proc_dict["labels"] == 'freesurferextra':
                        sub2proc._launch_freesurfer(force=force, 
                                                    extra_proc=supra.lower(),
                                                    cont_tech=cont_tech_freesurfer, 
                                                    cont_image=cont_image_freesurfer)
                        
                        
                        fsextra_files = glob(os.path.join(sub2proc.subjs_dir, sub2proc.subj_id, 'mri', '*' + atlas_parcs + '.mgz'))
                        if len(fsextra_files) == 0:
                            raise ValueError("The Freesurfer extra parcellation was not found.")
                        
                        elif len(fsextra_files) >  1:
                            lh_mgz_image = os.path.join(sub2proc.subjs_dir, sub2proc.subj_id, 'mri', 'lh.' + atlas_parcs + '.mgz')
                            rh_mgz_image = os.path.join(sub2proc.subjs_dir, sub2proc.subj_id, 'mri', 'rh.' + atlas_parcs + '.mgz')
                            lh_nii_image = os.path.join(deriv_dir, deriv_fold, path_cad, 'anat', fullid + '_hemi-L_atlas-' + atlas_str + '_dseg.nii.gz')
                            rh_nii_image = os.path.join(deriv_dir, deriv_fold, path_cad, 'anat', fullid + '_hemi-R_atlas-' + atlas_str + '_dseg.nii.gz')
                            
                            if not os.path.isfile(lh_nii_image):
                                dir_name = os.path.dirname(lh_nii_image)
                                dir_name = Path(dir_name)
                                dir_name.mkdir(parents=True, exist_ok=True)
                                    
                                if atlas_ref == 'conform':
                                    sub2proc._conform2native(mgz_conform=lh_mgz_image,
                                                            nii_native=lh_nii_image,
                                                            cont_image=cont_image_freesurfer,
                                                            cont_tech=cont_tech_freesurfer,
                                                            force=force)
                                else:
                                    lh_nii_image = lh_mgz_image
                                    
                            if not os.path.isfile(rh_nii_image):
                                dir_name = os.path.dirname(rh_nii_image)
                                dir_name = Path(dir_name)
                                dir_name.mkdir(parents=True, exist_ok=True)
                                
                                if atlas_ref == 'conform':
                                    sub2proc._conform2native(mgz_conform=rh_mgz_image,
                                                            nii_native=rh_nii_image,
                                                            cont_image=cont_image_freesurfer,
                                                            cont_tech=cont_tech_freesurfer,
                                                            force=force)
                                else:
                                    rh_nii_image = rh_mgz_image
                                    
                            lh_supra_parc = cltparc.Parcellation(parc_file=lh_nii_image)
                            lh_supra_parc.index  = self.supra_dict[supra][supra][atlas_code]["lh"]["index"]
                            lh_supra_parc.name   = self.supra_dict[supra][supra][atlas_code]["lh"]["name"]
                            lh_supra_parc.color  = self.supra_dict[supra][supra][atlas_code]["lh"]["color"]
                            lh_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['lh']['index'])
                            lh_supra_parc._export_colortable(out_file=lh_nii_image.replace('.nii.gz', '.lut'), lut_type="lut")
                            lh_supra_parc._export_colortable(out_file=lh_nii_image.replace('.nii.gz', '.tsv'), lut_type="tsv")
                                    
                            rh_supra_parc = cltparc.Parcellation(parc_file=rh_nii_image)
                            rh_supra_parc.index  = self.supra_dict[supra][supra][atlas_code]["rh"]["index"]
                            rh_supra_parc.name   = self.supra_dict[supra][supra][atlas_code]["rh"]["name"]
                            rh_supra_parc.color  = self.supra_dict[supra][supra][atlas_code]["rh"]["color"]
                            rh_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['rh']['index'])
                            rh_supra_parc._export_colortable(out_file=rh_nii_image.replace('.nii.gz', '.lut'), lut_type="lut")
                            rh_supra_parc._export_colortable(out_file=rh_nii_image.replace('.nii.gz', '.tsv'), lut_type="tsv")
                                    
                            
                        elif len(fsextra_files) ==  1:
                            mgz_image = fsextra_files[0]
                            nii_image = os.path.join(deriv_dir, deriv_fold, path_cad, 'anat', fullid + '_atlas-' + atlas_str + '_dseg.nii.gz')
                            if not os.path.isfile(nii_image):
                                dir_name = os.path.dirname(nii_image)
                                dir_name = Path(dir_name)
                                dir_name.mkdir(parents=True, exist_ok=True)
                                
                                if atlas_ref == 'conform':
                                    sub2proc._conform2native(mgz_conform=mgz_image,
                                                            nii_native=nii_image,
                                                            cont_image=cont_image_freesurfer,
                                                            cont_tech=cont_tech_freesurfer,
                                                            force=force)
                                else:
                                    nii_image = mgz_image
                            
                            tmp_parc = cltparc.Parcellation(parc_file=nii_image)
                            index, name, color = _mix_side_prop(self.supra_dict[supra][supra][atlas_code])
                            tmp_parc.index = index
                            tmp_parc.name = name
                            tmp_parc.color = color
                            
                            tmp_parc._export_colortable(out_file=nii_image.replace('.nii.gz', '.lut'), lut_type="lut")
                            tmp_parc._export_colortable(out_file=nii_image.replace('.nii.gz', '.tsv'), lut_type="tsv")
                            
                            # Left Hemisphere
                            if 'lh' in self.supra_dict[supra][supra][atlas_code].keys():
                                lh_supra_parc  = copy.deepcopy(tmp_parc)
                                lh_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['lh']['index'])
                                        
                            # Right Hemisphere
                            if 'rh' in self.supra_dict[supra][supra][atlas_code].keys():
                                rh_supra_parc  = copy.deepcopy(tmp_parc)
                                rh_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['rh']['index'])
                            
                            # Non-hemispheric structures
                            if 'mid' in self.supra_dict[supra][supra][atlas_code].keys():
                                mid_supra_parc = copy.deepcopy(tmp_parc)
                                mid_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['mid']['index'])    

                    else:
                    
                        if atlas_src == 'freesurfer':
                            nii_image = os.path.join(sub2proc.subjs_dir, sub2proc.subj_id, 'tmp', atlas_parcs + '.nii.gz')
                            mgz_image = os.path.join(sub2proc.subjs_dir, sub2proc.subj_id, 'mri', atlas_parcs + '.mgz')
                        
                        if 'aseg_parc' not in locals():
                            if not os.path.isfile(nii_image):
                                sub2proc._conform2native(mgz_conform=mgz_image,
                                                            nii_native=nii_image,
                                                            cont_image=cont_image_freesurfer,
                                                            cont_tech=cont_tech_freesurfer,
                                                            force=force)
                                files2del.append(nii_image)
                                
                            tmp_parc = cltparc.Parcellation(parc_file=nii_image)
                            
                        else:
                            tmp_parc = copy.deepcopy(aseg_parc)

                        # Left Hemisphere
                        if 'lh' in self.supra_dict[supra][supra][atlas_code].keys():
                            lh_supra_parc        = copy.deepcopy(tmp_parc)
                            lh_supra_parc.index  = self.supra_dict[supra][supra][atlas_code]["lh"]["index"]
                            lh_supra_parc.name   = self.supra_dict[supra][supra][atlas_code]["lh"]["name"]
                            lh_supra_parc.color  = self.supra_dict[supra][supra][atlas_code]["lh"]["color"]
                            lh_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['lh']['index'])
                            
                        # Right Hemisphere
                        if 'rh' in self.supra_dict[supra][supra][atlas_code].keys():
                            rh_supra_parc        = copy.deepcopy(tmp_parc)
                            rh_supra_parc.index  = self.supra_dict[supra][supra][atlas_code]["rh"]["index"]
                            rh_supra_parc.name   = self.supra_dict[supra][supra][atlas_code]["rh"]["name"]
                            rh_supra_parc.color  = self.supra_dict[supra][supra][atlas_code]["rh"]["color"]
                            rh_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['rh']['index'])

                        # Non-hemispheric structures
                        if 'mid' in self.supra_dict[supra][supra][atlas_code].keys():
                            mid_supra_parc        = copy.deepcopy(tmp_parc)
                            mid_supra_parc.index  = self.supra_dict[supra][supra][atlas_code]["mid"]["index"]
                            mid_supra_parc.name   = self.supra_dict[supra][supra][atlas_code]["mid"]["name"]
                            mid_supra_parc.color  = self.supra_dict[supra][supra][atlas_code]["mid"]["color"]
                            mid_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['mid']['index'])    
                        

                elif proc_dict["method"] == None:

                    # Running FIRST if it is needed
                    if atlas_code =='R':
                        fsl_outdir = os.path.join(deriv_dir, deriv_fold, path_cad, 'anat')
                        first_nii = os.path.join(str(fsl_outdir), fullid + '_atlas-' + atlas_str + '_dseg.nii.gz' )

                        if "first_parc" not in locals():
                            fsl_outdir = Path(fsl_outdir)
                            fsl_outdir.mkdir(parents=True, exist_ok=True)
                        
                            # Running the FIRST
                            _launch_fsl_first(t1,
                                                first_parc = first_nii,
                                                cont_tech = cont_tech_fsl, 
                                                cont_image = cont_image_fsl, 
                                                force=force)
                            first_parc = cltparc.Parcellation(parc_file=first_nii)
                        
                        tmp_parc = copy.deepcopy(first_parc)
                        index, name, color = _mix_side_prop(self.supra_dict[supra][supra][atlas_code])
                        tmp_parc.index = index
                        tmp_parc.name = name
                        tmp_parc.color = color
                        
                        tmp_parc._export_colortable(out_file=first_nii.replace('.nii.gz', '.lut'), lut_type="lut")
                        tmp_parc._export_colortable(out_file=first_nii.replace('.nii.gz', '.tsv'), lut_type="tsv")
                    # Left Hemisphere
                    if 'lh' in self.supra_dict[supra][supra][atlas_code].keys():
                        lh_supra_parc  = copy.deepcopy(tmp_parc)
                        lh_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['lh']['index'])
                            
                    # Right Hemisphere
                    if 'rh' in self.supra_dict[supra][supra][atlas_code].keys():
                        rh_supra_parc  = copy.deepcopy(tmp_parc)
                        rh_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['rh']['index'])
                    
                    # Non-hemispheric structures
                    if 'mid' in self.supra_dict[supra][supra][atlas_code].keys():
                        mid_supra_parc = copy.deepcopy(tmp_parc)
                        mid_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['mid']['index'])    
                    
                elif proc_dict["method"] == 'atlasbased':
                    t1_temp = proc_dict["reference"]
                    atlas = proc_dict["labels"]
                    atlas_type = proc_dict["type"]
                    
                    # Basename for transformations
                    spat_tf_ent = ent_dict_fullid.copy()
                    spat_tf_ent = cltbids._delete_entity(spat_tf_ent, ["space", "desc", "suffix", "extension"])
                    spat_tf_ent["from"] = "T1w"
                    spat_tf_ent["to"] = atlas_ref
                    spat_tf_ent["mode"] = "image"
                    spat_tf_ent["suffix"] = "xfm"
                    spat_tf_ent["extension"] = "mat"
                    xfm_base = os.path.join(deriv_dir, pipe_dict["outputs"]["transforms"], path_cad, 'anat', cltbids._entity2str(spat_tf_ent))
                    work_dir = os.path.join(deriv_dir, deriv_fold, path_cad, 'anat') 
                    out_parc_spam = os.path.join(work_dir, fullid + '_atlas-' + atlas_str + '_probseg.nii.gz')
                    out_parc_maxp = os.path.join(work_dir, fullid + '_atlas-' + atlas_str + '_dseg.nii.gz')

                    if not os.path.isfile(out_parc_maxp) or force:
                        work_dir.mkdir(parents=True, exist_ok=True) 
                        
                        # Detecting the side
                        sides_ids = list(self.supra_dict[supra][supra][atlas_code].keys())
                        sides_ids = sorted(sides_ids, key=lambda x: not ("lh" in x or "rh" in x))

                        # Masking the cerebellum from T1w image
                        tmp_t1 = t1
                        if supra == 'Cerebellum':
                            if self.parc_dict[supra]['name'] == 'SUIT':
                                tmp_t1 = os.path.join(str(work_dir), 'tmp_cerb_bs.nii.gz')
                                files2del.append(tmp_t1)
                                
                                cltimg.crop_image_from_mask(t1,
                                                                aseg_parc.data,
                                                                tmp_t1,
                                                                [7, 8, 16, 46, 47, 15, 16])
                            
                        if atlas_type == 'spam':
                            cltseg._abased_parcellation(t1, 
                                                    t1_temp, 
                                                    atlas, 
                                                    out_parc_spam, 
                                                    xfm_base,
                                                    cont_tech=cont_tech_ants,
                                                    cont_image=cont_image_ants)
                            
                            # Detecting the side
                            sides_id = list(self.supra_dict[supra][supra][atlas_code].keys())
                            
                            for side_cont, side in enumerate(sides_id):
                                vol_indexes = np.array(self.supra_dict[supra][supra][atlas_code][side]['index'])-1
                                tmp_par_file = os.path.join(work_dir, fullid + '_hemi-' + side + '_atlas-' + atlas_str + '_dseg.nii.gz')
                                files2del.append(tmp_par_file)
                                
                                tmp_parc_file = cltseg._spams2maxprob(out_parc_spam, 
                                                            prob_thresh=0.05, 
                                                            vol_indexes=vol_indexes,
                                                            maxp_name=tmp_par_file)
                                
                                tmp_parc = cltparc.Parcellation(parc_file=tmp_parc_file)
                                tmp_parc.index = vol_indexes + 1
                                tmp_parc.name  = self.supra_dict[supra][supra][atlas_code][side]['name']
                                tmp_parc.color = self.supra_dict[supra][supra][atlas_code][side]['color']
                                
                                if side in self.supra_dict[supra][supra]['F'].keys():
                                    aseg_code    = self.supra_dict[supra][supra]['F'][side]['index']
                                    tmp_parc._apply_mask(image_mask=aseg_parc, codes2mask=aseg_code)
                                
                                if side_cont == 0:  
                                    def_parc = copy.deepcopy(tmp_parc)
                                else:
                                    def_parc._add_parcellation(tmp_parc)
                                    
                                # Removing the temporal side images
                                if os.path.isfile(tmp_parc_file):
                                    os.remove(tmp_parc_file)
                                
                            def_parc._save_parcellation(out_file= out_parc_maxp, affine=def_parc.affine, save_lut=True, save_tsv=True)
                            
                        elif atlas_type == 'maxprob':
                            
                            cltseg._abased_parcellation(t1, 
                                                    t1_temp, 
                                                    atlas, 
                                                    out_parc_maxp, 
                                                    xfm_base, 
                                                    atlas_type='maxprob',
                                                    cont_tech=cont_tech_ants,
                                                    cont_image=cont_image_ants)

                    
                    tmp_parc = cltparc.Parcellation(parc_file=out_parc_maxp)
                    index, name, color = _mix_side_prop(self.supra_dict[supra][supra][atlas_code])
                    tmp_parc.index = index
                    tmp_parc.name = name
                    tmp_parc.color = color
                    
                    tmp_parc._export_colortable(out_file=out_parc_maxp.replace('.nii.gz', '.lut'), lut_type="lut")
                    tmp_parc._export_colortable(out_file=out_parc_maxp.replace('.nii.gz', '.tsv'), lut_type="tsv")
                    # Left Hemisphere
                    if 'lh' in self.supra_dict[supra][supra][atlas_code].keys():
                        lh_supra_parc  = copy.deepcopy(tmp_parc)
                        lh_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['lh']['index'])
                            
                    # Right Hemisphere
                    if 'rh' in self.supra_dict[supra][supra][atlas_code].keys():
                        rh_supra_parc  = copy.deepcopy(tmp_parc)
                        rh_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['rh']['index'])
                    
                    # Non-hemispheric structures
                    if 'mid' in self.supra_dict[supra][supra][atlas_code].keys():
                        mid_supra_parc = copy.deepcopy(tmp_parc)
                        mid_supra_parc._keep_by_code(codes2look=self.supra_dict[supra][supra][atlas_code]['mid']['index'])    

                
                # Appending the parcellations
                if "lh_supra_parc" in locals():
                    lh_supra_parc._rearange_parc()
                    
                    if 'F' in self.supra_dict[supra][supra].keys():
                        # Use the FreeSurfer parcellation to detect the voxels that are not in the lh_supra_parc
                        lh2refill_parc          = copy.deepcopy(aseg_parc)
                        lh2refill_parc._keep_by_code(codes2look=self.supra_dict[supra][supra]['F']['lh']['index'])
                        
                        # Find the voxels that are not in the lh_supra_parc and are in the lh2refill
                        ind = np.where((lh_supra_parc.data == 0) & (lh2refill_parc.data != 0))
                        lh2refill[ind] = 1
                        del lh2refill_parc
                    
                    # Add the parcellation to the global left subcortical parcellation                
                    lh_parc._add_parcellation(lh_supra_parc, append=True)
                    nlh_subc = len(lh_parc.index)
                    del lh_supra_parc
                    # lh_parc._save_parcellation(out_file= '/home/yaleman/lh_test.nii.gz', save_lut=True)
                    
                if "rh_supra_parc" in locals():
                    rh_supra_parc._rearange_parc()
                    
                    if 'F' in self.supra_dict[supra][supra].keys():
                        # Use the FreeSurfer parcellation to detect the voxels that are not in the lh_supra_parc
                        rh2refill_parc          = copy.deepcopy(aseg_parc)
                        rh2refill_parc._keep_by_code(codes2look=self.supra_dict[supra][supra]['F']['rh']['index'])
                        
                        # Find the voxels that are not in the lh_supra_parc and are in the lh2refill
                        ind = np.where((rh_supra_parc.data == 0) & (rh2refill_parc.data != 0))
                        rh2refill[ind] = 1
                        del rh2refill_parc
                    
                    # Add the parcellation to the global right subcortical parcellation
                    rh_parc._add_parcellation(rh_supra_parc, append=True)
                    nrh_subc = len(rh_parc.index)
                    del rh_supra_parc
                    # rh_parc._save_parcellation(out_file= '/home/yaleman/rh_test.nii.gz', save_lut=True)
                    
                if 'mid_supra_parc' in locals():
                    mid_supra_parc._rearange_parc()
                    mid_parc._add_parcellation(mid_supra_parc, append=True)
                    del mid_supra_parc
                    # mid_parc._save_parcellation(out_file= '/home/yaleman/mid_test.nii.gz', save_lut=True)

            # Detecting the number of regions 
            if "lh_parc" in locals():
                nlh_subc = len(lh_parc.index)
                    
            if "rh_parc" in locals():
                nrh_subc = len(rh_parc.index)
                    
            if 'mid_parc' in locals():
                nmid_subc = len(mid_parc.index)
            
            # if "WhiteMatter" in supra_names:
                #self.supra_dict[supra][supra][atlas_code]["mid"]["index"]
                
            date_time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
            if bool_ctx:
                
                # Atributes for the cortical parcellation
                atlas_names   = self.parc_dict["Cortical"]["parcels"]
                
                proc_dict     = self.parc_dict["Cortical"]["processing"]
                ctx_meth      = proc_dict["method"]
                
                nctx_parc     = len(self.parc_dict["Cortical"]["processing"]["labels"]["lh"])
                for c in np.arange(nctx_parc):
                    
                    # Temporal header lines
                    glob_header_info_tmp = copy.deepcopy(glob_header_info)
                    
                    ## -------- Cortical parcellation for the left hemisphere ---------------
                    # Creating the name for the output file
                    lh_in_parc = self.parc_dict["Cortical"]["processing"]["labels"]["lh"][c]
                    at_name = [s for s in atlas_names if s in lh_in_parc]
                    lh_out_annot = os.path.join(deriv_dir, 
                                                self.parc_dict["Cortical"]["deriv_surffold"], 
                                                path_cad, 'anat',
                                                fullid + '_hemi-L' + '_' + ''.join(at_name) + '_dseg.annot')
                    
                    ## -------- Cortical parcellation for the right hemisphere ---------------
                    # Creating the name for the output file
                    rh_in_parc = self.parc_dict["Cortical"]["processing"]["labels"]["rh"][c]
                    at_name = [s for s in atlas_names if s in rh_in_parc]
                    at_name = ''.join(at_name)
                    rh_out_annot = os.path.join(deriv_dir, 
                                                self.parc_dict["Cortical"]["deriv_surffold"], 
                                                path_cad, 'anat',     
                                                fullid + '_hemi-R' + '_' + at_name + '_dseg.annot')
                    
                    if ctx_meth == 'annot2indiv':
                        # Moving to individual space
                        sub2proc._annot2ind(ref_id=self.parc_dict["Cortical"]["processing"]["reference"], 
                                        hemi='lh', 
                                        fs_annot=lh_in_parc, 
                                        ind_annot=lh_out_annot, 
                                        cont_tech = cont_tech_freesurfer,
                                        cont_image=cont_image_freesurfer, 
                                        force=force)
                        
                        sub2proc._annot2ind(ref_id=self.parc_dict["Cortical"]["processing"]["reference"], 
                                        hemi='rh', 
                                        fs_annot=rh_in_parc, 
                                        ind_annot=rh_out_annot,
                                        cont_tech = cont_tech_freesurfer,
                                        cont_image=cont_image_freesurfer, 
                                        force=force)
                        
                    if ctx_meth == 'gcs2indiv':
                        # Moving to individual space
                        sub2proc._gcs2ind(fs_gcs=lh_in_parc, 
                                        hemi='lh', 
                                        ind_annot=lh_out_annot, 
                                        cont_tech=cont_tech_freesurfer,
                                        cont_image=cont_image_freesurfer,
                                        force=force)
                        
                        sub2proc._gcs2ind(fs_gcs=rh_in_parc, 
                                        hemi='rh', 
                                        ind_annot=rh_out_annot, 
                                        cont_tech=cont_tech_freesurfer,
                                        cont_image=cont_image_freesurfer,
                                        force=force)
                    
                    # Copying to the labels folder
                    temp_lh = os.path.join(sub2proc.subjs_dir, sub2proc.subj_id, 'label',  'lh.' + at_name + '.annot')
                    shutil.copyfile(lh_out_annot, temp_lh)
                    
                    # Copying to the labels folder
                    temp_rh = os.path.join(sub2proc.subjs_dir, sub2proc.subj_id, 'label',  'rh.' + at_name + '.annot')
                    shutil.copyfile(rh_out_annot, temp_rh)
                    
                    ## -------- Creating the volumetric parcellation ---------------
                    out_vol_dir = os.path.join(deriv_dir, self.parc_dict["Cortical"]["deriv_volfold"], path_cad, 'anat')
                    if growwm is None:
                        growwm = ['0']
                    
                    ent_dict = cltbids._str2entity(at_name)
                    if "scale" in ent_dict.keys():
                            scale_cad = 'Scale: {}'.format(ent_dict["scale"])
                    else:
                        scale_cad = None
                        
                    if "seg" in ent_dict.keys():
                        seg_cad = 'Segmentation: {}'.format(ent_dict["seg"])
                    else:
                        seg_cad = None

                    if scale_cad is not None or seg_cad is not None:
                        if scale_cad is not None and seg_cad is not None:
                            cad2add = '. ' + scale_cad + ' - ' + seg_cad 
                            
                        elif scale_cad is not None and seg_cad is None:
                            cad2add = '. ' + scale_cad 
                            
                        elif scale_cad is None and seg_cad is not None:
                            cad2add = '. ' + seg_cad
                        
                        glob_header_info_tmp[0] = glob_header_info_tmp[0] + cad2add
                        
                                
                    for ngrow in np.arange(len(growwm)):
                        if growwm[ngrow] == '0':
                            out_vol_name = fullid + '_' + at_name + '_dseg.nii.gz'
                        else:
                            
                            ent_dict = cltbids._str2entity(at_name)
                            if 'desc' in ent_dict.keys():
                                if bool_mixwm:
                                    ent_dict["desc"] = ent_dict["desc"] + 'grow' + str(growwm[ngrow]) + 'mm+mixwm'
                                else:
                                    ent_dict["desc"] = ent_dict["desc"] + 'grow' + str(growwm[ngrow]) + 'mm'
                                tmp_str = cltbids._entity2str(ent_dict)
                                out_vol_name = fullid + '_' + tmp_str + '_dseg.nii.gz'
                            else:
                                if bool_mixwm:
                                    out_vol_name = fullid + '_' + at_name + '_desc-grow' + str(growwm[ngrow]) + 'mm+mixwm_dseg.nii.gz'
                                else:
                                    out_vol_name = fullid + '_' + at_name + '_desc-grow' + str(growwm[ngrow]) + 'mm_dseg.nii.gz'

                        sub2proc._surf2vol(atlas=at_name, 
                                            out_vol=os.path.join(out_vol_dir, out_vol_name), 
                                            gm_grow=growwm[ngrow], 
                                            bool_mixwm = bool_mixwm,
                                            force=False, 
                                            bool_native=True, 
                                            color_table=['tsv', 'lut'],
                                            cont_tech=cont_tech_freesurfer,
                                            cont_image=cont_image_freesurfer)
                        
                        chim_parc_name = cltbids._replace_entity_value(out_vol_name, {"atlas": "chimera" + chim_code} ) 
                        chim_parc_file = os.path.join(str(chim_dir), chim_parc_name)
                        chim_parc_lut  = os.path.join(str(chim_dir), cltbids._replace_entity_value(chim_parc_name, {"extension": "lut"}))
                        chim_parc_tsv  = os.path.join(str(chim_dir), cltbids._replace_entity_value(chim_parc_name, {"extension": "tsv"}))
                        
                        # Creating the first part of the headers
                        part_header = ['# $Id: {} {} \n'.format(chim_parc_lut, date_time)]
                        
                        part_header.append('# Corresponding parcellation: {} \n'.format(chim_parc_file))
                        
                        lut_header = part_header + glob_header_info_tmp
                        lut_header = lut_header + ['\n']
                        lut_header.append('{:<4} {:<50} {:>3} {:>3} {:>3} {:>3}'.format("#No.", "Label Name:", "R", "G", "B", "A"))

                        if not os.path.isfile(chim_parc_file) or not os.path.isfile(chim_parc_lut) or not os.path.isfile(chim_parc_tsv)  or force:
                            # Creating the joined parcellation
                            ref_image = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)
                            chim_parc  = cltparc.Parcellation(parc_file=ref_image, affine=affine)
                            
                            ctx_parc = cltparc.Parcellation(parc_file=os.path.join(out_vol_dir, out_vol_name))
                            ctx_parc._remove_by_name(names2remove=['unknown', 'medialwall', 'corpuscallosum'])
                            
                            lh_ctx_parc = copy.deepcopy(ctx_parc)    
                            rh_ctx_parc = copy.deepcopy(ctx_parc)
                            
                            lh_ctx_parc._keep_by_name(names2look='ctx-lh-')
                            nlh_ctx = len(lh_ctx_parc.index)
                            rh_ctx_parc._keep_by_name(names2look='ctx-rh-')
                            nrh_ctx = len(rh_ctx_parc.index)
                            
                            # Detect the global White Matter
                            brain_wm_parc = copy.deepcopy(ctx_parc)
                            brain_wm_parc._keep_by_code(codes2look=[2, 41, 5001, 5002, 7, 46])
                            ind = np.where(brain_wm_parc.data != 0)
                            brain_wm_parc.data[ind] = 1
                            brain_wm_parc.index = [1]
                            brain_wm_parc.name = ['wm-brain-whitematter']
                            brain_wm_parc.color = ['#ffffff']
                            brain_wm_parc._rearange_parc(offset=2999)
                            brain_wm_parc.data[np.where(lh2refill)] = 3000
                            brain_wm_parc.data[np.where(rh2refill)] = 3000
                            
                            # White Matter for the Right Hemisphere
                            tmp_rh = cltmisc._filter_by_substring(ctx_parc.name, 'wm-rh-')
                            if tmp_rh:
                                rh_wm_parc = copy.deepcopy(ctx_parc)
                                rh_wm_parc._keep_by_name(names2look=tmp_rh)
                                rh_wm_parc._rearange_parc(offset=3000)
                            
                            # White Matter for the Left Hemisphere
                            tmp_lh = cltmisc._filter_by_substring(ctx_parc.name, 'wm-lh-')
                            if tmp_lh:
                                lh_wm_parc = copy.deepcopy(ctx_parc)
                                lh_wm_parc._keep_by_name(names2look=tmp_lh)
                                lh_wm_parc._rearange_parc(offset=3000 + nrh_ctx + nrh_subc)
                            
                            # Adding the right cortical parcellation to the final image
                            rh_ctx_parc._rearange_parc()
                            chim_parc._add_parcellation(rh_ctx_parc, append=True)
                            del rh_ctx_parc
                            
                            # Adding the right non-cortical parcellation to the final image
                            if "rh_parc" in locals():
                                chim_parc._add_parcellation(rh_parc, append=True)
                            
                            # Adding the left cortical parcellation to the final image
                            lh_ctx_parc._rearange_parc()
                            chim_parc._add_parcellation(lh_ctx_parc, append=True)
                            del lh_ctx_parc
                            
                            # Adding the left non-cortical parcellation to the final image
                            if "lh_parc" in locals():
                                chim_parc._add_parcellation(lh_parc, append=True)
                            
                            # Adding the regions that do not belong to any hemisphere to the final image
                            if "mid_parc" in locals():
                                chim_parc._add_parcellation(mid_parc, append=True)
                            
                            # Adding the white matter to the final image
                            chim_parc._add_parcellation(brain_wm_parc, append=False)
                            del brain_wm_parc
                                
                            if "rh_wm_parc" in locals():
                                chim_parc._add_parcellation(rh_wm_parc, append=False)
                                del rh_wm_parc
                                
                            if "lh_wm_parc" in locals():
                                chim_parc._add_parcellation(lh_wm_parc, append=False)
                                del lh_wm_parc
                            
                            # Adding the extra regions
                            if "extra_parc" in locals():
                                # Detecting if there is region overlap and removing it
                                tmp_extra = copy.deepcopy(extra_parc)  
                                mask = np.logical_and(tmp_extra.data  != 0, chim_parc.data != 0)
                                indexes = np.where(mask)
                                tmp_extra.data[indexes] = 0
                                tmp_extra._adjust_values()
                                chim_parc._add_parcellation(tmp_extra, append=False)
                                del tmp_extra
                                
                            # Saving the FINAL parcellation
                            chim_parc._save_parcellation(out_file=chim_parc_file, affine=affine, headerlines=lut_header, save_lut=True, save_tsv=True)
                            del chim_parc
            else:
                out_vol_name = fullid + '_dseg.nii.gz'

                chim_parc_name = cltbids._insert_entity(out_vol_name, {"atlas": "chimera" + chim_code} ) 
                chim_parc_file = os.path.join(str(chim_dir), chim_parc_name)
                chim_parc_lut  = os.path.join(str(chim_dir), cltbids._replace_entity_value(chim_parc_name, {"extension": "lut"}))
                chim_parc_tsv  = os.path.join(str(chim_dir), cltbids._replace_entity_value(chim_parc_name, {"extension": "tsv"}))
                
                if not os.path.isfile(chim_parc_file) or not os.path.isfile(chim_parc_lut) or not os.path.isfile(chim_parc_tsv) or force:
                    part_header = ['# $Id: {} {} \n'.format(chim_parc_lut, date_time)]
                    part_header.append('# Corresponding parcellation: {} \n'.format(chim_parc_file))
                        
                    lut_header = part_header + glob_header_info
                    lut_header = lut_header + ['\n']
                    lut_header.append('{:<4} {:<50} {:>3} {:>3} {:>3} {:>3}'.format("#No.", "Label Name:", "R", "G", "B", "A"))
                    
                    # Creating the joined parcellation
                    ref_image = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)
                    chim_parc  = cltparc.Parcellation(parc_file=ref_image, affine=affine)
                
                    # Adding the right non-cortical parcellation to the final image
                    if "rh_parc" in locals():
                        rh_parc._rearange_parc()
                        chim_parc._add_parcellation(rh_parc, append=True)
                    
                    # Adding the left non-cortical parcellation to the final image
                    if "lh_parc" in locals():
                        lh_parc._rearange_parc()
                        chim_parc._add_parcellation(lh_parc, append=True)
                    
                    # Adding the regions that do not belong to any hemisphere to the final image
                    if "mid_parc" in locals():
                        mid_parc._rearange_parc()
                        chim_parc._add_parcellation(mid_parc, append=True)
            
                    # Saving the FINAL parcellation
                    chim_parc._save_parcellation(out_file=chim_parc_file, affine=affine, headerlines=lut_header, save_lut=True, save_tsv=True)
                    del chim_parc

# Loading the JSON file containing the available parcellations
def _pipeline_info(pipe_json:str=None):
    """
    Load the JSON file containing the pipeline configuration.
    
    Parameters:
    ----------
    pipe_json : str
        JSON file containing the pipeline configuration dictionary.

    Returns:
    --------
    pipe_dict : dict
        Dictionary containing the pipeline information

    """
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Get the absolute of this file
    if pipe_json is None:
        
        pipe_json = os.path.join(cwd, 'config', 'pipe_config.json')
    else:
        if not os.path.isfile(pipe_json):
            raise ValueError("Please, provide a valid JSON file containing the pipeline configuration dictionary.")
    
    with open(pipe_json) as f:
        pipe_dict = json.load(f)
    
    return pipe_dict

def _set_templateflow_home(tflow_home: str='local'):
    """
    Setting up the templateflow home directory.
    
    Parameters:
    ----------
    tflow_home : str or Path
        Templatefow home directory.

    Returns:
    --------
    updated_tflow_home : str
        Updated templateflow home directory.

    """
    
    orig_tflow_home = os.getenv("TEMPLATEFLOW_HOME")
    if tflow_home != 'local':
        tflow_home = Path(tflow_home)
        tflow_home.mkdir(parents=True, exist_ok=True)

        if orig_tflow_home is None:
            if tflow_home.is_dir():
                updated_tflow_home = str(tflow_home)
                os.environ["TEMPLATEFLOW_HOME"] = updated_tflow_home
                tflow.update()
        else:   
            if not os.path.samefile(orig_tflow_home, tflow_home):
                if tflow_home.is_dir():
                    updated_tflow_home = str(tflow_home)
                    os.environ["TEMPLATEFLOW_HOME"] = updated_tflow_home
                    tflow.update()
                else:
                    updated_tflow_home = orig_tflow_home
    else:
        updated_tflow_home = orig_tflow_home
            
    return updated_tflow_home    

def _mix_side_prop(st_dict: dict, boolsort: bool = True):
    """
    Method to mix all the properties for a specific supra-region and a specific method.
    
    Parameters:
    ----------
    st_dict : dict
        Dictionary containing the properties separated by sides.

    Returns:
    --------
    pipe_dict : dict
        Dictionary containing the pipeline information

    """
    
    all_index = []
    all_name = []
    all_color = []
    for side in st_dict.keys():
        all_index = all_index + st_dict[side]["index"]
        all_name = all_name + st_dict[side]["name"]
        all_color = all_color + st_dict[side]["color"]
    
    if boolsort:
        # Sort the all_index and apply the order to all_name and all_color
        sort_index = np.argsort(all_index)
        all_index = [all_index[i] for i in sort_index]
        all_name = [all_name[i] for i in sort_index]
        all_color = [all_color[i] for i in sort_index]
        
    return all_index, all_name, all_color


# Loading the JSON file containing the available parcellations
def _load_parcellations_info(parc_json:str=None, supra_folder:str=None):
    """
    Load the JSON file containing the available parcellations
    
    Parameters:
    ----------
    parc_json : str
        JSON file containing the parcellation dictionary.

    supra_folder : str
        Folder containing the supraregions TSV files.

    Returns:
    --------
    parc_dict : dict
        Dictionary containing the parcellation information

    supra_dict : dict
        Dictionary containing the supraregions information

    """
    chim_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Get the absolute of this file
    if parc_json is None:
        parc_json = os.path.join(chim_dir, 'config', 'supraregions_dictionary.json')
    else:
        if not os.path.isfile(parc_json):
            raise ValueError("Please, provide a valid JSON file containing the parcellation dictionary.")
    
    with open(parc_json) as f:
        parc_dict = json.load(f)
        
    # Read all the tsv files 
    if supra_folder is None:
        spr_files = glob(os.path.join(chim_dir, 'config', 'supraregions', '*.tsv'))
    else:
        if os.path.isdir(supra_folder):
            spr_files = glob(os.path.join(supra_folder, '*.tsv'))
        else:
            raise ValueError("Please, provide a valid folder containing the supraregions TSV files.")
        
    supra_dict = {}
    for spr in spr_files:
        spr_name = os.path.basename(spr).split('.')[0]
        temp_df = pd.read_csv(spr, sep='\t')
        sp_ids = temp_df["supraregion"].unique().tolist()
        
        # Create a dictionary for each supraregion
        supra_dict[spr_name] = {}
        
        sint_dict = {} # Create a dictionary for each supraregion
        for sid in sp_ids:
            
            # Create a sub dataframe for each supraregion
            sub_df = temp_df.loc[temp_df["supraregion"] == sid]
            
            # For each supraregion id get the methods used to parcellate it
            st_methods = temp_df.loc[temp_df["supraregion"] == sid, "method"].tolist()
            st_methods = np.unique(st_methods).tolist()
            
            method_dict = {} # Create a dictionary for each method
            for mid in st_methods:
                
                # Create a sub dataframe for each method
                sub_df2 = sub_df.loc[sub_df["method"] == mid]
                
                # Get the hemispheres
                st_hemi = sub_df2["hemi"].tolist()
                
                # Get the unique hemispheres
                st_hemi = np.unique(st_hemi).tolist()
                
                hemi_dict = {} # Create a dictionary for each hemisphere
                for hemi in st_hemi:
                    
                    # Create a sub dataframe for each hemisphere
                    sub_df3 = sub_df2.loc[sub_df2["hemi"] == hemi]
                    
                    # Get the indexes of the regions
                    indexes = sub_df3["index"].tolist()
                    
                    # Get the names of the regions
                    names = sub_df3["name"].tolist()
                    
                    # Get the colors of the regions
                    colors = sub_df3["color"].tolist()
                    
                    # Create a dictionary for each hemisphere
                    temp_dict = {"index": indexes, "name": names, "color": colors}
                    
                    # Add the dictionary to the hemi_dict
                    hemi_dict[hemi] = temp_dict
                
                # Add the dictionary to the method_dict
                method_dict[mid] = hemi_dict
            
            # Add the dictionary to the supra_dict
            sint_dict[sid] = method_dict
        
        # Add the dictionary to the supra_dict
        supra_dict[spr_name] = sint_dict
        
        
    return parc_dict, supra_dict  


def _create_extra_regions_parc(aparc:str, offset:int=5000):
    """
    Create a parcellation object containing the extra regions. These parcellations 
    are included in the Aparc+aseg image and include the regions stored in the Axiliary*.tsv files.
    These files are stored in the config/supraregions folder.
    
    Parameters:
    ----------
    aparc : str
        Aparc+aseg image obrained by FreeSurfer.

    Returns:
    --------
    extra_parc : Parcellation object
        Parcellation object containing the extra regions.

    """
    chim_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Reading the auxiliary tsv file
    extra_tsv = glob(os.path.join(chim_dir, 'config', 'supraregions', 'Auxiliary*.tsv'))
    
    # Reading the tsv file
    if extra_tsv:
        
        # Do a loop and read and concatenate all the files
        for ind, et in enumerate(extra_tsv):
            temp_df = pd.read_csv(et, sep='\t')
            # Append pandas dataframes
            if ind == 0:
                extra_df = temp_df
            else:
                extra_df = extra_df.append(temp_df, ignore_index=True)

    else:
        raise ValueError("Please, provide a valid auxiliary TSV file.")
    
    # Reading the Aparc+aseg image
    aparc_parc = cltparc.Parcellation(parc_file=aparc)
    
    # Create a dataframe for the left hemisphere
    lh_df = extra_df.loc[extra_df["hemi"] == "lh"]
    
    # Sort the dataframe according to the name
    lh_df = lh_df.sort_values(by='name')
    
    # Create a dataframe for the right hemisphere
    rh_df = extra_df.loc[extra_df["hemi"] == "rh"]
    
    # Sort the dataframe according to the name
    rh_df = rh_df.sort_values(by='name')
    
    # Create a dictionary for the structures without hemispheres
    mid_df = extra_df.loc[extra_df["hemi"] == "mid"]
    
    # Sort the dataframe according to the name
    mid_df = mid_df.sort_values(by='name')
    
    # Create the parcellation for the left hemisphere
    lh_tmp_parc = copy.deepcopy(aparc_parc)
    lh_tmp_parc._keep_by_code(codes2look=lh_df["index"].tolist())
    lh_tmp_parc.index = lh_df["index"].tolist()
    lh_tmp_parc.name  = lh_df["name"].tolist()
    lh_tmp_parc.color = lh_df["color"].tolist()
    lh_tmp_parc._adjust_values()
    lh_tmp_parc._rearange_parc()
        
    # Create the parcellation for the right hemisphere
    rh_tmp_parc = copy.deepcopy(aparc_parc)
    rh_tmp_parc._keep_by_code(codes2look=rh_df["index"].tolist())
    rh_tmp_parc.index = rh_df["index"].tolist()
    rh_tmp_parc.name  = rh_df["name"].tolist()
    rh_tmp_parc.color = rh_df["color"].tolist()
    rh_tmp_parc._adjust_values()
    rh_tmp_parc._rearange_parc()
    
    # Create the parcellation for the structures without hemispheres
    mid_tmp_parc = copy.deepcopy(aparc_parc)
    mid_tmp_parc._keep_by_code(codes2look=mid_df["index"].tolist())
    mid_tmp_parc.index = mid_df["index"].tolist()
    mid_tmp_parc.name  = mid_df["name"].tolist()
    mid_tmp_parc.color = mid_df["color"].tolist()
    mid_tmp_parc._adjust_values()
    mid_tmp_parc._rearange_parc()
    
    # Unify the parcellations
    rh_tmp_parc._add_parcellation(lh_tmp_parc, append=True)
    rh_tmp_parc._add_parcellation(mid_tmp_parc, append=True)
    rh_tmp_parc._rearange_parc(offset=offset)
    
    return rh_tmp_parc

    
def _build_args_parser():

    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=52)

    from argparse import ArgumentParser

    p = argparse.ArgumentParser(formatter_class=cltmisc.SmartFormatter, description='\n Help \n')

    requiredNamed = p.add_argument_group('Required arguments')
    requiredNamed.add_argument('--regions', '-r', action='store_true', required=False,
                                help="R| List of available parcellations for each supra-region. \n"
                                "\n")

    requiredNamed.add_argument('--bidsdir', '-b', action='store', required=False, metavar='BIDSDIR', type=str, nargs=1,
                                help="R| BIDs dataset folder. \n"
                                "\n")
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
    
    requiredNamed.add_argument('--derivdir', '-d', action='store', required=False, metavar='DERIVDIR', type=str, nargs=1,
                                help="R| Folder containing the derivatives. If the folder does not exist it will be created \n"
                                    " If this option is not supplied the derivatives folder will be created inside the BIDs folder. \n"
                                "\n",
                                default=None)
    requiredNamed.add_argument('--freesurferdir', '-fr', action='store', required=False, metavar='FREESURFERDIR', type=str, nargs=1,
                                help="R| FreeSurfer subjects dir. If the folder does not exist it will be created. \n"
                                "\n",
                                default=None)
    requiredNamed.add_argument('--scale', '-s', action='store', required=False,
                                metavar='SCALE', type=str, nargs=1,
                                help="R| Scale identification.\n"
                                    " This option should be supplied for multi-resolution cortical parcellations (e.g. Lausanne or Schaeffer). \n"
                                    " If the scale is not specified. The parcellations will be generated for all the scales. \n"
                                    "\n", default=None)
    requiredNamed.add_argument('--seg', '-e', action='store', required=False,
                                metavar='SEG', type=str, nargs=1,
                                help="R| Segmentation identifier.\n"
                                    " This option should be supplied when a cortical parcellation have different versions (e.g. Shaeffer \n"
                                    " parcellation contains two different versions for the same scale: 7n and kong7n ) \n"
                                    "\n", default=None)
    requiredNamed.add_argument('--nthreads', '-n', action='store', required=False, metavar='NTHREADS', type=str, nargs=1,
                                help="R| Number of processes to run in parallel (default= Number of cores - 4). \n", default=['auto'])

    requiredNamed.add_argument('--growwm', '-g', action='store', required=False, metavar='GROWWM', type=str, nargs=1,
                                help="R| Grow of GM labels inside the white matter in mm. \n", default=None)

    requiredNamed.add_argument('--subjids', '-ids', action='store', required=False, metavar='FILE', type=str, nargs=1,
                                help="R| Subject IDs. Multiple subject ids can be specified separating them by a comma. \n"
                                    " The ID should be the basename of the T1-weighted image that will be ran. \n"
                                    " A txt file containing the IDs can be also used. \n"
                                    " Example of this file: \n"
                                    "   sub-00001_ses-0001_run-2 \n"
                                    "   sub-00001_ses-0003_run-1\n"
                                    "   sub-00001_ses-post_acq-mprage\n"
                                    " \n", default=None)
    
    requiredNamed.add_argument('--config', '-c', action='store', required=False, metavar='PIPECONFIG', type=str, nargs=1,
                                help="R| Pipeline configuration file. \n"
                                    " \n", default=None)
    
    requiredNamed.add_argument('--mergectx', '-mctx', action='store_true', required=False,
                                help="R| Join cortical white matter and cortical gray matter regions. \n"
                                    "\n", default=False)
    
    requiredNamed.add_argument('--force', '-f', action='store_true', required=False,
                                help="R| Overwrite the results. \n"
                                    "\n")
    
    p.add_argument('--verbose', '-v', action='store', required=False,
                    type=int, nargs=1,
                    help='verbosity level: 1=low; 2=debug')

    args = p.parse_args()

    global bids_dirs, supra_dict, deriv_dirs, fssubj_dirs, parcodes, pipe_json
    
    pipe_json = args.config
    
    if isinstance(args.config, list):
        pipe_json = args.config[0]
    
    if args.regions is True: 
        if args.bidsdir is None and args.parcodes is None:
            print(' ')
            mess = "Available parcellations for each supra-region"
            print('{}{}{}{}{}: '.format(bcolors.BOLD, bcolors.PURPLE, mess, bcolors.ENDC, bcolors.ENDC))  
            _print_availab_parcels()
            sys.exit()
        
        elif args.bidsdir is None or args.parcodes is None:
            print('--bidsdir and --parcodes are REQUIRED arguments')
            sys.exit()

    bids_dirs = args.bidsdir[0].split(sep=',')
    # Remove possible empty elements
    bids_dirs = [x for x in bids_dirs if x]
    
    for bids_dir in bids_dirs:
        if not os.path.isdir(bids_dir):
            print("\n")
            print("Please, supply a valid BIDs directory.")
            print("The supplied BIDs directory does not exist: {}".format(bids_dir))
            p.print_help()
            sys.exit()
            
    if args.derivdir is None:
        print('--derivdir is not supplied. ')
        print('The derivatives directory will be created in the corresponding BIDs directory.')
        deriv_dirs = []
        for bids_dir in bids_dirs:
            print('derivatives_dir: {}'.format(os.path.join(bids_dir, 'derivatives')))
            
            deriv_dir = Path(os.path.join(bids_dir, 'derivatives'))
            deriv_dir.mkdir(parents=True, exist_ok=True)
            
            deriv_dirs.append(str(deriv_dir))

    else:
        deriv_dirs = args.derivdir[0].split(sep=',')
        # Remove possible empty elements
        deriv_dirs = [x for x in deriv_dirs if x]
        
        if len(deriv_dirs) != len(bids_dirs):
            print("\n")
            print("The number of derivatives directories should be the same as the number of BIDs directories.")
            print('The first derivatives directory will be the same for all BIDs directories:')
            print('derivatives_dir: {}'.format(deriv_dirs[0]))
            
            deriv_dir = Path(deriv_dirs[0])
            
            # Create the folder if it does not exist
            deriv_dir.mkdir(parents=True, exist_ok=True)
            
            # Create a list of the same length as the number of BIDs directories
            deriv_dirs = [str(deriv_dir) for i in range(len(bids_dirs))]
        else:
            for deriv_dir in deriv_dirs:
                deriv_dir = Path(deriv_dir)
                
                # Create the folder if it does not exist
                deriv_dir.mkdir(parents=True, exist_ok=True)
                
    if args.freesurferdir is None:
        print('--freesurferdir is not supplied. ')
        
        if 'SUBJECTS_DIR' in os.environ:
            print('The FreeSurfer subjects directory will be the same for all derivatives directories.')
            print('We will use the enviroment variable SUBJECTS_DIR.')
            print('freesurfer_dir: {}'.format(os.environ["SUBJECTS_DIR"]))
            
            fssubj_dir = Path(os.environ["SUBJECTS_DIR"])
            fssubj_dir.mkdir(parents=True, exist_ok=True)
            fssubj_dirs = [str(fssubj_dir) for i in range(len(deriv_dirs))]
            
            
        else:
            print('The FreeSurfer subjects directory will be created in the following derivatives directory:')
            print('freesurfer_dir: {}'.format(os.path.join(deriv_dirs[0], 'freesurfer')))
            
            fssubj_dirs = []
            for deriv_dir in deriv_dirs:
                print('derivatives_dir: {}'.format(os.path.join(deriv_dir, 'freesurfer')))
                
                fssubj_dir = Path(os.path.join(deriv_dir, 'freesurfer'))
                fssubj_dir.mkdir(parents=True, exist_ok=True)
                
                fssubj_dirs.append(str(fssubj_dir))

    else:
        fssubj_dirs = args.freesurferdir[0].split(sep=',')
        # Remove possible empty elements
        fssubj_dirs = [x for x in fssubj_dirs if x]
        
        if len(fssubj_dirs) != len(deriv_dirs):
            print("\n")
            print("The number of freesurfer directories should be the same as the number of derivatives directories.")
            print('The FreeSurfer subjects directory  will be the same for all derivatives directories')
            print('freesurfer_dir: {}'.format(fssubj_dirs[0]))
            
            fssubj_dir = Path(fssubj_dirs[0])
            
            # Create the folder if it does not exist
            fssubj_dir.mkdir(parents=True, exist_ok=True)
            
            # Create a list of the same length as the number of BIDs directories
            fssubj_dirs = [str(fssubj_dir) for i in range(len(deriv_dirs))]
        else:
            for fssubj_dir in fssubj_dirs:
                fssubj_dir = Path(fssubj_dir)
                
                # Create the folder if it does not exist
                fssubj_dir.mkdir(parents=True, exist_ok=True)

    
    parcodes     = args.parcodes[0].split(sep=',')
    parc_dict, supra_dict = _load_parcellations_info()
    supra_reg_names = list(parc_dict.keys())
    n_supra = len(supra_reg_names)
    
    # Remove empty elements
    parcodes = [x for x in parcodes if x]
    for i, parcode in enumerate(parcodes):
        if len(parcode) != n_supra:
            parcode = parcode.ljust(n_supra, 'N')
            parcodes[i] = parcode
        
        # Checking if the code is correct
        bool_exit = False
        
        for ord, sp in enumerate(supra_reg_names):
            if parcode[ord] not in parc_dict[sp].keys() and parcode[ord] != 'N':
                bool_exit = True
                print(f"The parcellation code for the {sp} ({parcode[ord]}) is not correct.")
                _print_availab_parcels(sp)
                print(" ")
                print(f"The {sp} structures will not me included in the final parcellation")
                print(" ")
    
    return p

def _launch_fsl_first(t1:str, 
                        first_parc:str,
                        cont_tech:str = 'local', 
                        cont_image:str = None, 
                        force=False):
    """
    This function executes the FIRST subcortical parcellation.
    
    Parameters:
    ----------
    t1 : str
        T1-weighted image filename.

    first_parc : str
        Ouput name for the resulting parcellation.
    
    cont_tech : str
        Container technology (e.g. singularity, docker or local).
    
    cont_image : str
        Container image.
    
    force : bool
        Overwrite the results.

    Returns:
    --------

    """
    fsl_outdir = os.path.dirname(first_parc)
    
    if not os.path.isfile(first_parc) or force:
        fsl_outdir = Path(fsl_outdir)
        fsl_outdir.mkdir(parents=True, exist_ok=True)

        cmd_bashargs = ['run_first_all', '-i', t1, '-o', str(fsl_outdir) + os.path.sep + 'temp']
        cmd_cont = cltmisc._generate_container_command(cmd_bashargs, cont_tech, cont_image) # Generating container command
        subprocess.run(cmd_cont, stdout=subprocess.PIPE, universal_newlines=True) # Running container command
        
        cmd_bashargs = ['mv', os.path.join(str(fsl_outdir),'temp_all_fast_firstseg.nii.gz'), first_parc]
        cmd_cont = cltmisc._generate_container_command(cmd_bashargs, cont_tech, cont_image) # Generating container command
        subprocess.run(cmd_cont, stdout=subprocess.PIPE, universal_newlines=True) # Running container command
        
        cmd_bashargs = ['rm', '-rf', 'temp*']
        cmd_cont = cltmisc._generate_container_command(cmd_bashargs, cont_tech, cont_image) # Generating container command
        subprocess.run(cmd_cont, stdout=subprocess.PIPE, universal_newlines=True) # Running container command
    
def _print_availab_parcels(reg_name=None):
    """
    Print the available parcellations for each supra-region.
    
    Parameters:
    ----------
    reg_name : str
        Supra-region name. Default is None.
        
    Returns:
    --------
    
    """

    data, supra_dict = _load_parcellations_info()

    if reg_name is None:
        supra_keys = data.keys()
        parc_help = ''
        for sup in supra_keys:
            parc_opts = data[sup]
            parc_help = '{} "{}:\n"'.format(parc_help, sup)
            print('{}{}{}{}{}: '.format(bcolors.BOLD, bcolors.DARKCYAN, sup, bcolors.ENDC, bcolors.ENDC))           

            for opts in parc_opts:
                desc = data[sup][opts]["name"]
                cita = data[sup][opts]["citation"]
                parc_help = '{} "{}: {} {}\n"'.format(parc_help, opts, desc, cita)
                print('{}     {}{}: {}{}{}{}{} {}{}{}'.format(bcolors.OKGREEN, opts, bcolors.ENDC, 
                                                        bcolors.ITALIC, bcolors.DARKWHITE, desc, bcolors.ENDC, bcolors.ENDC, 
                                                        bcolors.OKYELLOW, cita, bcolors.ENDC))           
            print('')
    else:
        parc_opts = data[reg_name]
        print(' ')
        print('{}{}{}{}{}: '.format(bcolors.BOLD, bcolors.DARKCYAN, reg_name, bcolors.ENDC, bcolors.ENDC))  
        
        for opts in parc_opts:
            desc = data[reg_name][opts]["name"]
            cita = data[reg_name][opts]["citation"]
            print('{}     {}{}: {}{}{}{}{} {}{}{}'.format(bcolors.OKGREEN, opts, bcolors.ENDC, 
                                                        bcolors.ITALIC, bcolors.DARKWHITE, desc, bcolors.ENDC, bcolors.ENDC, 
                                                        bcolors.OKYELLOW, cita, bcolors.ENDC)) 
        print('')

# simple progress indicator callback function
def progress_indicator(future):
    """
    A simple progress indicator for the concurrent futures
    :param future: future object
    
    """
    global lock, n_subj, n_comp, pb, pb1, chim_code
    # obtain the lock
    with lock:
        # update the counter
        n_comp += 1
        # report progress
        # print(f'{tasks_completed}/{n_subj} completed, {n_subj-tasks_completed} remain.')
        # pb.update(task_id=pb1, description= f'[red]Completed {n_comp}/{n_subj}', completed=n_subj)
        pb.update(task_id=pb1, description= f'[red]{chim_code}: Finished ({n_comp}/{n_subj})', completed=n_comp) 


def chimera_parcellation(bids_dir:str, 
                        deriv_dir:str,
                        fssubj_dir:str,
                        code_dict:dict, 
                        t1s2run_file:str = None, 
                        growwm:list = ['0'],
                        mixwm:bool = False, 
                        nthreads:int = 1):
    """
    Preparing chimera to build the parcellations.
    
    Parameters
    ----------
    bids_dir : str
        The directory with the input dataset formatted according to the BIDS standard.
        
    deriv_dir : str
        The directory where the output files are stored.
        
    code_dict : dict
        Dictionary containing 3 main keys: code, scale and seg
        -- "code" can be either a string or a list of strings.
        -- "scale" can be either a string or a list of strings.
        -- "seg" can be either a string or a list of strings.
        
    t1s2run_file : str
        File containing the list of T1w images to be processed.
        
    growwm : list
        List of values, in mm, to grow the GM regions inside the WM.
        
    nthreads : int
        Number of threads to be used.
        
    """

    # Declaring global variables
    global pipe_json, pipe_dict, layout, pb, pb1, n_subj, n_comp, lock, chim_code
    
    ######## -- Reading the configuration dictionary  ------------ #
    pipe_dict = _pipeline_info(pipe_json=pipe_json)
    
    if t1s2run_file is None:
        # Selecting all the T1w images for each BIDS directory
        layout = BIDSLayout(bids_dir, validate=False, derivatives= False)
        t1s = layout.get( extension=['nii.gz', 'nii'], suffix='T1w', return_type='filename')
    else:
        if os.path.exists(t1s2run_file):
            with open(t1s2run_file) as file:
                t1s2run = [line.rstrip() for line in file]
        else:
            t1s2run = t1s2run_file.split(',')
        
        t1s = []
        for id in t1s2run:
            
            if not os.path.isfile(id):
                id_ent = cltbids._str2entity(id)
                if 'ses' in id_ent.keys():
                    path_cad = os.path.join(bids_dir, 'sub-' + id_ent['sub'], 'ses-' + id_ent['ses'], 'anat')
                else:
                    path_cad = os.path.join(bids_dir, 'sub-' + id_ent['sub'], 'anat')
                
                if 'suffix' not in id_ent.keys():
                    id_ent['suffix'] = 'T1w'
                
                if 'extension' not in id_ent.keys():
                    id_ent['extension'] = 'nii.gz'
                
                t1_temp = os.path.join(path_cad, cltbids._entity2str(id_ent))
                if os.path.isfile(t1_temp):
                    t1s.append(t1_temp)         

    chim_codes = code_dict["code"]

    n_parc = len(chim_codes)
    n_subj = len(t1s)

    with Progress() as pb:
        pb2 = pb.add_task('[green]Parcellation: ', total=n_parc)

        # Loop around each parcellation
        for p, chim_code in enumerate(chim_codes):
            
            # Creating and configuring the Chimera object
            chim_obj = Chimera(parc_code=chim_code,
                                scale=code_dict["scale"], 
                                seg = code_dict["seg"])
            
            # Creating the color table
            chim_obj._create_table()
            
            # Configuring and downloading the templates
            chim_obj._prepare_templates(fssubj_dir=fssubj_dir)
            
                        
            # create a lock for the counter
            lock = Lock()

            # Completed subjects
            n_comp = 0
            failed = []
            
            pb.update(task_id=pb2, description= f'[green]Parcellation: {chim_code} ({p+1}/{n_parc})', completed=p+1)
            
            # print("Parcellation: % d"% (p+1), "of % d"% (n_parc))
            if nthreads == 1:
                
                pb1 = pb.add_task(f'[red]Processing: Subject ({1}/{n_subj}) ', total=n_subj)
                for i, t1 in enumerate(t1s):
                    # ent_dict = layout.parse_file_entities(t1)

                    t1_name = os.path.basename(t1)
                    temp = t1_name.split("_")
                    full_id = '_'.join(temp[:-1])
                    pb.update(task_id=pb1, description= f'[red]{chim_code}: {full_id} ({i+1}/{n_subj})', completed=i+1) 

                    chim_obj._build_parcellation(t1, bids_dir, deriv_dir, fssubj_dir, growwm, mixwm)
                
            else:
                start_time = time.perf_counter()
                
                # create a progress bar for the subjects
                pb1 = pb.add_task(f'[red]Processing: Subject ({1}/{n_subj}) ', total=n_subj)

                # Adjusting the number of threads to the number of subjects
                if n_subj < nthreads:
                    nthreads = n_subj
                    
                # start the thread pool
                with ThreadPoolExecutor(nthreads) as executor:
                    # send in the tasks
                    # futures = [executor.submit(_build_parcellation, t1s[i],
                    # bids_dir, deriv_dir, parccode, growwm, mixwm) for i in range(n_subj)]
                    
                    futures = [executor.submit(chim_obj._build_parcellation, t1s[i], bids_dir, deriv_dir, fssubj_dir, growwm, mixwm) for i in range(n_subj)]
                    
                    # register the progress indicator callback
                    for future in futures:
                        future.add_done_callback(progress_indicator)
                    # wait for all tasks to complete


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

    global bids_dirs, deriv_dirs, fssubj_dirs, parcodes
    
        
    if args.scale is not None:
        scale_id = args.scale[0].split(sep=',')
        
        # Remove empty elements
        scale_id = [x for x in scale_id if x]
        
    else: 
        scale_id = None
        
    if args.seg is not None:
        seg_id = args.seg[0].split(sep=',')
        
        # Remove empty elements
        seg_id = [x for x in seg_id if x]
    else:
        seg_id = None
        
    # Create dictionary with the code, scale and segmentation
    code_dict = {"code": parcodes, "scale": scale_id, "seg": seg_id}
        
    
    if args.subjids is not None:
        t1s2run_file = args.subjids[0]
    else:
        t1s2run_file = None
    
    if args.growwm is not None:
        growwm = args.growwm[0].split(sep=',')
        
        # Remove empty elements
        growwm = [x for x in growwm if x]
    else:
        growwm = ['0']
    
    mixwm = args.mergectx
        
    # Detecting the number of cores to be used
    ncores = os.cpu_count()
    nthreads = args.nthreads[0]
    if nthreads == 'auto':
        nthreads = ncores
        if nthreads > 4:
            nthreads = nthreads - 4
        else:
            nthreads = 1
    else:
        nthreads     = int(args.nthreads[0])
        
    for i, bids_dir in enumerate(bids_dirs):
        
        deriv_dir = deriv_dirs[i]
        fssubj_dir = fssubj_dirs[i]
        chimera_parcellation(bids_dir, 
                        deriv_dir,
                        fssubj_dir,
                        code_dict, 
                        t1s2run_file, 
                        growwm,
                        mixwm, 
                        nthreads)
            

if __name__ == "__main__":
    main()
