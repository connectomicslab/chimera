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
from rich.progress import Progress


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

        cwd = os.path.dirname(os.path.abspath(__file__))
        chim_dir = os.path.dirname(cwd)

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
                    meth_dict['Parcels'] = parcel_names
            
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
        
        cltfree._create_fsaverage_links(fssubj_dir, fsavg_dir=None, refsubj_name=self.parc_dict["Cortical"]["reference"])
        
        
        # Detecting the base directory
        cwd = os.path.dirname(os.path.abspath(__file__))
        chim_dir = os.path.dirname(cwd)
        
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
                    
                    if not ctx_parc_lh_annot or not ctx_parc_rh_annot:
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
                    ref_img = tflow.get(atlas_ref, desc=None, resolution=1, suffix='T1w', extension='nii.gz')
                    
                    # Getting the thalamic nuclei spams 
                    parc_img = tflow.get(atlas_ref, desc=None, resolution=1,atlas=atlas_str, suffix=atlas_type, extension='nii.gz')
                    
                if supra in self.supra_dict.keys():
                    meth_dict = self.parc_dict[supra]
                    st_dict = self.supra_dict[supra][supra][meth_dict["code"]]
                    if len(st_dict) == 1:
                            bs_noctx_codes = bs_noctx_codes + st_dict['none']['index']
                            bs_noctx_names = bs_noctx_names + st_dict['none']['name']
                            bs_noctx_colors = bs_noctx_colors + st_dict['none']['color']
                            
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
    
    def _build_parcellation1(self, t1:str, bids_dir:str, 
                            deriv_dir:str = None,
                            fssubj_dir:str = None,
                            growwm:Union[str, int] = None, 
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
        global pipe_dict, layout
        
        # Detecting the base directory
        cwd = os.path.dirname(os.path.abspath(__file__))
        chim_dir = os.path.dirname(cwd)
        
        if not os.path.isfile(t1):
            raise ValueError("Please provide a valid T1 image.")

        # Detecting the entities of the T1 image
        ent_dict = layout.parse_file_entities(t1)
        
        # Getting the entities from the name 
        anat_dir = os.path.dirname(t1)
        t1_name = os.path.basename(t1)
        temp_entities = t1_name.split('_')[:-1]
        fullid = "_".join(temp_entities)
        
        if 'session' in ent_dict.keys():
            path_cad       = "sub-" + ent_dict["subject"] + os.path.sep + "ses-" + ent_dict["session"]
        else:
            path_cad       = "sub-" + ent_dict["subject"]
            
        # Creating Chimera directories
        if deriv_dir is None:
            chim_dir = os.path.join(bids_dir, 'derivatives', 'chimera', path_cad)
        else:
            chim_dir = os.path.join(deriv_dir,'chimera', path_cad)
        
        # Create the Chimera directory if it does not exist
        chim_dir = Path(chim_dir)
        chim_dir.mkdir(parents=True, exist_ok=True)
        
        # Create the anat directory if it does not exist
        tmp_dir = os.path.join(str(chim_dir), 'tmp')
        tmp_dir = Path(tmp_dir)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        
        ######## ----------- Detecting FREESURFER_HOME directory ------------- #
        fshome_dir = os.getenv('FREESURFER_HOME')
        fslut_file = os.path.join(fshome_dir, 'FreeSurferColorLUT.txt')
        
        ######## ------------- Reading FreeSurfer color lut table ------------ #
        st_codes, st_names, st_colors = cltparc.Parcellation.read_luttable(fslut_file)
        
        ######## ----- Running FreeSurfer if it was not previously computed ------ #
    
        sub2proc = cltfree.FreeSurferSubject(fullid, subjs_dir=fssubj_dir)
        supra_names = list(self.parc_dict.keys())
        
        cont_tech  = pipe_dict["packages"]["freesurfer"]["cont_tech"]
        cont_image = pipe_dict["packages"]["freesurfer"]["container"]
        
        if 'F' in self.parc_code:
            sub2proc._launch_freesurfer(force=force, 
                                        cont_tech=cont_tech, 
                                        cont_image=cont_image)
        bool_ctx = False
        
        # This is done to avoid the cortical parcellation in the first step. 
        # We are going to create parcellations Right, Left and Middle parcellations 
        # that will be added to the different cortical configurations
        
        if 'Cortical' in supra_names:
            # Remove it from the list
            supra_names.remove('Cortical')
            bool_ctx = True
        
        for supra in supra_names:
            
            # Getting the information of the common atributes
            atlas_code    = self.parc_dict[supra]["code"]
            atlas_str     = self.parc_dict[supra]["atlas"]
            atlas_desc    = self.parc_dict[supra]["description"]
            atlas_cita    = self.parc_dict[supra]["citation"]
            atlas_src     = self.parc_dict[supra]["source"]
            atlas_ref     = self.parc_dict[supra]["reference"]
            
            supra_dict = self.supra_dict[supra][supra][atlas_code]
            
            proc_dict = self.parc_dict[supra]["processing"]
            
            
            if proc_dict['method'] == 'comform2native':
                # Running the conformation to native space
                sub2proc._conform2native(ref_id=atlas_ref, 
                                        force=force, 
                                        cont_tech=cont_tech, 
                                        cont_image=cont_image)

        
        
            # # Selecting the source and downloading the parcellation
            # if atlas_src == 'templateflow':

            #     # Reference space
            #     ref_img = tflow.get(atlas_ref, desc=None, resolution=1, suffix='T1w', extension='nii.gz')
                
            #     # Getting the thalamic nuclei spams 
            #     parc_img = tflow.get(atlas_ref, desc=None, resolution=1,atlas=atlas_str, suffix=atlas_type, extension='nii.gz')
                
            # if supra in self.supra_dict.keys():
            #     meth_dict = self.parc_dict[supra]
            #     st_dict = self.supra_dict[supra][supra][meth_dict["code"]]
            #     if len(st_dict) == 1:
            #             bs_noctx_codes = bs_noctx_codes + st_dict['none']['index']
            #             bs_noctx_names = bs_noctx_names + st_dict['none']['name']
            #             bs_noctx_colors = bs_noctx_colors + st_dict['none']['color']
                        
            #     elif len(st_dict) == 2:
            #         lh_noctx_codes = lh_noctx_codes + st_dict['lh']['index']
            #         rh_noctx_codes = rh_noctx_codes + st_dict['rh']['index']
                    
            #         lh_noctx_names = lh_noctx_names + st_dict['lh']['name']
            #         rh_noctx_names = rh_noctx_names + st_dict['rh']['name']
                    
            #         lh_noctx_colors = lh_noctx_colors + st_dict['lh']['color']
            #         rh_noctx_colors = rh_noctx_colors + st_dict['rh']['color']
        
        
        if bool_ctx:
            
            # Atributes for the cortical parcellation
            atlas_names   = self.parc_dict["Cortical"]["parcels"]
            
            proc_dict     = self.parc_dict["Cortical"]["processing"]
            ctx_meth      = proc_dict["method"]
            
            nctx_parc     = len(self.parc_dict["Cortical"]["processing"]["labels"]["lh"])
            for c in np.arange(nctx_parc):
                
                ## -------- Cortical parcellation for the left hemisphere ---------------
                # Creating the name for the output file
                lh_in_parc = self.parc_dict["Cortical"]["processing"]["labels"]["lh"][c]
                at_name = [s for s in atlas_names if s in lh_in_parc]
                lh_out_annot = os.path.join(deriv_dir, 
                                            self.parc_dict["Cortical"]["deriv_surffold"], 
                                            path_cad, 
                                            fullid + '_hemi-L' + '_' + ''.join(at_name) + '_dseg.annot')
                
                ## -------- Cortical parcellation for the right hemisphere ---------------
                # Creating the name for the output file
                rh_in_parc = self.parc_dict["Cortical"]["processing"]["labels"]["rh"][c]
                at_name = [s for s in atlas_names if s in rh_in_parc]
                at_name = ''.join(at_name)
                rh_out_annot = os.path.join(deriv_dir, 
                                            self.parc_dict["Cortical"]["deriv_surffold"], 
                                            path_cad, 
                                            fullid + '_hemi-R' + '_' + at_name + '_dseg.annot')
                
                if ctx_meth == 'annot2indiv':
                    # Moving to individual space
                    sub2proc._annot2ind(ref_id=self.parc_dict["Cortical"]["processing"]["reference"], 
                                    hemi='lh', 
                                    fs_annot=lh_in_parc, 
                                    ind_annot=lh_out_annot, 
                                    cont_tech = cont_tech,
                                    cont_image=cont_image, 
                                    force=force)
                    
                    sub2proc._annot2ind(ref_id=self.parc_dict["Cortical"]["processing"]["reference"], 
                                    hemi='rh', 
                                    fs_annot=rh_in_parc, 
                                    ind_annot=rh_out_annot,
                                    cont_tech = cont_tech,
                                    cont_image=cont_image, 
                                    force=force)
                    
                if ctx_meth == 'gcs2indiv':
                    # Moving to individual space
                    sub2proc._gcs2ind(fs_gcs=lh_in_parc, 
                                    hemi='lh', 
                                    ind_annot=lh_out_annot, 
                                    cont_tech=cont_tech,
                                    cont_image=cont_image,
                                    force=force)
                    
                    sub2proc._gcs2ind(fs_gcs=rh_in_parc, 
                                    hemi='rh', 
                                    ind_annot=rh_out_annot, 
                                    cont_tech=cont_tech,
                                    cont_image=cont_image,
                                    force=force)
                
                # Copying to the labels folder
                temp_lh = os.path.join(sub2proc.subjs_dir, sub2proc.subj_id, 'label',  'lh.' + at_name + '.annot')
                shutil.copyfile(lh_out_annot, temp_lh)


                # Copying to the labels folder
                temp_rh = os.path.join(sub2proc.subjs_dir, sub2proc.subj_id, 'label',  'rh.' + at_name + '.annot')
                shutil.copyfile(rh_out_annot, temp_rh)
                
                ## -------- Creating the volumetric parcellation ---------------
                out_vol_dir = os.path.join(deriv_dir, self.parc_dict["Cortical"]["deriv_volfold"], path_cad)
                if growwm is not None:
                    for ngrow in np.arange(len(growwm)):
                        if growwm[ngrow] == '0':
                            out_vol_name = fullid + '_' + at_name + '_dseg.nii.gz'
                        else:
                            
                            ent_dict = cltbids._str2entity(at_name)
                            if 'desc' in ent_dict.keys():
                                ent_dict["desc"] = ent_dict["desc"] + 'grow' + str(growwm[ngrow]) + 'mm'
                                tmp_str = cltbids._entity2str(ent_dict)
                                out_vol_name = fullid + '_' + tmp_str + '_dseg.nii.gz'
                            else:
                                out_vol_name = fullid + '_' + at_name + '_desc_grow' + str(growwm[ngrow]) + 'mm_dseg.nii.gz'

                        sub2proc._surf2vol(atlas=at_name, 
                                            out_vol=os.path.join(out_vol_dir, out_vol_name), 
                                            gm_grow=growwm[ngrow], 
                                            force=force, bool_native=True, 
                                            color_table=['tsv', 'lut'])
                        
                else:
                    out_vol_name = fullid + '_' + at_name + '_dseg.nii.gz'
                    sub2proc._surf2vol(atlas=at_name, 
                                        out_vol=os.path.join(out_vol_dir, out_vol_name), 
                                        gm_grow=growwm[ngrow], 
                                        force=force, bool_native=True)


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
    cwd = os.path.dirname(cwd)
    # Get the absolute of this file
    if pipe_json is None:
        pipe_json = os.path.join(cwd, 'config', 'pipe_config.json')
    else:
        if not os.path.isfile(pipe_json):
            raise ValueError("Please, provide a valid JSON file containing the pipeline configuration dictionary.")
    
    with open(pipe_json) as f:
        pipe_dict = json.load(f)
    
    return pipe_dict

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
    cwd = os.path.dirname(os.path.abspath(__file__))
    chim_dir = os.path.dirname(cwd)
    
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

def _build_args_parser():

    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=52)

    from argparse import ArgumentParser

    p = argparse.ArgumentParser(formatter_class=cltmisc.SmartFormatter, description='\n Help \n')

    requiredNamed = p.add_argument_group('Required arguments')
    requiredNamed.add_argument('--regions', '-r', action='store_true', required=False,
                                help="R| List of available parcellations for each supra-region. \n"
                                "\n")

    requiredNamed.add_argument('--bidsdir', '-b', action='store', required=True, metavar='BIDSDIR', type=str, nargs=1,
                                help="R| BIDs dataset folder. \n"
                                "\n")
    requiredNamed.add_argument('--parcodes', '-p', action='store', required=True,
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
                                help="R| FreeSurfer subjects dir. If the folder is not supplied \n"
                                    " If the folder does not exist it will be created. \n"
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

    requiredNamed.add_argument('--txt2filter', '-x', action='store', required=False, metavar='FILE', type=str, nargs=1,
                                help="R| File containing the basename of the NIFTI images that will be ran. \n"
                                    "   This file is useful to tun Chimera, only, on certain T1s in case of multiple T1s \n"
                                    " for the same session.\n"
                                    " Example of this file: \n"
                                    "   sub-00001_ses-0001_run-2 \n"
                                    "   sub-00001_ses-0003_run-1\n"
                                    "   sub-00001_ses-post_acq-mprage\n"
                                    " \n", default=None)
    requiredNamed.add_argument('--force', '-f', action='store_true', required=False,
                                help="R| Overwrite the results. \n"
                                    "\n")
    p.add_argument('--verbose', '-v', action='store', required=False,
                    type=int, nargs=1,
                    help='verbosity level: 1=low; 2=debug')

    args = p.parse_args()

    global bids_dirs, deriv_dirs, fssubj_dirs, parcodes
    
    if args.bidsdir is None or args.parcodes is None :
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
    
    if args.regions is True:
        print('Available parcellations for each supra-region:')
        _print_availab_parcels()
        sys.exit()

    return p

def _print_availab_parcels(reg_name=None):

    data, supra_dict = _load_parcellations_info()

    if reg_name is None:
        supra_keys = data.keys()
        parc_help = ''
        for sup in supra_keys:
            parc_opts = data[sup]
            parc_help = '{} "{}:\n"'.format(parc_help, sup)
            print(sup + ':')
            for opts in parc_opts:
                desc = data[sup][opts]["Name"]
                cita = data[sup][opts]["citation"]
                parc_help = '{} "{}: {} {}\n"'.format(parc_help, opts, desc, cita)
                print('     {}: {} {}'.format(opts, desc, cita))
            print('')
    else:
        parc_opts = data[reg_name]
        print(reg_name + ':')
        for opts in parc_opts:
            desc = data[reg_name][opts]["Name"]
            cita = data[reg_name][opts]["citation"]
            print('     {}: {} {}'.format(opts, desc, cita))
        print('')

# Search the value inside a vector
def search(values, st_tolook):
    ret = []
    for v in st_tolook:
        index = values.index(v)
        ret.append(index)
    return ret


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




def _compute_abased_thal_parc(t1, vol_tparc, deriv_dir, pathcad, fullid, aseg_nii, out_str):
    """
    Compute atlas-based thalamic parcellation. 
    
    Parameters:
    ----------
    t1 : str
        T1 image file
        
    vol_tparc : str
        Output thalamic parcellation file
        
    deriv_dir : str
        Derivative directory
        
    pathcad : str
        Path to the CAD directory
        
    fullid : str
        Full subject id
        
    aseg_nii : str
        ASEG file
        
    out_str : str
        Output string
        
    Returns:
    --------
    mial_thalparc : str
        MIAL thalamic parcellation file
        
    
    """
    

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
        
        cmd_bashargs = ['antsRegistrationSyN.sh', '-d', '3', '-f', t1_temp, '-m', t1, '-t', 's',
                        '-o', defFile + 'desc-t12mni_']

        cltmisc._generate_container_command(cmd_bashargs, cont_tech, cont_image) # Generating container command
        # subprocess.run(,
        #                 stdout=subprocess.PIPE, universal_newlines=True)

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
        _conform2native(thal_mgz, vol_tparc, fssubj_dir, fullid)

        out_parc = [vol_tparc]

    elif parcid == 'amygdala' or  parcid == 'hippocampus':
        # Running Hippocampal and Amygdala parcellation
        process = subprocess.run(
            ['segmentHA_T1.sh', fullid, fssubj_dir],
            stdout=subprocess.PIPE, universal_newlines=True)

        # Moving Hippocampal and amygdala parcellation to native space
        lh_mgz = os.path.join(fssubj_dir, fullid, 'mri', 'lh.hippoAmygLabels-T1.v21.mgz')
        lh_gz = os.path.join(volatlas_dir, fullid + '_space-orig_hemi-L_desc-' + out_str + '_dseg.nii.gz')
        _conform2native(lh_mgz, lh_gz, fssubj_dir, fullid)

        rh_mgz = os.path.join(fssubj_dir, fullid, 'mri', 'rh.hippoAmygLabels-T1.v21.mgz')
        rh_gz = os.path.join(volatlas_dir, fullid + '_space-orig_hemi-R_desc-' + out_str + '_dseg.nii.gz')
        _conform2native(rh_mgz, rh_gz, fssubj_dir, fullid)
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
        _conform2native(hypo_mgz, hypo_gz, fssubj_dir, fullid)
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
        _conform2native(bs_mgz, bs_gz, fssubj_dir, fullid)
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



    ######## ------------- Detecting FreeSurfer Subjects Directory  ------------ #
    fshome_dir = os.getenv('FREESURFER_HOME')
    fssubj_dir = os.path.join(deriv_dir, 'freesurfer')
    os.environ["SUBJECTS_DIR"] = fssubj_dir

    if not os.path.isdir(fssubj_dir):
        print("The freesurfer subjects directory is not inside the derivative folder.")
        sys.exit()

    ## ================ Starting the global color lut table
    lut_lines = ['{:<4} {:<40} {:>3} {:>3} {:>3} {:>3} \n \n'.format("#No.", "Label Name:", "R", "G", "B", "A")]
    parc_desc_lines = ["# Parcellation code: " + parccode]


    ##### ========== Selecting the cortical parcellation ============== #####
    atlas_str     = parc_dict["Cortical"][parccode[0]]["atlas"]
    atlas_desc    = parc_dict["Cortical"][parccode[0]]["description"]
    atlas_cita    = parc_dict["Cortical"][parccode[0]]["citation"]
    atlas_type    = parc_dict["Cortical"][parccode[0]]["type"]
    atlas_names   = parc_dict["Cortical"][parccode[0]]["parcels"]
    atlas_surfloc = parc_dict["Cortical"][parccode[0]]["deriv_surffold"]
    atlas_volloc  = parc_dict["Cortical"][parccode[0]]["deriv_volfold"]
    surfatlas_dir = os.path.join(deriv_dir, atlas_surfloc, path_cad, 'anat')
    volatlas_dir = os.path.join(deriv_dir, atlas_volloc, path_cad, 'anat')
    parc_desc_lines.append(atlas_desc + ' ' + atlas_cita)
    

    # Finding FreeSurfer folder
    fs_indivdir = os.path.join(fssubj_dir, fullid)

    if not os.path.isfile(os.path.join(fs_indivdir,'mri', 'aparc+aseg.mgz')):
        _launch_freesurfer(t1, fssubj_dir, fullid)

    if os.path.isfile(os.path.join(fs_indivdir,'mri', 'aparc+aseg.mgz')):
        # Finding the cortical parcellations
        out_sparc = glob(surfatlas_dir + os.path.sep + fullid + '*' + atlas_str + '*.annot')

        if len(out_sparc) != len(atlas_names)*2:
            print("The selected cortical parcellation (" +  parc_dict["Cortical"][parccode[0]]["Name"] + ") is not computed.")
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
                out_annot = os.path.join(surfatlas_dir, fullid + '_hemi-L_space-orig_' + atlas + '_dseg.label.annot')
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

                        out_annot = os.path.join(surfatlas_dir, fullid + '_hemi-L_space-orig_' + atlas_str + '_dseg.label.annot')
                        subprocess.run(['cp', ind_annot, out_annot], stdout=subprocess.PIPE, universal_newlines=True)


                out_sparc.append(out_annot) # Annot in individual space (Atlases subject's directory)

                # 2. Right Hemisphere
                ind_annot = os.path.join(fssubj_dir, fullid, 'label', 'rh.' + atlas + '.annot') # Annot in individual space (freesurfer subject's directory)
                out_annot = os.path.join(surfatlas_dir, fullid + '_hemi-R_space-orig_' + atlas + '_dseg.label.annot')
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
                            atlas_str = "atlas-Desikan2006"
                        elif atlas == "aparc.a2009s":
                            atlas_str = "atlas-Destrieux2009"

                        out_annot = os.path.join(surfatlas_dir, fullid + '_hemi-R_space-orig_' + atlas_str + '_dseg.label.annot')
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
                aseg_parc = cltparc.Parcellation(parc_file=vol_cparc[0])

            else:
                print("Error: Cannot create the parcellation because there are missing files.\n")
                sys.exit(1)

        if 'aseg_parc' in locals():
            dim = aseg_parc.data.shape
        else:
            aseg_parc = cltparc.Parcellation(parc_file=vol_cparc[0])
            dim = aseg_parc.data.shape.shape

        outparc_lh = np.zeros((dim[0], dim[1], dim[2]), dtype='uint16')  # Temporal parcellation for the left hemisphere
        outparc_rh = np.zeros((dim[0], dim[1], dim[2]), dtype='uint16')  # Temporal parcellation for the right hemisphere

        # Loading FIRST parcellation
        if parccode[1] == 'R' or parccode[2] == 'R' or parccode[3] == 'R' or parccode[4] == 'R' or parccode[7] == 'R':
            atlas_str = parc_dict["Subcortical"]["R"]["atlas"]
            atlas_desc = parc_dict["Subcortical"]["R"]["description"]
            atlas_volloc = parc_dict["Subcortical"]["R"]["deriv_volfold"]
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

            first_parc = cltparc.Parcellation(parc_file=first_parc)
            first_parc = aseg_parc.get_fdata()


        
        
        


        ##### ========== Selecting Subcortical parcellation ============== #####
        if parccode[1] in parc_dict["Subcortical"].keys():
            atlas_str = parc_dict["Subcortical"][parccode[1]]["atlas"]
            atlas_desc = parc_dict["Subcortical"][parccode[1]]["description"]
            atlas_cita = parc_dict["Subcortical"][parccode[1]]["citation"]
            atlas_volloc = parc_dict["Subcortical"][parccode[1]]["deriv_volfold"]
            volatlas_dir = os.path.join(deriv_dir, atlas_volloc, path_cad, 'anat')
            parc_desc_lines.append(atlas_desc + ' ' + atlas_cita)
        else:
            print("Incorrect subcortical parcellation. The available parcellations are: ")
            _print_availab_parcels("Subcortical")
            sys.exit(1)

        # Get the regions index 
        subc_codesl = supra_dict["Subcortical"]["Subcortical"][parccode[1]]["lh"]["index"]
        subc_namesl = supra_dict["Subcortical"]["Subcortical"][parccode[1]]["lh"]["name"]
        subc_colorsl = supra_dict["Subcortical"]["Subcortical"][parccode[1]]["lh"]["color"]
        
        subc_codesr = supra_dict["Subcortical"]["Subcortical"][parccode[1]]["rh"]["index"]
        subc_namesr = supra_dict["Subcortical"]["Subcortical"][parccode[1]]["rh"]["name"]
        subc_colorsr = supra_dict["Subcortical"]["Subcortical"][parccode[1]]["rh"]["color"]
        
        if parccode[1] == 'F':
            
            # Creating the left hemisphere
            outparc_lh = copy.deepcopy(aseg_parc)
            outparc_lh._keep_by_code(codes2look=subc_codesl, rearrange=True)
            lh_lut = cltparc.Parcellation.write_luttable(codes=subc_codesl, names=subc_namesl, colors=subc_colorsl, headerlines=[" "])
                        
            # Creating the right hemisphere
            outparc_rh = copy.deepcopy(aseg_parc)
            outparc_rh._keep_by_code(codes2look=subc_codesr, rearrange=True)
            rh_lut = cltparc.Parcellation.write_luttable(codes=subc_codesr, names=subc_namesr, colors=subc_colorsr, headerlines=[" "])


        elif parccode[1] == 'R':  # TODO
            # Creating the left hemisphere
            outparc_lh = copy.deepcopy(first_parc)
            outparc_lh._keep_by_code(codes2look=subc_codesl, rearrange=True)
            
            # Creating the right hemisphere
            outparc_rh = copy.deepcopy(first_parc)
            outparc_rh._keep_by_code(codes2look=subc_codesr, rearrange=True)


        ##### ========== Selecting Thalamic parcellation ============== #####
        try:
            atlas_str = parc_dict["Thalamus"][parccode[2]]["atlas"]
            atlas_desc = parc_dict["Thalamus"][parccode[2]]["description"]
            atlas_cita = parc_dict["Thalamus"][parccode[2]]["citation"]
            atlas_volloc = parc_dict["Thalamus"][parccode[2]]["deriv_volfold"]
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
            atlas_str = parc_dict["Amygdala"][parccode[3]]["atlas"]
            atlas_desc = parc_dict["Amygdala"][parccode[3]]["description"]
            atlas_cita = parc_dict["Amygdala"][parccode[3]]["citation"]
            atlas_volloc = parc_dict["Amygdala"][parccode[3]]["deriv_volfold"]
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
            atlas_str = parc_dict["Hippocampus"][parccode[4]]["atlas"]
            atlas_desc = parc_dict["Hippocampus"][parccode[4]]["description"]
            atlas_cita = parc_dict["Hippocampus"][parccode[4]]["citation"]
            atlas_volloc = parc_dict["Hippocampus"][parccode[4]]["deriv_volfold"]
            volatlas_dir = os.path.join(deriv_dir, atlas_volloc, path_cad, 'anat')
            parc_desc_lines.append(atlas_desc + ' ' + atlas_cita)

            if parccode[4] == 'H':
                atlas_str_ig = parc_dict["Hippocampus"]["I"]["atlas"]
                atlas_volloc_ig = parc_dict["Hippocampus"]["I"]["deriv_volfold"]
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
            atlas_str_ig = parc_dict["Hippocampus"]["I"]["atlas"]

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
            atlas_str = parc_dict["Hypothalamus"][parccode[5]]["atlas"]
            atlas_desc = parc_dict["Hypothalamus"][parccode[5]]["description"]
            atlas_cita = parc_dict["Hypothalamus"][parccode[5]]["citation"]
            atlas_volloc = parc_dict["Hypothalamus"][parccode[5]]["deriv_volfold"]
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
            atlas_str = parc_dict["Cerebellum"][parccode[6]]["atlas"]
            atlas_desc = parc_dict["Cerebellum"][parccode[6]]["description"]
            atlas_cita = parc_dict["Cerebellum"][parccode[6]]["citation"]
            atlas_volloc = parc_dict["Cerebellum"][parccode[6]]["deriv_volfold"]
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
            atlas_str = parc_dict["Brainstem"][parccode[7]]["atlas"]
            atlas_desc = parc_dict["Brainstem"][parccode[7]]["description"]
            atlas_cita = parc_dict["Brainstem"][parccode[7]]["citation"]
            atlas_volloc = parc_dict["Brainstem"][parccode[7]]["deriv_volfold"]
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
            atlas_str = parc_dict["GyralWM"][parccode[8]]["atlas"]
            atlas_desc = parc_dict["GyralWM"][parccode[8]]["description"]
            atlas_cita = parc_dict["Cortical"][parccode[8]]["citation"]
            atlas_volloc = parc_dict["GyralWM"][parccode[8]]["deriv_volfold"]
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

def test(name):
    
    
    # Create folders with the name of the task
    os.system('cp ' + name + ' /home/yaleman/Test/')
    time.sleep(1)


def code2table(code:str, 
            lut_file:str=None, boolsave:bool=False):
    """
    Extract the LUT information from the LUT file.
    
    Parameters
    ----------
    code : str
        The code of the parcellation.

    lut_file : str
        The LUT file.
        
    Returns
    -------
    st_codes_lut : list
        List of codes.

    st_names_lut : list
        List of names.
        
    st_colors_lut : list
        List of colors.
        
    """
    
    # Reading the parcellations dictionary
    parc_dict, supra_dict = _load_parcellations_info()
    

    ######## ------------ Selecting the templates  ------------ #
    if pipe_dict["templates"]["reference"]["tool"] == "templateflow":
        ###########################################################
        ################ Downloading the templates ################
        ###########################################################
        
        # Setting templateflow home directory
        tflow_dir = pipe_dict["packages"]["templateflow"]["home_dir"]
        
        if tflow_dir == "local":
            
            tflow_dir = os.environ.get('TEMPLATEFLOW_HOME')
            
            if tflow_dir is None:
                # Setting the templateflow home directory in the same directory as the script
                temp_dir = os.path.dirname(os.path.realpath(__file__))
                
                # Select the directory before 
                temp_dir = os.path.dirname(temp_dir)
                tflow_dir = os.path.join(temp_dir, 'templateflow')
        
        # Create the directory if it does not exist using the library Path
        tflow_dir = Path(tflow_dir)
        
        # If the directory does not exist create the directory and if it fails because it does not have write access send an error
        try:
            tflow_dir.mkdir(parents=True, exist_ok=True)
        except PermissionError:
            print("The TemplateFlow directory does not have write access.")
            sys.exit()
            
        if os.path.isdir(tflow_dir):        
            # Setting the templateflow home directory
            os.environ["TEMPLATEFLOW_HOME"] = tflow_dir.as_posix()
            reload(api)
            
        else:
            print("The TemplateFlow directory does not exist.")
            sys.exit()
            
        # Getting the templates
        # Reference space
        temp_cad = pipe_dict["templates"]["reference"]["space"]
        t1_temp = tflow.get(temp_cad, desc=None, resolution=1, suffix='T1w', extension='nii.gz')
        
        # Getting the thalamic nuclei spams 
        atlas_cad = pipe_dict["templates"]["spams"]["atlas"]
        thal_spam = tflow.get(temp_cad, desc=None, resolution=1,atlas=atlas_cad, suffix='probseg', extension='nii.gz')
        
    else:
        t1_temp = pipe_dict["templates"]["reference"]["space"]
        if not os.path.isfile(t1_temp):
            print("The template file does not exist.")
            sys.exit()
        
        thal_spam = pipe_dict["templates"]["spams"]["atlas"]
        if not os.path.isfile(thal_spam):
            print("The thalamic atlas file does not exist.")
            sys.exit()
        
        temp_cad = "CustomSpace"
        atlas_cad = "CustomParc"
    
    

def chimera_parcellation(bids_dir:str, 
                        deriv_dir:str,
                        fssubj_dir:str,
                        code_dict:dict, 
                        t1s2run_file:str = None, 
                        growwm:list = [], 
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
    global pipe_dict, layout, pb, pb1, n_subj, n_comp, lock, chim_code
    
    ######## -- Reading the configuration dictionary  ------------ #
    pipe_dict = _pipeline_info()
    
    # Selecting all the T1w images for each BIDS directory
    layout = BIDSLayout(bids_dir, validate=False)
    t1s = layout.get( extension=['nii.gz', 'nii'], suffix='T1w', return_type='filename')

    # Filtering the T1w images to be processed
    if os.path.isfile(t1s2run_file):
        t1s = cltmisc._select_ids_from_file(t1s, t1s2run_file)
    else:
        t12run = t1s2run_file.split(',')
        t1s = [s for s in t1s if any(xs in s for xs in t12run)]
    
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
                
                pb1 = pb.add_task(f'[red]Processing: Parcellation {chim_code} ({p + 1}/{n_parc}) ', total=n_subj)
                for i, t1 in enumerate(t1s):
                    # ent_dict = layout.parse_file_entities(t1)
                    
                    t1_name = os.path.basename(t1)
                    temp = t1_name.split("_")
                    full_id = '_'.join(temp[:-1])
                    chim_obj._build_parcellation1(t1, bids_dir, deriv_dir, fssubj_dir, growwm)
                    pb.update(task_id=pb1, description= f'[red]{chim_code}: {full_id} ({i+1}/{n_subj})', completed=i+1) 
                
            else:
                start_time = time.perf_counter()
                
                # create a progress bar for the subjects
                pb1 = pb.add_task(f'[red]Processing: Parcellation {chim_code} ({p + 1}/{n_parc}) ', total=n_subj)

                # Adjusting the number of threads to the number of subjects
                if n_subj < nthreads:
                    nthreads = n_subj
                    
                # start the thread pool
                with ThreadPoolExecutor(nthreads) as executor:
                    # send in the tasks
                    # futures = [executor.submit(_build_parcellation, t1s[i],
                    # bids_dir, deriv_dir, parccode, growwm) for i in range(n_subj)]
                    
                    futures = [executor.submit(test, t1s[i]) for i in range(n_subj)]
                    
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
        
    
    if args.txt2filter is not None:
        t1s2run_file = args.txt2filter[0]
    else:
        t1s2run_file = ''
    
    if args.growwm is not None:
        growwm = args.growwm[0].split(sep=',')
        
        # Remove empty elements
        growwm = [x for x in growwm if x]
    else:
        growwm = ['0']
        
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
                        nthreads)
            

if __name__ == "__main__":
    main()
