#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Chimera object.

This is the refactored, decomposed counterpart of the ``Chimera`` class that
originally lived in :mod:`chimera.chimera`.  The class itself is now a thin
facade: each heavy method delegates to a dedicated module

* :mod:`chimera.template_preparation` -- ``prepare_templates``
* :mod:`chimera.region_table`         -- ``create_table`` / ``export_table`` /
                                          ``build_lut_header``
* :mod:`chimera.parcellation_builder` -- ``build_parcellation``

The ``__init__`` logic is preserved verbatim; the only addition is an optional
``pipe_dict`` that is threaded to the delegated helpers instead of the
module-level global the original relied on.
"""

import os
from typing import Union

import clabtoolkit.misctools as cltmisc

from .config_manager import load_parcellations_info
from .template_preparation import prepare_templates as _prepare_templates
from .region_table import (
    build_region_table as _build_region_table,
    export_table as _export_table,
    build_lut_header as _build_lut_header,
)
from .parcellation_builder import ParcellationBuilder


class Chimera:
    """Class to create and work with Chimera objects.

    Parameters
    ----------
    parc_code : str
        Parcellation code.

    parc_dict_file : str
        Parcellation dictionary file.

    supra_folder : str
        Folder containing the supraregions TSV files.

    pipe_dict : dict
        Pipeline configuration dictionary.  When provided it is used by the
        template-preparation and parcellation-building steps; otherwise those
        steps resolve it lazily via ``_pipeline_info``.

    Returns
    -------
    Chimera object
    """

    def __init__(
        self,
        parc_code: str,
        scale: Union[str, list] = None,
        seg: Union[str, list] = None,
        parc_dict_file: str = None,
        supra_folder: str = None,
        pipe_dict: dict = None,
    ):
        """Initialize the Chimera class.

        Notes
        -----
        The parcellation code should be a string with the following format:
        "SFMIHISIFN", where each letter corresponds to a supra-region.
        """

        # Detecting the base directory
        chim_dir = os.path.dirname(os.path.abspath(__file__))

        # Rise an error if the parcellation code is not provided
        if parc_code is None:
            raise ValueError("Please provide a parcellation code")

        if parc_dict_file is not None:
            if not os.path.isfile(parc_dict_file):
                raise ValueError("The parcellation dictionary file does not exist")
        else:
            parc_dict_file = os.path.join(
                chim_dir, "config", "supraregions_dictionary.json"
            )

        if supra_folder is not None:
            if not os.path.isdir(supra_folder):
                raise ValueError(
                    "The the folder containing the supra-regions TSV files is not valid"
                )
            else:
                self.suprafolder = supra_folder
        else:
            self.suprafolder = os.path.join(chim_dir, "config", "supraregions")

        self.parc_dict, self.supra_dict = load_parcellations_info(
            parc_json=parc_dict_file, supra_folder=supra_folder
        )

        ####  Filtering the parcellation dictionary according to the parcellation code ####
        supra_names = list(self.parc_dict.keys())

        temp_dict = {}
        for i in range(len(parc_code)):
            if parc_code[i] in self.parc_dict[supra_names[i]].keys():

                # Defining the dictionary for the method
                meth_dict = {}
                meth_dict["code"] = parc_code[i]  # Add the parcellation code

                # Append the information of the parcellation using that method
                meth_dict.update(self.parc_dict[supra_names[i]][parc_code[i]])

                # Filtering the parcellation names by the scale and segmentation
                if i == 0:
                    parcel_names = meth_dict["parcels"]

                    # Filtering the parcellation names by the scale
                    if scale is not None:
                        if isinstance(scale, list):
                            # If the scale is a list do a loop over the elements and
                            # verify if the scale contains the string '_scale-'
                            scale_tmp = []
                            for sc in scale:
                                if "_scale-" not in sc:
                                    scale_tmp.append("_scale-" + sc)

                        elif isinstance(scale, str):
                            if "_scale-" not in scale:
                                scale_tmp = "_scale-" + scale

                        # Detect if the word scale is on any of the strings in parcel_names
                        if [s for s in parcel_names if "scale" in s]:
                            parcel_names = cltmisc.filter_by_substring(
                                parcel_names, scale_tmp, bool_case=False
                            )

                    # Filtering the parcellation names by the segmentation
                    if seg is not None:
                        if isinstance(seg, list):
                            # If the seg is a list do a loop over the elements and
                            # verify if the seg contains the string '_seg-'
                            seg_tmp = []
                            for sc in seg:
                                if "_seg-" not in sc:
                                    seg_tmp.append("_seg-" + sc)

                        elif isinstance(seg, str):
                            if "_seg-" not in seg:
                                seg_tmp = "_seg-" + seg

                        parcel_names = cltmisc.filter_by_substring(
                            parcel_names, seg_tmp, bool_case=False
                        )

                    # Saving the new parcels names
                    meth_dict["parcels"] = parcel_names

                # Adding the dictionary to the temp_dict
                temp_dict[supra_names[i]] = meth_dict

            else:
                # Print a message that the parcellation code is not present in the dictionary
                print(
                    "The parcellation code {} is not present in the dictionary for the supra-region {}.".format(
                        parc_code[i], supra_names[i]
                    )
                )

        self.parc_dict = temp_dict
        self.parc_code = parc_code
        self.scale = scale
        self.seg = seg
        self.pipe_dict = pipe_dict

    def prepare_templates(self, fssubj_dir: str = None):
        """Prepare the templates for the Chimera parcellation.

        See :func:`chimera.template_preparation.prepare_templates`.
        """
        _prepare_templates(self, fssubj_dir=fssubj_dir, pipe_dict=self.pipe_dict)

    def create_table(
        self,
        wm_index_offset: int = 3000,
        reg2rem: Union[list, str, None] = None,
    ):
        """Create the table of regions produced by the Chimera parcellation.

        See :func:`chimera.region_table.build_region_table`.
        """
        _build_region_table(
            self,
            wm_index_offset=wm_index_offset,
            reg2rem=reg2rem,
            pipe_dict=self.pipe_dict,
        )

    def export_table(self, out_basename: str = None, format: Union[list, str] = "tsv"):
        """Export the table of regions to a TSV or a LUT file.

        See :func:`chimera.region_table.export_table`.
        """
        _export_table(self, out_basename=out_basename, format=format)

    def build_lut_header(self):
        """Build the header of the LUT file.

        See :func:`chimera.region_table.build_lut_header`.
        """
        return _build_lut_header(self)

    def build_parcellation(
        self,
        t1: str,
        bids_dir: str,
        deriv_dir: str = None,
        fssubj_dir: str = None,
        growwm: Union[str, int] = None,
        bool_mixwm: bool = False,
        force: bool = False,
    ):
        """Build the parcellation for the selected parcellation code.

        See :class:`chimera.parcellation_builder.ParcellationBuilder`.
        """
        ParcellationBuilder(
            self,
            t1=t1,
            bids_dir=bids_dir,
            deriv_dir=deriv_dir,
            fssubj_dir=fssubj_dir,
            growwm=growwm,
            bool_mixwm=bool_mixwm,
            force=force,
            pipe_dict=self.pipe_dict,
        ).run()
