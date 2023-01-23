# __CHIMERA__: An open source framework for combining multiple parcellations

Creating multi-source parcellations of the human brain is a fundamental task at several steps of the MRI analysis research workflow. __Chimera__ facilitates this otherwise difficult operation with an intuitive and flexible interface for humans and machines, thereby assisting in the construction of sophisticated and more reliable processing pipelines.
This repository contains the source code and atlases needed by __Chimera__.

##### Parcellations fusion
Chimera defines nine different supra-regions (cortex, basal ganglia, thalamus, amygdala, hippocampus, hypothalamus, cerebellum, brainstem and white-matter). Basal ganglia includes only the regions that are not labeled as supra-regions. Subdivisions in each supra-region will be populated with the parcellation
information of a single source. The available parcellation sources per supra-region, as well as one corresponding parcellation name, and a one-character unique identifier are configured in a JSON (JavaScript Object Notation) file. 
__Chimera code__: A sequence of nine one-character identifiers (one per each supra-region) unambiguosly denotes a single instance of combined parcellation (Figure. 1B). Given the sequence of nine identifier characters, Chimera selects the atlas and/or applies the corresponding methodology to obtain the parcellation for each supra-region. These supra-region-specific parcellations are finally integrated to obtain the combined volumetric parcellation for each input subject, as well as its corresponding tab-separated values table of labels, region names, and rendering colors for visualization.
Chimera uses FreeSurfer to map cortical templates from fsaverage to individual space. It also applies different methods to obtain the hippocampal subfields and brainstem parcellations as well as the thalamic, amygdala and hypothalamic nuclei segmentations. FIRST and ANTs are also used for segmenting subcortical structures and thalamic nuclei respectively.
![](Figure1.png)

##### Requirements
Required python packages: 
- [pandas], [pybids],  [numpy], [nibabel], [time], [os], [pathlib], [argparse], [sys], [csv], [subprocess]

Required image processing packages: 
- [FreeSurfer (version>7.2.0)], [FSL],  [ANTs] 
---

##### Options:
Brief description of input options:

| Option | Description |
| ---------- | ------ |
| `--regions`, `-r` | List available parcellations for each supra-region.|
| `--bidsdir`, `-b` | BIDs dataset folder. |
| `--derivdir`, `-d` | Derivatives folder. |
| `--parcodes`, `-p` | Sequence of nine one-character identifiers (one per each supra-region). |
| `--growwm`, `-g` | Grow of GM labels inside the white matter (mm). |
| `--t1file`, `-t` | File containing basename of T1w images that will be ran. |
| `--force`, `-f` | Overwrite the results. |
| `--verbose`, `-v` | Verbose (**0**, **1** or **2**). |
| `--help`, `-h` | Help. |

***
##### Usage
General command line to use **Chimera**:
```sh
    $ python chimera_parcellation.py -b <BIDs directory> -d <Derivatives directory> -p <Chimera code>
```

##### Simple examples

1. Running **Chimera** for 3 different parcellation codes (LFMFIIFIF,SFMFIIFIF,CFMFIIFIF). This will obtain the combined parcellations for all the T1-weighted images inside the BIDs dataset.

```sh
    $ python chimera_parcellation.py -b <BIDs directory> -d <Derivatives directory> -p LFMFIIFIF,SFMFIIFIF,CFMFIIFI
```
2. Running **Chimera** for T1-weighted images included in a txt file: 

```sh
    $ python chimera_parcellation.py -b <BIDs directory> -d <Derivatives directory> -p LFMFIIFIF -t <t1s.txt>
```
Example of **t1s.txt** file
|   sub-00001_ses-0001_run-2
|   sub-00001_ses-0003_run-1
|   sub-00001_ses-post_acq-mprage

3. Cortical volumes will grow 0 and 2 mm respectively inside the white matter for the selected cortical parcellations. 
```sh
    $ python chimera_parcellation.py -b <BIDs directory> -d <Derivatives directory> -p LFMFIIFIF -g 0,2
```

## Main files in the repository
1. __chimera_parcellation.py__: Main python library for performing **Chimera** parcellations. 
2. __parcTypes.json__: JSON file especifying the available parcellation sources per supra-region.
3. **ANNOT_atlases** and **GCS_atlases**: Folder containing cortical atlases in *.annot* and *.gcs* file formats.
4. **mni_icbm152_t1_tal_nlin_asym_09c**: Folder containing the reference atlas used by the MIAL atlas-based thalamic parcellation method. The atlas is referenced in standard [MNI] (_Montreal Neurological Institute_) space with a high resolution T1 weighted image (__ICBM 2009c Nonlinear Asymmetric__ ). 
5. **thalamic_nuclei_MIALatlas**: Folder containing the spatial probabilistic maps of 14 thalamic nuclei.

#### Parcellations and methodologies for each supra-region
#### 1. Cortical parcellation
| Code | Citation |Code | Citation |
| ------ | ------ |----------- | ---------- |
|`D` | Desikan et al, 2006 | `X` | Destrieux et al, 2009 |
|`T` | Klein and Tourville, 2012|`B` | Fan et al, 2016 |
|`R` | Broadmann, 1909 |`C` | Campbell, 1905 |
|`K` | Kleist, 1934 |`L` | Symmetric version of Cammoun et al, 2012 |
|`H` | Glasser et al, 2016 |`S` | Schaefer et al, 2018 |
|`M` | Smith et al, 1907 |`V` | von Economo and Koskinas, 1925 |
|`Y` | Yeo et al, 2011 |`F` | Flechsig, 1920 |

#### 2. Basal ganglia parcellation
| Code | Citation |Code | Citation |
| ------ | ------ |----------- | ---------- |
|`F` | Fischl et al, 2002 |`R`| Patenaude et al, 2011 |

#### 3. Thalamus parcellation
| Code | Citation |Code | Citation |
| ------ | ------ |----------- | ---------- |
| `F` | Fischl et al, 2002 |`I`|  Iglesias et al, 2018|
|`M` | Najdenovska and Alemán-Gómez et al, 2018 |`R`| Patenaude et al, 2011|

#### 4. Amygdala parcellation
| Code | Citation |Code | Citation |
| ------ | ------ |----------- | ---------- |
|`F`| Fischl et al, 2002|`I`| Saygin et al, 2017|
|`R`| Patenaude et al, 2011

#### 5. Hippocampus parcellation
| Code | Citation |Code | Citation |
| ------ | ------ |----------- | ---------- |
|`F`| Fischl et al, 2002|`I`| Iglesias et al, 2015|
|`H`| Iglesias et al, 2015|`I`|Patenaude et al, 2011|

#### 6. Hypothalamus parcellation 
| Code | Citation |Code | Citation |
| ------ | ------ |----------- | ---------- |
|`F`| Based on in-house protocol|`I`| Billot et al, 2020|

#### 7. Cerebellum parcellation 
| Code | Citation |
| ------ | ------ |
|`F`| Fischl et al, 2002|

#### 8. Brainstem parcellation 
| Code | Citation |Code | Citation |
| ------ | ------ |----------- | ---------- |
|`F`| Fischl et al, 2002|`I`| Iglesias et al, 2015|
|`R`| Patenaude et al, 2011|

#### 9. Gyral white matter parcellation 
| Code | Citation |
| ------ | ------ |
|`F`| Cortical (Depends on the cortical parcellation)|

##### Results

## License
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

[FreeSurfer (version>7.2.0)]: <https://surfer.nmr.mgh.harvard.edu/>
[FSL]: <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki>
[ANTs]: <http://stnava.github.io/ANTs/>
   [Nifti-1]: <https://www.nitrc.org/docman/view.php/26/204/TheNIfTI1Format2004.pdf>
   [MNI]: <https://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152NLin2009> 
   [subprocess]: <https://docs.python.org/3/library/subprocess.html>
   [numpy]:<https://numpy.org/>
   [nibabel]:<https://nipy.org/nibabel/>
   [time]:<https://docs.python.org/3/library/time.html>
   [os]:<https://docs.python.org/3/library/os.html>
   [pathlib]:<https://docs.python.org/3/library/pathlib.html>
   [argparse]:<https://docs.python.org/3/library/argparse.html>
   [sys]:<https://docs.python.org/3/library/sys.html>
   [csv]:<https://docs.python.org/3/library/csv.html>
   [pybids]:<https://bids-standard.github.io/pybids/>
   [pandas]:<https://pandas.pydata.org/>
   
 
