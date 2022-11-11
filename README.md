# Dynamic Causal Modelling (DCM) of fMRI & M/EEG
This repository consists of demo scripts for Dynamic Causal Modelling (DCM) of M/EEG data from the Wakeman & Henson (2015) open dataset. Connectivity of the face-processing network is studied using a single subject's DCM, as well as group DCM on 16 subjects using a hierarchical Bayesian framework called Parametric Empirical Bayes (PEB). SPM scripts in MATLAB and materials in this repository were used for demonstration during [COGNESTIC (Cognitive Neuroscience Skills Training in Cambridge)](https://imaging.mrc-cbu.cam.ac.uk/methods/COGNESTIC2022) held at the MRC Cognition and Brain Sciences Unit (CBU), University of Cambridge in September 2022.

## Setup/How-to
In order to get everything up and running, you need to do two things: prepare your environment, and prepare data.

### Environment
1. Ensure you have a fully functional copy of MATLAB installed on your PC/cluster. Please use a recent version of MATLAB, SPM supports all versions up to R2021b. All scripts in this repository were tested on MATLAB R2018a.
2. [Install SPM12 by following instructions at the official website.](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) After installation, launch MATLAB and in the command window type `spm` and press the return key. This should launch the SPM GUI, confirming a successful installation.
3. Clone this repository on your system by navigating to a folder on your system and running `git clone git@github.com:pranaysy/cognestic22_multimodal_dcm.git`. Alternatively if you would not like to use git, navigate to the top of this page, click on the green 'Code' button and select 'Download ZIP'. This will download a copy of the entire repository on your system. Extract this to a folder on your system.
4. Navigate to the folder containing the cloned/downloaded repository, and add the 'code' folder to MATLAB's path recursively by running: `addpath(genpath('code'))`

### Dataset
There are two ways to proceed with the DCM analyses presented here: you could either start with the raw multimodal dataset and process it yourself, or obtain processed data ready for fitting DCMs.

#### Raw data
The raw BIDS-formatted dataset can be obtained from [OpenNeuro](https://openneuro.org/datasets/ds000117) using tools like [datalad](https://www.datalad.org/) or [git-annex.](https://git-annex.branchable.com/). This is about 85GB for both fMRI and MEG. Download the data into the `data/` folder. The raw data can now be processed for subsequent DCM analyses.

##### fMRI
Processing of raw fMRI is done as per Appendix 2 in the Supplementary of Henson et al (2019). Once processed, volumes-of-interest (VOI) need to be extracted after concatenated multiple runs per subject. The concatentation and VOI extraction can be done via the batch interface or using scripts provided.

##### M/EEG
Processing of the raw M/EEG data is done according to Henson et al (2019), and a [convenience script is provided](https://github.com/pranaysy/cognestic22_multimodal_dcm/blob/main/code/meg/spm_master_script_data_preprocessing.m) for processing data in preparation for this tutorial. The script differs from the one in Henson et al (2019) in three ways:
  1. Baseline correction is done during epoching in this tutorial, but is skipped in the paper.
  2. Robust averaging is done per condition in this tutorial, whereas a simple average is taken in the paper.
  3. Specify contrasts to create two conditions: one for combining famous and unfamiliar faces into one condition, Faces and another for Scrambled.
  4. The script for this tutorial stops at forward modelling, with an additional line of code to force-write the estimated leadfield to a file, whereas the leadfield is estimated only during inversion in the paper.

#### Processed data
The pre-processing steps can be skipped entirely and tutorial-ready data can be directly obtained from figshare.
##### fMRI
  - Option A: Smoothed normalised images, condition onsets and motion parameters for 15 subjects can be downloaded [from figshare here.](https://figshare.com/articles/dataset/fMRI_Data/20936143) Subject 10 had fewer scans than the remaining subjects, so was excluded from this data. The data need to be processed further by concatenating across runs, reparametrizing conditions, and restimating the GLM per subject followed by extraction of VOI time courses.
  - Option B: DCM-ready data with VOI time courses and concatenated SPM.mat files can be downloaded [from figshare here.](https://figshare.com/articles/dataset/Face_processing_M_EEG_data_for_Dynamic_Causal_Modelling/21333996)
##### M/EEG
Artefact-corrected grand-averaged data with forward models for all 16 subjects can be downloaded [from figshare here.](https://figshare.com/articles/dataset/Face_processing_M_EEG_data_for_Dynamic_Causal_Modelling_Faces_vs_Scrambled_/21342066)

Once processed DCM-ready data for both/either fMRI and M/EEG have been prepared, the `data/` folder should have the structure indicated in [./data/folder_structure.txt](https://github.com/pranaysy/cognestic22_multimodal_dcm/blob/main/data/filelist.txt)

### List of modified SPM functions
  1. [`spm_fmri_concatenate.m`](https://github.com/pranaysy/cognestic22_multimodal_dcm/blob/main/code/modified_spm_functions/spm_fmri_concatenate.m): Function for concatenating fMRI data across runs. Updated to return the estimated SPM struct as a variable, original version only writes this variable to the disk and does not return a variable. The variable can be used by subsequent functions in a script.
  2. [`spm_dcm_fit.m`](https://github.com/pranaysy/cognestic22_multimodal_dcm/blob/main/code/modified_spm_functions/spm_dcm_fit.m): Function for estimating single or multiple DCMs. Updated so that it parallelizes by default, unless explictly disabled. The original function does not parallelize by default, which prevents SPM's batch functions to make use of parallel processing directly. SPM's batch functionality can be used for parallel processing even with the original function, but that approach will not make use of parallel processing if the estimating is done with optimized iterative fitting with PEB using `spm_dcm_peb_fit.m`.
  3. [`spm_dcm_peb_fit.m`](https://github.com/pranaysy/cognestic22_multimodal_dcm/blob/main/code/modified_spm_functions/spm_dcm_peb_fit.m): Function for estimating multiple DCMs with optimized iterative fitting using PEB. Updated so that it parallelizes by default, unless explictly disabled. 
  4. [`spm_dcm_peb_review.m`](https://github.com/pranaysy/cognestic22_multimodal_dcm/blob/main/code/modified_spm_functions/spm_dcm_peb_review.m): Function for reviewing estimated PEB structures and averaged models (BMAs). Updated to fix a minor bug when attempting to plot the 'matrix' form of parameters. The original code relied on variables present only in an estimated DCM/GCM, and would fail when using a DCM/GCM that has not been estimated yet. 
  5. [`spm_dcm_erp.m`](https://github.com/pranaysy/cognestic22_multimodal_dcm/blob/main/code/modified_spm_functions/spm_dcm_erp.m): Function for estimating a single DCM for MEEG evoked responses. Updated to have a higher limit of 64 steps on the maximum number of iterations during the fitting process. The original function capped the number of iterations to 64 steps. The increased number of steps allows all models to converge, although all models on our data converge in fewer steps, but still greater than 64 in most cases.
  6. [`spm_dcm_erp_dipfit.m`](https://github.com/pranaysy/cognestic22_multimodal_dcm/blob/main/code/modified_spm_functions/spm_dcm_erp_dipfit.m): Function for preparing the spatial model of DCM for subsequent inversion. Updated to have 2 modes per source ROI when using the distributed 'IMG' or 'Imaging' spatial observer model. The original function defaulted to 6 modes, which we believe is not best suited for our ROIs (OFA and FFA) which are very small and close to each other, and therefore may share redundant information. We find that reducing the number of modes to 2 or 3 improves the fitted time courses at each source.
  7. [`spm_L_priors.m`](https://github.com/pranaysy/cognestic22_multimodal_dcm/blob/main/code/modified_spm_functions/spm_L_priors.m): Function for setting the priors on the spatial observer model. Updated to disable contributions from non-pyramidal populations in the Jansen-Rit-based 'ERP' neuronal model. The original function allowed for weak but non-zero contributions of the spiny stellate cell populations and the inhibitory interneuronal populations to the output of each source, which then get mapped to sensor activity by the spatial observer model. 

## References
1. Wakeman, D.G. & Henson, R.N. (2015). A multi-subject, multi-modal human neuroimaging dataset. Sci. Data 2:150001 https://doi.org/10.1038/sdata.2015.1
2. Henson RN, Abdulrahman H, Flandin G and Litvak V (2019) Multimodal Integration of M/EEG and f/MRI Data in SPM12. Front. Neurosci. 13:300. https://doi.org/10.3389/fnins.2019.00300
3. Yadav, Pranay; Henson, Rik (2022): Face processing fMRI data for Dynamic Causal Modelling. figshare. Dataset. https://doi.org/10.6084/m9.figshare.21333996.v2 
4. Yadav, Pranay; Henson, Rik (2022): Face processing M/EEG data for Dynamic Causal Modelling (Faces vs Scrambled). figshare. Dataset. https://doi.org/10.6084/m9.figshare.21342066.v1 
