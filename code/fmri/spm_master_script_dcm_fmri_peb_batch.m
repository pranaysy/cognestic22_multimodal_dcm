%---------------------------------------------------------------------------------------
% Group Dynamic Causal Modelling of the Face Perception Network using MEG
%---------------------------------------------------------------------------------------
% This script consists of SPM and MATLAB code for fitting Dynamic Causal Models on fMRI
% time courses. All analyses covered were presented as a tutorial at COGNESTIC-22 in
% September 2022 at the MRC Cognition and Brain Sciences Unit. The script covers
% specification of a single DCM, replication of this specified DCM to multiple subjects
% (called a GCM) and fitting this 'group' DCM in parallel. Further, various ways to
% perform inference at the group level using a hierarchical Bayesian framework called
% Parametric Empirical Bayes (PEB) are also demonstrated. These include greedy search of
% nested models, binary model comparison of nested models, and comparing families of
% nested models. Lastly, the inclusion of subject-level covariates for inference at the
% group level are also demo'ed in brief.

% Note that this script uses the 'batch' interface exposed by SPM, which in turn is
% built on MATLAB's batching functionality.

% Sections in this script are organized in the same order as covered in the tutorial.

% Authored in September 2022 by:
%   Rik Henson - rik.henson@mrc-cbu.cam.ac.uk
% With help from:
%   Pranay Yadav - pranay.yadav@mrc-cbu.cam.ac.uk

%---------------------------------------------------------------------------------------
% Data Sources
%---------------------------------------------------------------------------------------

% A. Processed DCM-ready data can be obtained from:
%   Henson, Rik (2022): fMRI VOI data. figshare. Dataset.
%   https://doi.org/10.6084/m9.figshare.21270222.v1 
% If you are starting with this data, then skip the steps 'COMBINE' and 'VOI' since the
% data are already concatenated across runs, & consist of extracted VOIs. After 'SETUP'
% jump directly to 'DCM'

% B. Processed volumes per run per subject can be obtained from:
%   Henson, Rik (2022): fMRI Data. figshare. Dataset.
%   https://doi.org/10.6084/m9.figshare.20936143.v2 
% These need to be processed further in order to fit DCMs. Steps include concatenation
% across runs, followed by extraction of time courses from Volumes-of-Interest (VOI).
% The sections 'COMBINE' and 'VOI' in this script cover these steps respectively.

% C. Alternatively, raw data can be obtained from:
%   Wakeman, D.G. & Henson, R.N. (2015). A multi-subject, multi-modal human neuroimaging
%   dataset. Sci. Data 2:150001 https://doi.org/10.1038/sdata.2015.1
% Process this data as per the tutorial instructed in Appendix 2 of:
%   Henson RN, Abdulrahman H, Flandin G and Litvak V (2019) Multimodal Integration of
%   M/EEG and f/MRI Data in SPM12. Front. Neurosci. 13:300.
%   https://doi.org/10.3389/fnins.2019.00300
% After processing, concatenate across runs and extract VOIs as described in (B) to
% obtain DCM-ready data.

%---------------------------------------------------------------------------------------
%                                                            
%      .d8888b.           888                      
%     d88P  Y88b          888                      
%     Y88b.               888                      
%      "Y888b.    .d88b.  888888 888  888 88888b.  
%         "Y88b. d8P  Y8b 888    888  888 888 "88b 
%           "888 88888888 888    888  888 888  888 
%     Y88b  d88P Y8b.     Y88b.  Y88b 888 888 d88P 
%      "Y8888P"   "Y8888   "Y888  "Y88888 88888P"  
%                                         888      
%                                         888      
%                                         888                                               
%                                                                        
%---------------------------------------------------------------------------------------
% Set up the MATLAB & SPM environment with necessary paths and variables

%---------------------------------------------------------------------------------------
% STEP 1 
%---------------------------------------------------------------------------------------
clear

% Add your local installation of SPM12 to MATLAB Path
SPM12PATH = '/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'
addpath(SPM12PATH);

%---------------------------------------------------------------------------------------
% STEP 2: Configure & launch SPM 
%---------------------------------------------------------------------------------------

% Initialize SPM
spm('asciiwelcome');
spm_jobman('initcfg'); % Allows batch operations inside a script
spm('defaults','fmri');
spm_get_defaults('cmdline',true);

spm fmri

%---------------------------------------------------------------------------------------
% STEP 3: Variables for folders
%---------------------------------------------------------------------------------------

% Specify root working directory 
base_dir = '/imaging/henson/Wakeman/multimodal_dcm';
addpath(fullfile(base_dir, 'code')) % Add scripts & functions to workspace

% Assign operational directories to variables
rawpth = fullfile(base_dir, 'data'); % contains raw BIDS data
derpth = fullfile(rawpth,'derivatives','SPM12');
% above will exist if you have run preprocessing scripts for Henson et al, 2019, 
% or else download and extract fmri_data.tar.gz from Figshare link

scrpth = fullfile(base_dir,'code');

% Add all provided scripts to MATLAB path
% Needed for parfor version of spm_dcm_peb_fit and spm_fmri_concatenate
addpath(genpath(scrpth)) 

%---------------------------------------------------------------------------------------
% STEP 4: Variables for data 
%---------------------------------------------------------------------------------------

% If you have all raw data...
% BIDS   = spm_BIDS(rawpth);
% subs   = spm_BIDS(BIDS,'subjects', 'task','facerecognition');
% runs = spm_BIDS(BIDS,'runs', 'modality','func', 'type','bold', 'task','facerecognition'); 

% ...else just re-specify
subs = compose('%02g', [1:9, 11:16]); % subject 10 had fewer scans in last run
nsub   = numel(subs)
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);
% 
% runs = {}; for r = [1:9]; runs{end+1} = sprintf('%02d',r); end
% nrun = numel(runs)

% ...alternatively, infer from folder structure

% Specify subjects, runs for VOI extraction
%subdir = dir(fullfile(derpth, 'sub*'));
%subdir = {subdir.name};
%subs = strrep(subdir, 'sub-', '');
%nsub   = numel(subs);

runs = compose('%02g', 1:9);
nrun = numel(runs);

nscan = repmat(208,1,nrun); % sub-15, run-01 has 209 scans, so ignore last
TR = 2;

%---------------------------------------------------------------------------------------
% STEP 5: Launch parallel pool
%---------------------------------------------------------------------------------------

numworkers = nsub; % Number of workers for distributed computing (depends on system)
if numworkers > 0
    delete(gcp('nocreate')) % Shut down any existing pool
    P=cbupool(numworkers, '--mem-per-cpu=4G --time=12:00:00 --nodelist=node-j10');
    parpool(P, P.NumWorkers);
%     parpool(numworkers);
end

% Proceed to the next step if you need to extract VOI files, else jump to 'DCM' section

%---------------------------------------------------------------------------------------
% 
%      .d8888b.                         888      d8b                   
%     d88P  Y88b                        888      Y8P                   
%     888    888                        888                            
%     888         .d88b.  88888b.d88b.  88888b.  888 88888b.   .d88b.  
%     888        d88""88b 888 "888 "88b 888 "88b 888 888 "88b d8P  Y8b 
%     888    888 888  888 888  888  888 888  888 888 888  888 88888888 
%     Y88b  d88P Y88..88P 888  888  888 888 d88P 888 888  888 Y8b.     
%      "Y8888P"   "Y88P"  888  888  888 88888P"  888 888  888  "Y8888
%                                                                        
%---------------------------------------------------------------------------------------
% Combine conditions and the concatenate runs for DCM

%---------------------------------------------------------------------------------------
% STEP 1: Create concatenated trial definition and movement parameter files
%---------------------------------------------------------------------------------------

% Define conditions
conds = [];
conds.names{1} = 'All'; 
conds.names{2} = 'Faces';    % Famous and Nonfamous

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comments pending
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(derpth)
for s = 1:nsub    
%     try mkdir(subdir{s}); end; 
    outdir = fullfile(derpth,subdir{s}, 'fmri')
%     try mkdir(outdir); end;
    
    for t = 1:2
        conds.durations{t} = 0;
        conds.onsets{t} = [];
    end
    
    time_so_far = 0; volfiles = {}; movepar = [];
    for r = 1:length(runs)       
        volfiles{r} = spm_select('ExtFPList',outdir,sprintf('^swsub-.*run-%s_bold\\.nii$',runs{r}),[1:nscan(r)]);
   
        trlfile = fullfile(outdir,sprintf('sub-%s_run-%s_spmdef.mat',subs{s},runs{r}));
        d = load(trlfile);
        conds.onsets{1} = [conds.onsets{1}; sort([d.onsets{1}; d.onsets{2}; d.onsets{3}]) + time_so_far];
        conds.onsets{2} = [conds.onsets{2}; sort([d.onsets{1}; d.onsets{2}]) + time_so_far];
 
        time_so_far = time_so_far + nscan(r)*TR;
        
        d = load(spm_select('FPList',outdir,sprintf('^rp.*run-%s.*\\.txt$',runs{r})));
        d = d(1:nscan(r),:);
        movepar = [movepar; d];
    end
    volsfile = fullfile(outdir,sprintf('sub-%s_run-concat_volfiles.mat',subs{s}));
    save(volsfile,'volfiles');
    condfile = fullfile(outdir,sprintf('sub-%s_run-concat_spmdef.mat',subs{s}));
    save(condfile,'-struct','conds');
    movefile = fullfile(outdir,sprintf('rp_sub-%s_run-concat_spmdef.txt',subs{s}));
    save(movefile,'movepar','-ascii');
end

% Create concatenated SPM.mat file
%----------------------------------------------------------------------

% Care: if run in parfor, will hang if SPM.mat files already exist, so make
% sure delete them before re-running
parfor (s = 1:nsub, numworkers)
    if numworkers > 0
        spm_jobman('initcfg');
        spm('defaults', 'fmri');
        spm_get_defaults('cmdline',true);
    end
    
    % Specify concatenated model
    %----------------------------------------------------------------------    
    jobfile = {fullfile(scrpth,'fmri', 'batch_stats_fmri_concatenated_specify_job.m')};

    outdir = fullfile(derpth,subdir{s}, 'fmri')
    inputs  = {};
    inputs{1} = cellstr(outdir);

    volfiles = load(fullfile(outdir,sprintf('sub-%s_run-concat_volfiles.mat',subs{s})));
    inputs{2} = cellstr(cat(1,volfiles.volfiles{:}));
    condfile = fullfile(outdir,sprintf('sub-%s_run-concat_spmdef.mat',subs{s}));
    inputs{3} = cellstr(condfile);
    movefile = fullfile(outdir,sprintf('rp_sub-%s_run-concat_spmdef.txt',subs{s}));
    inputs{4} = cellstr(movefile);
    spm_jobman('run', jobfile, inputs{:});
    
    % Call spm_fmri_concatenate to update SPM files for concatenated runs
    %----------------------------------------------------------------------    
    cd(outdir)
    spmfile = fullfile(outdir,'SPM.mat');
    SPM = spm_fmri_concatenate(spmfile, nscan); % Assumes nscan already in memory from above (could save instead)
    delete('SPM_backup.mat');
    
    % Estimate new model
    %----------------------------------------------------------------------    
    jobfile = {fullfile(scrpth, 'fmri','batch_stats_fmri_concatenated_estimate_job.m')};
    inputs  = {};
    inputs{1} = cellstr(spmfile);
    spm_jobman('run', jobfile, inputs{:});
end

%---------------------------------------------------------------------------------------
%
%     888     888  .d88888b. 8888888 
%     888     888 d88P" "Y88b  888   
%     888     888 888     888  888   
%     Y88b   d88P 888     888  888   
%      Y88b d88P  888     888  888   
%       Y88o88P   888     888  888   
%        Y888P    Y88b. .d88P  888   
%         Y8P      "Y88888P" 8888888
%                                                                        
%---------------------------------------------------------------------------------------
% Create VOI files (ie timeseries for ROIs)
%==========================================================================

ROI_names = {'lOFA','rOFA','lFFA','rFFA'};
ROI_coord = {[-38, -86, -14],[+36, -86, -10],[-42, -56, -20],[+42, -52, -14]};

parfor (s = 1:nsub, numworkers)
    if numworkers > 0
        spm_jobman('initcfg');
        spm('defaults', 'fmri');
        spm_get_defaults('cmdline',true);
    end
    
    outdir = fullfile(derpth,subdir{s}, 'fmri')
    cd(outdir)
    spmfile = fullfile(outdir,'SPM.mat');
    
    jobfile = {fullfile(scrpth,'fmri', 'batch_VOI_job.m')};  
    for r = 1:numel(ROI_names)
        inputs  = {};
        inputs{1} = cellstr(spmfile);
        inputs{2} = ROI_names{r};
        inputs{3} = ROI_coord{r};
        spm_jobman('run', jobfile, inputs{:});
    end
end

%---------------------------------------------------------------------------------------
%                                                 
%     8888888b.   .d8888b.  888b     d888 
%     888  "Y88b d88P  Y88b 8888b   d8888 
%     888    888 888    888 88888b.d88888 
%     888    888 888        888Y88888P888 
%     888    888 888        888 Y888P 888 
%     888    888 888    888 888  Y8P  888 
%     888  .d88P Y88b  d88P 888   "   888 
%     8888888P"   "Y8888P"  888       888
%
%---------------------------------------------------------------------------------------
% Specify DCM for sub-01

% Start from here if you have SPM.mat and VOI* files ready for each subject.

%---------------------------------------------------------------------------------------
% STEP 1: Data
%---------------------------------------------------------------------------------------
ref_sub = 1; % eg first subject
indir = fullfile(derpth,subdir{ref_sub}, 'fmri');
load(fullfile(indir,'SPM.mat'));
        
DCM = [];

% Run this again if VOI already present
ROI_names = {'lOFA','rOFA','lFFA','rFFA'};

for r = 1:numel(ROI_names)
    load(fullfile(indir,sprintf('VOI_%s_1.mat',ROI_names{r})),'xY');
    DCM.xY(r) = xY;
end      
      
DCM.n = length(DCM.xY);      % number of regions
DCM.v = length(DCM.xY(1).u); % number of time points

DCM.Y.dt  = SPM.xY.RT;
DCM.Y.X0  = DCM.xY(1).X0;
for i = 1:DCM.n
    DCM.Y.y(:,i)  = DCM.xY(i).u;
    DCM.Y.name{i} = DCM.xY(i).name;
end
        
DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v); % Models autocorrelation

DCM.U.dt   = SPM.Sess(1).U(1).dt;
        
for u = 1:2 
    DCM.U.u(:,u)  = SPM.Sess(1).U(u).u((32+1):end); % DCM allows for 2 TRs before first stimulus
    DCM.U.name{u} = SPM.Sess(1).U(u).name{1};
    DCM.U.idx(u,:) = [u 1]; 
end
        
DCM.delays = repmat(SPM.xY.RT,DCM.n,1)/2;
DCM.TE     = 0.03;  % 30ms on CBU

%---------------------------------------------------------------------------------------
% STEP 2: Options
%---------------------------------------------------------------------------------------

DCM.options.nonlinear  = 0;
DCM.options.two_state  = 0;
DCM.options.stochastic = 0;
DCM.options.centre     = 0; % centre regressors (default = 0)
DCM.options.nograph    = 1;
DCM.options.maxnodes   = 8;
DCM.options.maxit      = 128;
DCM.options.hidden     = [];
DCM.options.induced    = 0;
        
%DCM.M.options = struct(); % needed else crashes in estimation

%---------------------------------------------------------------------------------------
% STEP 3: Connectivity
%---------------------------------------------------------------------------------------

% Full model

DCM.a = [
     1     1     1     0
     1     1     0     1
     1     0     1     1
     0     1     1     1
];
DCM.c = [
     1     0
     1     0
     0     0
     0     0
];
DCM.b(:,:,1) = zeros(DCM.n,DCM.n); % Corresponding to 'All'
DCM.b(:,:,2) = [                   % Corresponding to 'Faces'
     1     1     1     0
     1     1     0     1
     1     0     1     1
     0     1     1     1
];
DCM.d        = zeros(DCM.n,DCM.n,0); % needed else crashes in estimation

%---------------------------------------------------------------------------------------
% STEP 4: Save
%---------------------------------------------------------------------------------------

fmri_fits = fullfile(base_dir, 'fits', 'batch_script', 'fmri');
outfile_full = fullfile(fmri_fits, 'templates', 'DCMs', 'DCM_Full.mat');
save(outfile_full,'DCM');

%---------------------------------------------------------------------------------------
% STEP 5: Estimate
%---------------------------------------------------------------------------------------

% Estimate DCM and save
DCM = spm_dcm_estimate(DCM);
save(fullfile(fmri_fits, 'DCM_Full.mat'),'DCM');

% Review and evaluate model fit
spm_dcm_review(DCM);
spm_dcm_fmri_check(DCM);

%---------------------------------------------------------------------------------------
%
%      .d8888b.                             d8b  .d888          
%     d88P  Y88b                            Y8P d88P"           
%     Y88b.                                     888             
%      "Y888b.   88888b.   .d88b.   .d8888b 888 888888 888  888 
%         "Y88b. 888 "88b d8P  Y8b d88P"    888 888    888  888 
%           "888 888  888 88888888 888      888 888    888  888 
%     Y88b  d88P 888 d88P Y8b.     Y88b.    888 888    Y88b 888 
%      "Y8888P"  88888P"   "Y8888   "Y8888P 888 888     "Y88888 
%                888                                        888 
%                888                                   Y8b d88P 
%                888                                    "Y88P"
%
%---------------------------------------------------------------------------------------
% Specify GCM (DCM for each subject and model)

%---------------------------------------------------------------------------------------
% STEP 1: Prepare job environment
%---------------------------------------------------------------------------------------

% Name this GCM (will be used as name of GCM file, with 'GCM_' appended)
GCMname = 'Full';

% Path to template folder: This is where specified GCM will be saved
templatedir = fullfile(fmri_fits, 'templates', 'GCMs', GCMname);
% Note: innermost folder name is same as GCM file name

% Path to job file for executing this operation
jobfile = {fullfile(scrpth,'fmri', 'batch_create_gcm_job.m')}; 

%---------------------------------------------------------------------------------------
% STEP 2: Specify inputs
%---------------------------------------------------------------------------------------

% Prepare inputs, 9 total, as per order listed in the job file
inputs  = {};

% Input #1: Path to output directory containing template GCM with 'Full' model
inputs{1} = cellstr(templatedir);

% Input #2: Name of GCM
inputs{2} = GCMname;

% Input #3: Path to 'Full' template DCM (specified earlier)
inputs{3} = {outfile_full};

% Input #4: Path to alternative DCMs (empty here, we'll use later for model space)
inputs{4} = ''; %{fullfile(derpth,subdir{ref_sub},'DCM_self.mat')};

% Input #5: Paths to SPM.mat file for each subject
for s = 1:nsub
    inputs{5}{s,1} = fullfile(derpth,subdir{s},'fmri','SPM.mat');
end

% Input #6-9: Paths to VOL*.mat files for each subject
for r = 1:numel(ROI_names)
    files = {};
    for s = 1:nsub
        inputs{5+r}{s,1} = fullfile(derpth,subdir{s},'fmri',sprintf('VOI_%s_1.mat',ROI_names{r}));
    end
end

%---------------------------------------------------------------------------------------
% STEP 3: Execute
%---------------------------------------------------------------------------------------

spm_jobman('run', jobfile, inputs{:});

%---------------------------------------------------------------------------------------
%                                                                                         
%     8888888888         888    d8b                        888            
%     888                888    Y8P                        888            
%     888                888                               888            
%     8888888   .d8888b  888888 888 88888b.d88b.   8888b.  888888 .d88b.  
%     888       88K      888    888 888 "888 "88b     "88b 888   d8P  Y8b 
%     888       "Y8888b. 888    888 888  888  888 .d888888 888   88888888 
%     888            X88 Y88b.  888 888  888  888 888  888 Y88b. Y8b.     
%     8888888888 88888P'  "Y888 888 888  888  888 "Y888888  "Y888 "Y8888
%
%---------------------------------------------------------------------------------------
% Estimate specified DCMs for all subjects using the batch interface (may take time)

%---------------------------------------------------------------------------------------
% STEP 1: Prepare job environment
%---------------------------------------------------------------------------------------

% Path to job file for executing this operation
jobfile = {fullfile(scrpth,'fmri','batch_fit_gcm_job.m')}; 

% Path to GCM template specified in previous step
GCMfile = fullfile(templatedir, ['GCM_' GCMname '.mat']);

%---------------------------------------------------------------------------------------
% STEP 2: Specify inputs
%---------------------------------------------------------------------------------------

% Prepare inputs, 4 total, as per order listed in the job file
inputs  = {};

% Input #1: Path to GCM template specified in previous step
inputs{1} = cellstr(GCMfile);

% Input #2: Path to folder where estimated GCM should be stored
inputs{2} = cellstr(fmri_fits);

% Input #3: Name of estimated GCM
inputs{3} = GCMname;

% Input #4: Type of estimation: 1 for single-pass, 2 for iterative optimization 
inputs{4} = 1; % 2 for iterative PEB fit

%---------------------------------------------------------------------------------------
% STEP 3: Execute
%---------------------------------------------------------------------------------------

spm_jobman('run', jobfile, inputs{:});

%---------------------------------------------------------------------------------------
% STEP 4: Review
%---------------------------------------------------------------------------------------

% Load and review estimated GCM to evaluate model fit
load(fullfile(fmri_fits,['GCM_' GCMname '.mat']))
spm_dcm_fmri_check(GCM);

%---------------------------------------------------------------------------------------
%
%     8888888b.  8888888888 888888b.   
%     888   Y88b 888        888  "88b  
%     888    888 888        888  .88P  
%     888   d88P 8888888    8888888K.  
%     8888888P"  888        888  "Y88b 
%     888        888        888    888 
%     888        888        888   d88P 
%     888        8888888888 8888888P"
%
%---------------------------------------------------------------------------------------
% Estimate PEB on fitted GCMs

%---------------------------------------------------------------------------------------
% STEP 1: Prepare job environment
%---------------------------------------------------------------------------------------

% Path to estimated GCM specified in previous step
GCMfile = fullfile(fmri_fits,['GCM_' GCMname '.mat']);

% Path to job file for executing this operation
jobfile = {fullfile(scrpth,'fmri','batch_fit_peb_job.m')}; 

%---------------------------------------------------------------------------------------
% STEP 2: Specify inputs
%---------------------------------------------------------------------------------------

% Prepare inputs, 3 total, as per order listed in the job file
inputs  = {};

% Input #1: Name of PEB: Keep same as GCM's name. 'PEB_' will be appended to filename
inputs{1} = GCMname;

% Input #2: Path to estimated GCM from previous step
inputs{2} = cellstr(GCMfile);

% Input #3: Column of index in GCM, since we only have one column, leave it as 1
inputs{3} = 1;

% Note, the job file has 'B' matrix specified, and review is turned off

%---------------------------------------------------------------------------------------
% STEP 3: Execute
%---------------------------------------------------------------------------------------
spm_jobman('run', jobfile, inputs{:});
 
%---------------------------------------------------------------------------------------
% STEP 4: Review
%---------------------------------------------------------------------------------------

% Load estimated PEB and GCM to review parameters
load(fullfile(fmri_fits,['GCM_' GCMname '.mat']))
load(fullfile(fmri_fits,['PEB_' GCMname '.mat']))
spm_dcm_peb_review(PEB, GCM)

%---------------------------------------------------------------------------------------
%
%      .d8888b.                                    888      
%     d88P  Y88b                                   888      
%     Y88b.                                        888      
%      "Y888b.    .d88b.   8888b.  888d888 .d8888b 88888b.  
%         "Y88b. d8P  Y8b     "88b 888P"  d88P"    888 "88b 
%           "888 88888888 .d888888 888    888      888  888 
%     Y88b  d88P Y8b.     888  888 888    Y88b.    888  888 
%      "Y8888P"   "Y8888  "Y888888 888     "Y8888P 888  888
%
%---------------------------------------------------------------------------------------
% Perform greedy search over full model space (B-matrix) with Bayesian model reduction

%---------------------------------------------------------------------------------------
% STEP 1: Prepare job environment
%---------------------------------------------------------------------------------------

% Path to estimated GCM and PEB from previous steps
GCMfile = fullfile(fmri_fits,['GCM_' GCMname '.mat']);
PEBfile = fullfile(fmri_fits,['PEB_' GCMname '.mat']);

% Path to job file for executing this operation
jobfile = {fullfile(scrpth,'fmri', 'batch_peb_bmr_search_job.m')}; 

%---------------------------------------------------------------------------------------
% STEP 2: Specify inputs
%---------------------------------------------------------------------------------------

% Prepare inputs, 2 total, as per order listed in the job file
inputs  = {};

% Input #1: Path to estimated PEB 
inputs{1} = cellstr(PEBfile);

% Input #2: Path to estimated GCM 
inputs{2} = cellstr(GCMfile); % This only needed for plotting - doesn't restrict to two models

%---------------------------------------------------------------------------------------
% STEP 3: Execute
%---------------------------------------------------------------------------------------
spm_jobman('run', jobfile, inputs{:});

%---------------------------------------------------------------------------------------
%
%      .d8888b.                                                           
%     d88P  Y88b                                                          
%     888    888                                                          
%     888         .d88b.  88888b.d88b.  88888b.   8888b.  888d888 .d88b.  
%     888        d88""88b 888 "888 "88b 888 "88b     "88b 888P"  d8P  Y8b 
%     888    888 888  888 888  888  888 888  888 .d888888 888    88888888 
%     Y88b  d88P Y88..88P 888  888  888 888 d88P 888  888 888    Y8b.     
%      "Y8888P"   "Y88P"  888  888  888 88888P"  "Y888888 888     "Y8888  
%                                       888                               
%                                       888                               
%                                       888
%
%---------------------------------------------------------------------------------------
% Perform Bayesian Model Comparison (BMC) for Full vs reduced Self-only models

%---------------------------------------------------------------------------------------
% STEP 1: Specify a reduced Self-only model
%---------------------------------------------------------------------------------------

% Load full model
load(outfile_full)

% Switch off between-region connections
DCM.b(:,:,2) = [
     1     0     0     0
     0     1     0     0
     0     0     1     0
     0     0     0     1
];

% Save reduced model as a template 
outfile_self = fullfile(fmri_fits, 'templates', 'DCMs','DCM_Self.mat');
save(outfile_self,'DCM');

%---------------------------------------------------------------------------------------
% STEP 2: Prepare job environment
%---------------------------------------------------------------------------------------

% Name this GCM (will be used as name of GCM file, with 'GCM_' appended)
GCMname = 'Full_vs_Self';

% Path to template folder: This is where specified GCM will be saved
templatedir = fullfile(fmri_fits, 'templates', 'GCMs', GCMname);
% Note: innermost folder name is same as GCM file name

% Path to job file for executing this operation
jobfile = {fullfile(scrpth,'fmri','batch_create_gcm_job.m')}; 

%---------------------------------------------------------------------------------------
% STEP 3: Specify inputs
%---------------------------------------------------------------------------------------

% Prepare inputs, 9 total, as per order listed in the job file
inputs  = {};

% Input #1: Path to output directory containing template GCM with 'Full' model
inputs{1} = cellstr(templatedir);

% Input #2: Name of GCM
inputs{2} = GCMname;

% Input #3: Path to 'Full' template DCM (specified earlier)
inputs{3} = {outfile_full}; %cellstr(fullfile(outpth,subdir{ref_sub},'DCM_full.mat'));

% Input #4: Path to nested 'Self'-only DCMs (can add more!)
inputs{4} = {outfile_self};%{fullfile(outpth,subdir{ref_sub},'DCM_self.mat')};

% Input #5: Paths to SPM.mat file for each subject
for s = 1:nsub
    inputs{5}{s,1} = fullfile(derpth,subdir{s},'fmri','SPM.mat');
end

% Input #6-9: Paths to VOL*.mat files for each subject
for r = 1:numel(ROI_names)
    files = {};
    for s = 1:nsub
        inputs{5+r}{s,1} = fullfile(derpth,subdir{s},'fmri',sprintf('VOI_%s_1.mat',ROI_names{r}));
    end
end

%---------------------------------------------------------------------------------------
% STEP 4: Execute
%---------------------------------------------------------------------------------------
spm_jobman('run', jobfile, inputs{:});

%---------------------------------------------------------------------------------------
% STEP 5: Setup BMC
%---------------------------------------------------------------------------------------

% Path to job file
jobfile = {fullfile(scrpth,'fmri','batch_peb_bmc_job.m')}; 

% Prepare inputs, 2 total, as per order listed in the job file
inputs  = {};

% Input #1: Path to estimated PEB file
inputs{1} = cellstr(fullfile(fmri_fits,sprintf('PEB_%s.mat','Full')));

% Input #2: Path to template GCM file with model space: 'Full' and 'Self'
inputs{2} = cellstr(fullfile(fmri_fits,'templates', 'GCMs', GCMname, sprintf('GCM_%s.mat',GCMname)));

% Execute
spm_jobman('run', jobfile, inputs{:});

%---------------------------------------------------------------------------------------
%
%     888     888         888 d8b      888          888            
%     888     888         888 Y8P      888          888            
%     888     888         888          888          888            
%     Y88b   d88P 8888b.  888 888  .d88888  8888b.  888888 .d88b.  
%      Y88b d88P     "88b 888 888 d88" 888     "88b 888   d8P  Y8b 
%       Y88o88P  .d888888 888 888 888  888 .d888888 888   88888888 
%        Y888P   888  888 888 888 Y88b 888 888  888 Y88b. Y8b.     
%         Y8P    "Y888888 888 888  "Y88888 "Y888888  "Y888 "Y8888
%
%---------------------------------------------------------------------------------------
%% Check DCM reproduces activations
% 
% load(fullfile(outpth,'GCM_full_self_fit.mat'));
% m = 1; % Cannot estimate model 2 because batch does not update GCM
% cols = [1:2]; 
% betas = []; dcm_betas = []; pvar = [];
% for s = 1:nsub
%     outdir = fullfile(outpth,subdir{s})
%     for r = 1:numel(ROI_names)        
%         %% Get GLM betas 
%         load(fullfile(outdir,'SPM.mat'));      
%         VOI =  load(fullfile(outdir,sprintf('VOI_%s_1.mat',ROI_names{r})),'Y');
%         tmp = pinv(SPM.xX.xKXs.X)*VOI.Y;
%         betas(s,r,:) = tmp(cols);
%         
%         %% Get DCM betas
%         load(GCM{s,m})
%         Y2 = DCM.y(:,r); % fitted data
%         Y2 = Y2/DCM.Y.scale;
%         tmp = pinv(SPM.xX.xKXs.X)*Y2;
%         dcm_betas(s,r,:) = tmp(cols);
%         
%         PSS   = sum(DCM.y(:,r).^2);
%         RSS   = sum(DCM.R(:,r).^2);
%         pvar(s,r) = PSS/(PSS + RSS);        
%     end
%     allF(s) = DCM.F;
% end
% figure,boxplot(pvar)
% %figure,boxplot(allF)
% mean(pvar)
% %mean(allF)
% 
% cw = [0 1]; % Faces only
% T=[]; p=[]; figure; clf
% for r = 1:numel(ROI_names)   
%     for c=1:size(cw,1)
%          dc(:,c) = squeeze(betas(:,r,:))*cw(c,:)';
%          T(c) = mean(dc(:,c))/(std(dc(:,c))/sqrt(size(dc,1)));
%     end  
%     p = t2p(T,size(dc,1)-1,2);
%     
%     fprintf('%s\n',ROI_names{r})
%     disp([cw T' p'])
%     
%     for c=1:size(cw,1)
%          dc(:,c) = squeeze(dcm_betas(:,r,:))*cw(c,:)';
%          T(c) = mean(dc(:,c))/(std(dc(:,c))/sqrt(size(dc,1)));
%     end  
%     p = t2p(T,size(dc,1)-1,2);
%     
%     fprintf('DCM: %s\n',ROI_names{r})
%     disp([cw T' p'])
%     
%     subplot(2,2,r),
%     mbeta = mean(squeeze(betas(:,r,:)))';
%     mbeta = mbeta/mean(mbeta)
%     mdcm_beta = mean(squeeze(dcm_betas(:,r,:)))';  
%     mdcm_beta = mdcm_beta/mean(mdcm_beta)
%     bar([mbeta mdcm_beta])
%     legend('real data','DCM fitted');
%     Title = ROI_names;
%     
%     title(Title{r},'Interpreter','none')
% end
% 
% b=2;
% for r = 1:numel(ROI_names)   
%     corr(betas(:,r,b),dcm_betas(:,r,b))
% end
