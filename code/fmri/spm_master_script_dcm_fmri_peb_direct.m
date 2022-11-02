%---------------------------------------------------------------------------------------
% Group Dynamic Causal Modelling of the Face Perception Network using fMRI
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

% Note that this script directly calls relevant SPM functions for executing tasks, and
% does not rely on the batch system exposed by SPM/MATLAB.

% Sections in this script are organized in the same order as covered in the tutorial.

% Authored in September 2022 by:
%   Rik Henson - rik.henson@mrc-cbu.cam.ac.uk
% With help from:
%   Pranay Yadav - pranay.yadav@mrc-cbu.cam.ac.uk

%---------------------------------------------------------------------------------------
% Data Sources
%---------------------------------------------------------------------------------------

% A. Processed DCM-ready data can be obtained from:
%   Yadav, Pranay; Henson, Rik (2022): Face processing fMRI data for Dynamic Causal
%   Modelling. figshare. Dataset. https://doi.org/10.6084/m9.figshare.21333996.v2 
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
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')

%---------------------------------------------------------------------------------------
% STEP 2: Configure & launch SPM 
%---------------------------------------------------------------------------------------

% Initialize SPM
spm('asciiwelcome');
spm_jobman('initcfg'); % Allows batch operations inside a script
spm_get_defaults('cmdline',true);
spm('defaults','fmri');

spm fmri

%---------------------------------------------------------------------------------------
% STEP 3: Variables for folders
%---------------------------------------------------------------------------------------

% Specify root working directory 
base_dir = '/imaging/henson/Wakeman/cognestic22_multimodal_dcm'; % Change this to yours

% Sub-directory containing scripts
srcpth = fullfile(base_dir,'code');
addpath(genpath(srcpth)) % Add scripts & functions to workspace

% Assign operational directories to variables
derpth = fullfile(base_dir, 'data','derivatives','SPM12');
% above will exist if you have run preprocessing scripts for Henson et al, 2019, 
% or else download and extract data from Figshare link

% All fits go in this directory
fits_dir = fullfile(base_dir, 'fits', 'direct_script', 'fmri');

%---------------------------------------------------------------------------------------
% STEP 4: Variables for data 
%---------------------------------------------------------------------------------------

% If you started with raw data uncomment the following lines...
% BIDS   = spm_BIDS(fullfile(base_dir, 'data'));
% subs   = spm_BIDS(BIDS,'subjects', 'task','facerecognition');
% runs = spm_BIDS(BIDS,'runs', 'modality','func', 'type','bold', 'task','facerecognition'); 

% ...else just re-specify
subs = compose('%02g', [1:9, 11:16]); % subject 10 had fewer scans in last run
nsub   = numel(subs);
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);
runs = compose('%02g', 1:9);
nrun = numel(runs);

% ...alternatively, infer from folder structure

% Specify subjects, runs for VOI extraction
%subdir = dir(fullfile(derpth, 'sub*'));
%subdir = {subdir.name};
%subs = strrep(subdir, 'sub-', '');
%nsub   = numel(subs);

nscan = repmat(208,1,nrun); % sub-15, run-01 has 209 scans, so ignore last
TR = 2;

%---------------------------------------------------------------------------------------
% STEP 5: Launch parallel pool
%---------------------------------------------------------------------------------------

numworkers = nsub; % Number of workers for distributed computing (depends on system)
if numworkers > 0
    delete(gcp('nocreate')) % Shut down any existing pool
    % Initialize and launch a parallel pool (only if running at the CBU)
    P=cbupool(numworkers, '--mem-per-cpu=4G --time=12:00:00 --nodelist=node-j10');
    parpool(P, P.NumWorkers);
    % Run the following line to initialize pool if script is not being run at the CBU
    % parpool(numworkers); 
end

% Proceed to the next step if you need to extract VOI files, else jump to 'DCM' section
% if you already have VOI timecourses.

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

% Could just load separate-run SPM.mat files, concatenate and re-save, but
% recreate below in case not been estimated before

% Care: if run in parfor, will hang if SPM.mat files already exist, so make
% sure delete them before re-running

% Parallellize across subjects
parfor (s = 1:nsub, numworkers)
    if numworkers > 0
        spm_jobman('initcfg');
        spm('defaults', 'fmri');
        spm_get_defaults('cmdline',true);
    end
    
    % Prepare SPM structure for this subject with concatenated info
    SPM = [];
    SPM.xY.RT  = 2;
    SPM.nscan = sum(nscan);
    SPM.xBF.T = 16;
    SPM.xBF.T0 = 8;
    SPM.xBF.name = 'hrf';
    SPM.xBF.UNITS = 'secs';
    SPM.xBF.Volterra = 1;
    
    % Prepare conditions struct
    conds = [];
    conds.names{1} = 'All';
    conds.names{2} = 'Faces';    % Famous and Nonfamous
    for t = 1:2; conds.onsets{t} = []; end
    
    % Initialize variables for tracking across runs
    time_so_far = 0; volfiles = {}; movepar = [];
    
    % Specify output directory for this subject
    outdir = fullfile(derpth,subdir{s}, 'func')   
    cd(outdir) 
      
    % Loop over runs
    for r = 1:length(runs)       
        
        % Get files with volumes for this run
        volfiles{r} = spm_select('ExtFPList',fullfile(derpth,subdir{s},'func'),sprintf('^swsub-.*run-%s_bold\\.nii$',runs{r}),[1:nscan(r)]);
    
        % Get files with trial info for this run
        trlfile = fullfile(derpth,subdir{s},'func',sprintf('sub-%s_run-%s_spmdef.mat',subs{s},runs{r}));
        d = load(trlfile);
        
        % Read and assign onsets to conditions struct
        conds.onsets{1} = [conds.onsets{1}; sort([d.onsets{1}; d.onsets{2}; d.onsets{3}]) + time_so_far];
        conds.onsets{2} = [conds.onsets{2}; sort([d.onsets{1}; d.onsets{2}]) + time_so_far];
        
        % Keep track of time elapsed across runs
        time_so_far = time_so_far + nscan(r)*TR;
        
        % Get files with movement info for this run
        d = load(spm_select('FPList',fullfile(derpth,subdir{s},'func'),sprintf('^rp.*run-%s.*\\.txt$',runs{r})));
        d = d(1:nscan(r),:);
        
        % Append movement parameters 
        movepar = [movepar; d];
        
    end
    
    % Specify concatenated model
    %----------------------------------------------------------------------    
    
    % Iterate over conditions and update SPM struct
    for c = 1:length(conds.names)
        SPM.Sess(1).U(c).name{1} = conds.names{c};
        SPM.Sess(1).U(c).ons = conds.onsets{c};
        SPM.Sess(1).U(c).dur = 0;
        SPM.Sess(1).U(c).P.name = 'none';
    end
    
    % Update SPM struct with movement parameters
    SPM.Sess(1).C.C = movepar;
    SPM.Sess(1).C.name = {'x';'y';'z';'pitch';'roll';'yaw'};
    
    % Update SPM struct with volume files
    SPM.xY.P = strvcat(volfiles{:});
    SPM.xX.K.HParam = 128;  
    SPM.xVi.form = 'AR(0.2)';
    SPM.xGX.iGXcalc = 'none';
    
    % Set up new GLM with concatenated data
    SPM = spm_fmri_spm_ui(SPM);
    
    % Call spm_fmri_concatenate to update SPM files for concatenated runs
    %----------------------------------------------------------------------     
    SPM = spm_fmri_concatenate('SPM.mat', nscan); 
    delete('SPM_backup.mat');
    
    % Estimate new model
    %----------------------------------------------------------------------    
    SPM = spm_spm(SPM);
    
    % Evaluate Effects of Interest contrast
    try SPM = rmfield(SPM,'xCon'); end;
    c = [eye(2) zeros(2,size(SPM.xX.X,2)-2)]';
    SPM.xCon(1)   = spm_FcUtil('Set','Effects of Interest','F','c',c,SPM.xX.xKXs);
    spm_contrasts(SPM);
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

% Specify names and locations of ROIs
ROI_names = {'lOFA','rOFA','lFFA','rFFA'};
ROI_coord = {[-38, -86, -14],[+36, -86, -10],[-42, -56, -20],[+42, -52, -14]};
rad = 10; % radius

% Preload SPMs just so parfor can work without loading
allSPMs = {};
for s = 1:nsub
    tmp = load(fullfile(derpth,subdir{s}, 'func', 'SPM.mat'));
    allSPMs{s} = tmp.SPM;
end

% Loop in parallell over subjects
parfor (s = 1:nsub, numworkers) % can't parallelise with needing to load SPM
    if numworkers > 0
        spm_jobman('initcfg');
        spm('defaults', 'fmri');
        spm_get_defaults('cmdline',true);
    end
    
    % Output directory for subject
    outdir = fullfile(derpth,subdir{s}, 'func')
    SPM = allSPMs{s};
%     spmfile = fullfile(outdir,'SPM.mat');
%     load(spmfile)
    
    % Define coordinate space
    [x,y,z] = ndgrid(1:SPM.xVol.DIM(1),1:SPM.xVol.DIM(2),1:SPM.xVol.DIM(3));
    XYZ     = [x(:),y(:),z(:)]'; 
    XYZmm   = SPM.xVol.M(1:3,:) * [XYZ;ones(1,size(XYZ,2))];
    Q       = ones(1,size(XYZmm,2));   
    
    % Prepare SPM structure for VOIs
    xSPM = [];
    xSPM.Ic        = 1;
    xSPM.n         = 0;
    xSPM.u         = 0.001;
    xSPM.k         = 0;
    xSPM.title     = '';
    xSPM.Im        = [];    
    xSPM.XYZmm = XYZmm;
    xSPM.XYZ   = XYZ;
    xSPM.M     = SPM.xVol.M; % irrelevant here
    
    % Iterate over ROIs
    for r = 1:length(ROI_names)
        
        % Initialize VOI structure
        xY = [];
        xY.Ic = 1; xY.Sess = 1;
        xY.xyz = []'; xY.def = 'mask';
        
        xSPM.thresDesc = 'none'; % (reset because overwritten below)

        % Within a sphere from coordinates..
        voi1           = zeros(SPM.xVol.DIM');
        voi1(sum((XYZmm - ROI_coord{r}'*Q).^2) <= rad^2) = 1;
        
        % ...and shows significant task effect
        voi2           = zeros(SPM.xVol.DIM');
        xSPM.swd       = outdir;
        [SPM1,xSPM]    = spm_getSPM(xSPM);
        voi2(sub2ind(SPM1.xVol.DIM',xSPM.XYZ(1,:),xSPM.XYZ(2,:),xSPM.XYZ(3,:))) = 1;
        voi2(spm_sample_vol(voi2, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0) = 1;
        
        % then combine and write out mask image...
        voi = voi1.*voi2;
        Vm = struct('fname', fullfile(outdir, ['VOI_' ROI_names{r} '_mask' spm_file_ext]), ...
            'dim',     SPM.xVol.DIM', ...
            'dt',      [spm_type('uint8') spm_platform('bigend')], ...
            'mat',     SPM.xVol.M, ...
            'pinfo',   [1 0 0]', ...
            'descrip', 'VOI');
        Vm = spm_write_vol(Vm,voi);

        % Extract first eigenvariate
        xY.name    = ROI_names{r};    
        xY.spec    = Vm;   
        SPM.swd    = outdir;
        [Y,xY]     = spm_regions(xSPM,SPM,[],xY);
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

ref_sub = 15; % eg first subject
outdir = fullfile(derpth,subdir{ref_sub},'func');
load(fullfile(outdir,'SPM.mat'));

% Initialize empty DCM structure
DCM = [];

% Specify names of ROIs in order
ROI_names = {'lOFA','rOFA','lFFA','rFFA'};

% Populate VOIs for each ROI
for r = 1:numel(ROI_names)
    load(fullfile(outdir,sprintf('VOI_%s_1.mat',ROI_names{r})),'xY');
    DCM.xY(r) = xY;
end      
      
% Specify descriptions of data
DCM.n = length(DCM.xY);      % number of regions
DCM.v = length(DCM.xY(1).u); % number of time points

DCM.Y.dt  = SPM.xY.RT;
DCM.Y.X0  = DCM.xY(1).X0;

% Add time courses for each VOI
for i = 1:DCM.n
    DCM.Y.y(:,i)  = DCM.xY(i).u;
    DCM.Y.name{i} = DCM.xY(i).name;
end

% Model autocorrelation
DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v);

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
        
%---------------------------------------------------------------------------------------
% STEP 3: Connectivity
%---------------------------------------------------------------------------------------

% Full model

DCM.a = [
%    lOFA rOFA lFFA rFFA
    [  1    1    1    0  ];   % lOFA
    [  1    1    0    1  ];   % rOFA
    [  1    0    1    1  ];   % lFFA
    [  0    1    1    1  ];   % rFFA    
];
DCM.c = [
%     All  Faces  
    [  1    0  ];   % lOFA
    [  1    0  ];   % rOFA
    [  0    0  ];   % lFFA
    [  0    0  ];   % rFFA
];
DCM.b(:,:,1) = zeros(DCM.n,DCM.n);
DCM.b(:,:,2) = [
%    lOFA rOFA lFFA rFFA
    [  1    0    1    0  ];   % lOFA
    [  0    1    0    1  ];   % rOFA
    [  1    0    1    0  ];   % lFFA
    [  0    1    0    1  ];   % rFFA    
];
DCM.d        = zeros(DCM.n,DCM.n,0); % needed else crashes in estimation

%---------------------------------------------------------------------------------------
% STEP 4: Save
%---------------------------------------------------------------------------------------

outfile_full = fullfile(fits_dir, 'templates', 'DCMs', 'DCM_Full.mat');
save(outfile_full,'DCM');

%---------------------------------------------------------------------------------------
% STEP 5: Estimate (Optional)
%---------------------------------------------------------------------------------------

% Estimate reference subject's DCM and save
DCM = spm_dcm_estimate(DCM);
save(fullfile(fits_dir, sprintf('DCM_Full_sub-%02g.mat', ref_sub)),'DCM');

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

models = {'Full'}; % Could specify multiple models here

% Load template models
model = load(outfile_full); % This will load the 'Full' model
DCM_Full = model.DCM;

GCM = {};
for s = 1:nsub
    
    % Path to this subject's folder under derivatives
    outdir = fullfile(derpth,subdir{s}, 'func');
    
    % Iterate over models (we only have one here)
    for m = 1:length(models)
        
        % Specify output file for DCM specification
        DCMfile = fullfile(outdir,['DCM_' models{m} '.mat']);
        
        % Make copy of full model
        DCM = DCM_Full;
        
        % Update VOI data
        for r = 1:numel(ROI_names)
            load(fullfile(outdir,sprintf('VOI_%s_1.mat',ROI_names{r})),'xY');
            DCM.xY(r) = xY;
        end 
        
        % Set up DCM.Y
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end

        % Update onsets    
        load(fullfile(outdir,'SPM.mat'));
        for u = 1:2
            DCM.U.u(:,u)  = SPM.Sess(1).U(u).u((32+1):end); % DCM allows for 2 TRs before first stimulus
            DCM.U.name{u} = SPM.Sess(1).U(u).name{1};
        end

        % If some alternatives to update B (could update A,C too!)
        if m>1
            DCM.b(:,:,2) = altpar(m).B;
        end
        
        % Save this DCM specified for this subject
        save(DCMfile,'DCM');
        
        % Append to GCM array
        GCM{s,m} = DCMfile;
        
    end
end

% Save GCM
save(fullfile(fits_dir, 'templates', 'GCMs', 'Full', 'GCM_Full.mat'),'GCM');    

%---------------------------------------------------------------------------------------
% OUTPUTS
%---------------------------------------------------------------------------------------
% Running this section will produce the following output in the folder 'fits_dir'
% 1. GCM specification file called 'GCM_Full.mat' under fits_dir/templates/GCMs/Full
%       This is the GCM array with full DCM models which has not yet been fitted

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
% Estimate specified DCMs for all subjects using direct function calls

GCM_fit = spm_dcm_fit(GCM(:,1),(numworkers > 0));

% This GCM_fit cell array can be saved directly but will take up a lot of space
% save(['GCM_' models{1}],'GCM'); % This saves big file with all estimated
% values, but to mimic batch interface, update the DCM files in the subject
% directories instead (note files not appended with "_m0001" etc like in
% batched version)

% Instead iterate over models in GCM and save them individually
for s=1:size(GCM,1) 
    DCM = GCM_fit{s,1}; % GCM_fit contains estimated DCMs in an array
    save(GCM{s,1},'DCM') % GCM contains paths to DCM files in an array
end

% Clear the GCM_fit variable as it takes up a lot of space in memory
clear GCM_fit

% Save a copy of estimated GCM for consistency with steps in the batch pipeline
save(fullfile(fits_dir, 'GCM_Full.mat'), 'GCM')

% Check fits
spm_dcm_fmri_check(GCM);

%---------------------------------------------------------------------------------------
% OUTPUTS
%---------------------------------------------------------------------------------------
% Running this section will produce the following outputs in the folder 'fits_dir'
% 1. Estimated GCM file called 'GCM_Full.mat' under fits_dir
%       This is the GCM array consisting of paths to files with fitted DCM models, one
%       row per subject
% 2. Additionally, DCM fits for each subject are also stored in data/derivatives/SPM12
% in their respective directories with the filename 'DCM_Full_m0001.mat'.
%       These are the same as the estimated GCM file above but are one per subject (15
%       total) instead of an array like GCM. Useful for quick inspection in the DCM GUI. 

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
M = []; M.Q = 'all'; % Random effects over all parameters in field 'B'
PEB = spm_dcm_peb(GCM, M, {'B'});
save(fullfile(fits_dir, 'PEB_Full'), 'PEB')

% Review fitted PEB
spm_dcm_peb_review(PEB,GCM);

%---------------------------------------------------------------------------------------
% OUTPUTS
%---------------------------------------------------------------------------------------
% Running this section will produce the following output in the folder 'fits_dir'
% 1. Estimated PEB file called 'PEB_Full.mat' under fits_dir
%       This is the group PEB estimated from the fitted GCM in the previous section

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
% STEP 1: Bayesian Model Selection
%---------------------------------------------------------------------------------------
[BMA, BMR] = spm_dcm_peb_bmc(PEB);

% Write to disk
save(fullfile(fits_dir, 'BMA_search_PEB_Full'), 'BMA')

%---------------------------------------------------------------------------------------
% STEP 2: Review estimated BMA
%---------------------------------------------------------------------------------------
spm_dcm_peb_review(BMA, GCM)

%---------------------------------------------------------------------------------------
% OUTPUTS
%---------------------------------------------------------------------------------------
% Running this section will produce the following output in the folder 'fits_dir'
% 1. BMA file called 'BMA_search_PEB_Full.mat' under fits_dir
%       This is the BMA obtained after averaging reduced models that contribute
%       significantly to model evidence.

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
% STEP 1: Bayesian Model Selection
%---------------------------------------------------------------------------------------

% Load full model
model = load(outfile_full);
DCM_Full = model.DCM;

% Make a copy of DCM for editing
DCM = DCM_Full;

% Switch off between-region connections
DCM.b(:,:,2) = [
     1     0     0     0
     0     1     0     0
     0     0     1     0
     0     0     0     1
];

% Save reduced model as a template 
outfile_self = fullfile(fits_dir, 'templates', 'DCMs','DCM_Self.mat');
save(outfile_self,'DCM');

DCM_Self = DCM;

% Construct model space with 'full' and 'self' models
GCM = {DCM_Full, DCM_Self};

% Save model space
save(fullfile(fits_dir,'templates', 'GCMs', 'Full_vs_Self', 'GCM_Full_vs_Self.mat'), 'GCM')

% Perform reduction of the nested 'self' model and estimate BMA
[BMA, BMR] = spm_dcm_peb_bmc(PEB, GCM);

% Write to disk
save(fullfile(fits_dir, 'BMA_PEB_Full_vs_Self'), 'BMA')

%---------------------------------------------------------------------------------------
% STEP 2: Review estimated BMA
%---------------------------------------------------------------------------------------
spm_dcm_peb_review(BMA, GCM)

%---------------------------------------------------------------------------------------
% OUTPUTS
%---------------------------------------------------------------------------------------
% Running this section will produce the following outputs in the folder 'fits_dir'
% 1. Model space as a GCM array called 'GCM_Full_vs_Self.mat' under the folder
% fits_dir/templates/GCMs/Full_vs_Self
%       This is a single row array with two columns, the first of which corresponds to
%       the full model while the second one corresponds to the nested self-only model
% 2. BMA file called 'BMA_PEB_Full_sv_Self.mat' under fits_dir
%       This is the BMA obtained after taking a weighted average of the full and reduced
%       self models. 

%---------------------------------------------------------------------------------------
%
%     8888888888                     d8b 888 d8b                   
%     888                            Y8P 888 Y8P                   
%     888                                888                       
%     8888888  8888b.  88888b.d88b.  888 888 888  .d88b.  .d8888b  
%     888         "88b 888 "888 "88b 888 888 888 d8P  Y8b 88K      
%     888     .d888888 888  888  888 888 888 888 88888888 "Y8888b. 
%     888     888  888 888  888  888 888 888 888 Y8b.          X88 
%     888     "Y888888 888  888  888 888 888 888  "Y8888   88888P'
%
%---------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------
% STEP 1: Define model space
%---------------------------------------------------------------------------------------

% As usual, our model space is a GCM
GCM = {};

% Load the 'Full' template DCM we specified earlier
% load(fullfile(fits_dir, 'templates', 'DCMs', 'DCM_Full.mat')) % Direct path
model = load(outfile_full); % We already assigned the direct path to this variable
DCM_Full = model.DCM;

% Remove priors if present, they interfere with the internal model comparison code
if isfield(DCM_Full, 'M')
    DCM_Full = rmfield(DCM_Full, 'M');
end

% 1. Model F+B+S (= Full model)
GCM{1, 1} = DCM_Full;

% 2. Model F+S (= No-backward, with Self)
% Get full DCM specification
DCM = DCM_Full;

% Switch off Backward connections in B-matrix
DCM.b(:,:,2) =  [
%    lOFA rOFA lFFA rFFA
    [  1    1    0    0  ];   % lOFA
    [  1    1    0    0  ];   % rOFA
    [  1    0    1    1  ];   % lFFA
    [  0    1    1    1  ];   % rFFA
];

% Append to GCM
GCM{1, 2} = DCM;

% 3. Model B+S (= No-forward, with Self)
% Get full DCM specification
DCM = DCM_Full;

% Switch off Forward connections in B-matrix
DCM.b(:,:,2) =  [
%    lOFA rOFA lFFA rFFA
    [  1    1    1    0  ];   % lOFA
    [  1    1    0    1  ];   % rOFA
    [  0    0    1    1  ];   % lFFA
    [  0    0    1    1  ];   % rFFA
];

% Append to GCM
GCM{1, 3} = DCM;

% 4. Model S (= Neither forward nor backward, only Self)
% Get full DCM specification
DCM = DCM_Full;

% Switch off Forward and Backward connections in B-matrix
DCM.b(:,:,2) =  [
%    lOFA rOFA lFFA rFFA
    [  1    1    0    0  ];   % lOFA
    [  1    1    0    0  ];   % rOFA
    [  0    0    1    1  ];   % lFFA
    [  0    0    1    1  ];   % rFFA
];

% Append to GCM
GCM{1, 4} = DCM;

% 5. Model F+B (= Both forward and backward, without Self)
% Get full DCM specification
DCM = DCM_Full;

% Switch off Self connections in B-matrix
DCM.b(:,:,2) =  [
%    lOFA rOFA lFFA rFFA
    [  0    1    1    0  ];   % lOFA
    [  1    0    0    1  ];   % rOFA
    [  1    0    0    1  ];   % lFFA
    [  0    1    1    0  ];   % rFFA
];

% Append to GCM
GCM{1, 5} = DCM;

% 6. Model F (= No-backward, without Self)
% Get full DCM specification
DCM = DCM_Full;

% Switch off Backward connections in B-matrix
DCM.b(:,:,2) =  [
%    lOFA rOFA lFFA rFFA
    [  0    1    0    0  ];   % lOFA
    [  1    0    0    0  ];   % rOFA
    [  1    0    0    1  ];   % lFFA
    [  0    1    1    0  ];   % rFFA
];

% Append to GCM
GCM{1, 6} = DCM;

% 7. Model B (= No-forward, without Self)
% Get full DCM specification
DCM = DCM_Full;

% Switch off Forward and Self connections in B-matrix
DCM.b(:,:,2) =  [
%    lOFA rOFA lFFA rFFA
    [  0    1    1    0  ];   % lOFA
    [  1    0    0    1  ];   % rOFA
    [  0    0    0    1  ];   % lFFA
    [  0    0    1    0  ];   % rFFA
];

% Append to GCM
GCM{1, 7} = DCM;

% 8. Model with no F/B/S (= Null model)
% Get full DCM specification
DCM = DCM_Full;

% Switch off Forward, Backward & Self connections in B-matrix
DCM.b(:,:,2) =  [
%    lOFA rOFA lFFA rFFA
    [  0    1    0    0  ];   % lOFA
    [  1    0    0    0  ];   % rOFA
    [  0    0    0    1  ];   % lFFA
    [  0    0    1    0  ];   % rFFA
];

% Append to GCM
GCM{1, 8} = DCM;

% Save model space
gcm_families_file = fullfile(fits_dir, 'templates', 'GCMs', 'Families', 'GCM_ModelSpace8.mat');
save(gcm_families_file, 'GCM')

% Visualize model space
figure;
for k=1:8
    subplot(2,4,k);
    imagesc(GCM{1, k}.b(:,:,2) + 0.75*fliplr(eye(length(ROI_names))));
    colormap(gray)
    caxis([0, 1])
    title(sprintf('Model %02d', k))
    axis square
end

%---------------------------------------------------------------------------------------
% STEP 2: Load estimated PEB and perform BMR of model space
%---------------------------------------------------------------------------------------

% Load estimated PEB from file
load(fullfile(fits_dir, 'PEB_Full.mat'))

% Bayesian Model Reduction (BMR) and comparison of models
[BMA, BMR] = spm_dcm_peb_bmc(PEB, GCM);

% Save BMA and BMR
outfile = fullfile(fits_dir, 'BMA_BMR_Families.mat');
save(outfile, 'BMA', 'BMR')

%---------------------------------------------------------------------------------------
% STEP 3: Group models into families and compare
%---------------------------------------------------------------------------------------
% Now partition the model space into families and perform comparisons at the level of
% families to test hypotheses about modulation of connection groups due to faces

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESIS 1
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Are any between-region connections modulated regardless of self-connections?
% Family 1: Models 1, 2, 3 & 5, 6, 7 have at least one forward or backward connection
% Family 2: Models 4 and 8 have no F/B connections
families = [1, 1, 1, 2, 1, 1, 1, 2];
[BMAf, fam] = spm_dcm_peb_bmc_fam(BMA, BMR, families, 'NONE');
% Family 1 has overwhelming evidence (~1) -> between-region connections are modulated

% Save this family-wise comparison
outfile = fullfile(fits_dir, 'BMC_Families_BetweenRegion.mat');
save(outfile, 'BMAf', 'fam')

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESIS 2a
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Are any forward connections modulated regardless of backward or self-connections?
% Family 1: Models 1, 2, 5, 6 have at least one forward connection
% Family 2: Models 3, 4, 7, 8 have no forward connection
families = [1, 1, 2, 2, 1, 1, 2, 2];
[BMAf, fam] = spm_dcm_peb_bmc_fam(BMA, BMR, families, 'NONE');
% Family 1 has significantly greater evidence (>0.95)

% Save this family-wise comparison
outfile = fullfile(fits_dir, 'BMC_Families_BetweenRegion_Forward.mat');
save(outfile, 'BMAf', 'fam')

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESIS 2b
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Are any backward connections modulated regardless of forward or self-connections?
% Family 1: Models 1, 3, 5, 7 have at least one backward connection
% Family 2: Models 2, 4, 6, 8 have no backward connection
families = [1, 2, 1, 2, 1, 2, 1, 2];
[BMAf, fam] = spm_dcm_peb_bmc_fam(BMA, BMR, families, 'NONE');
% Family 1 has overwhelming evidence (~1)

% Save this family-wise comparison
outfile = fullfile(fits_dir, 'BMC_Families_BetweenRegion_Backward.mat');
save(outfile, 'BMAf', 'fam')

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HYPOTHESIS 3
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Are any self connections modulated regardless of forward or backward connections?
% Family 1: Models 1, 2, 3, 4 have at least one self connection
% Family 2: Models 5, 6, 7, 8 have no self connection
families = [1, 1, 1, 1, 2, 2, 2, 2];
[BMAf, fam] = spm_dcm_peb_bmc_fam(BMA, BMR, families, 'NONE');
% Family 1 has overwhelming evidence (~1)

% Save this family-wise comparison
outfile = fullfile(fits_dir, 'BMC_Families_Self.mat');
save(outfile, 'BMAf', 'fam')

%---------------------------------------------------------------------------------------
% OUTPUTS
%---------------------------------------------------------------------------------------
% Running this batch job will produce the following output in the folder 'fits_dir'
% 1. GCM model space file 'GCM_ModelSpace8.mat' in fits_dir/templates/GCMs/Families
%       This file consists of the 8 models we specified as columns of the GCM cell array
% 2. BMA and BMR in the file 'BMA_BMR_Families.mat' in fits_dir
%       This file consists of both BMA and BMR variables which represent the average
%       over all 8 models and the 8 reduced models respectively.
% 3. Four files in fits_dir, one for each hypothesis and family-wise comparison:
%       i. 'BMC_Families_BetweenRegion.mat': Modulation of any between-region connections
%       ii. 'BMC_Families_BetweenRegion_Forward.mat': Modulation of any forward connections
%       iii. 'BMC_Families_BetweenRegion_Backward.mat': Modulation of any backward connections
%       iv. 'BMC_Families_Self.mat': Modulation of any self-connections

%---------------------------------------------------------------------------------------
%
%      .d8888b.                                     d8b          888                     
%     d88P  Y88b                                    Y8P          888                     
%     888    888                                                 888                     
%     888         .d88b.  888  888  8888b.  888d888 888  8888b.  888888 .d88b.  .d8888b  
%     888        d88""88b 888  888     "88b 888P"   888     "88b 888   d8P  Y8b 88K      
%     888    888 888  888 Y88  88P .d888888 888     888 .d888888 888   88888888 "Y8888b. 
%     Y88b  d88P Y88..88P  Y8bd8P  888  888 888     888 888  888 Y88b. Y8b.          X88 
%      "Y8888P"   "Y88P"    Y88P   "Y888888 888     888 "Y888888  "Y888 "Y8888   88888P'
%
%---------------------------------------------------------------------------------------
% We demonstrate the inclusion of covariates for 2nd-level (group) inference with PEB.
% The dataset includes ages of participants, and while we do not anticipate any effect
% of age on modulation of connections due to faces, we illustrate specification of age
% as a covariate in the PEB design matrix for inference.

%---------------------------------------------------------------------------------------
% STEP 1: Setup Design Matrix with Covariates
%---------------------------------------------------------------------------------------

% Define covariates, and assign appropriate labels
PEB_name = 'Age';
covariate_name = 'Age';
covariate_values = [31, 25, 30, 26, 23, 26, 31, 26, 29, 24, 24, 25, 24, 30, 25]';
%covariate_name = 'Sex';
%covariate_values = [0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0]'; % 1=Female

% Mean-center the covariate (Optional)
covariate_values = covariate_values - mean(covariate_values);

% Design Matrix
M.X = [ones([length(covariate_values), 1]), covariate_values]; % First covariate is group mean
M.Xnames = {'Commonalities', 'Age'};
M.Q = 'all'; % Random effects over all parameters

% Load fitted GCM
load(fullfile(fits_dir, 'GCM_Full'))

%---------------------------------------------------------------------------------------
% STEP 2: Fit 2nd-level PEB
%---------------------------------------------------------------------------------------
PEB = spm_dcm_peb(GCM, M, {'B'});

% Write to disk
save(fullfile(fits_dir, 'PEB_Age'), 'PEB')

%---------------------------------------------------------------------------------------
% STEP 3: Review
%---------------------------------------------------------------------------------------
spm_dcm_peb_review(PEB, GCM)

% Greedy search for nested models on this PEB can be done to perform inference and
% identify which connections modulated by faces are affected by age.

%---------------------------------------------------------------------------------------
% OUTPUTS
%---------------------------------------------------------------------------------------
% Running this batch job will produce the following outputs in the folder 'fits_dir'
% 1. Estimated PEB file called 'PEB_Age.mat' under fits_dir
%       This is the group PEB estimated with age as a covariate.

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

% Check DCM reproduces activations
% 
% load(fullfile(fits_dir,'GCM_Full.mat'));
% m = 1; % Cannot estimate model 2 because batch does not update GCM
% cols = [1:2]; 
% betas = []; dcm_betas = []; pvar = [];
% for s = 1:nsub
%     outdir = fullfile(derpth,subdir{s}, 'fmri', 'stat_concat')
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
