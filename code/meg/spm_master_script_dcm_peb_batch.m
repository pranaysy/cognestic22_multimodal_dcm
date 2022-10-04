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

% Add SPM12 to MATLAB Path
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')

%---------------------------------------------------------------------------------------
% STEP 2 
%---------------------------------------------------------------------------------------

% Initialize SPM
spm('asciiwelcome');
spm_jobman('initcfg'); % Allows batch operations inside a script
%spm_get_defaults('cmdline',true);
spm('defaults','EEG');

%---------------------------------------------------------------------------------------
% STEP 3 
%---------------------------------------------------------------------------------------

% Specify root working directory 
base_dir = '/imaging/henson/Wakeman/cognestic22_multimodal_dcm/'; % Change this to yours
addpath(genpath(fullfile(base_dir, 'code'))) % Add scripts & functions to workspace

% All fits go in this directory
fits_dir = fullfile(base_dir, 'fits', 'batch_script', 'meg');

% All code is present in this directory
code_dir = fullfile(base_dir, 'code', 'meg');

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
% Configure DCM model with data and options

% Name of this DCM model (Used for each subject's fits, and downstream by GCM/PEB)
name = 'DCM_Full';
DCM.name = name;

%---------------------------------------------------------------------------------------
% STEP 1: Setup analysis options
%---------------------------------------------------------------------------------------

% Specify modality
DCM.xY.modality = 'MEG';

% Set up DCM analysis type
DCM.options.analysis = 'ERP';   % Analyze evoked responses
DCM.options.model    = 'ERP';   % Neuronal temporal model: Extended Jansen-Rit model
DCM.options.spatial  = 'ECD';   % Spatial observation model: ECD

% Set up preprocessing parameters and analysis options
DCM.options.Nmodes   = 8;       % Number of modes of Leadfield for data selection
DCM.options.CVA      = 0;       % Optimize modes of Leadfield
DCM.options.h        = 1;       % Number of DCT components for detrending (1 is mean)
DCM.options.han      = 1;       % Hanning Window Taper
DCM.options.onset    = 64;      % Selection of onset (prior mean) for input stimulus
DCM.options.dur      = 16;      % Duration of onset (prior sd) for input stimulus
DCM.options.D        = 1;       % Downsampling (decimation of time series by a factor)
DCM.options.multiC   = 0;       % Multiple input vectors for multiple stimuli
DCM.options.location = 0;       % Optimize dipole locations
DCM.options.symmetry = 1;       % Lock orientation of dipoles across hemispheres

%---------------------------------------------------------------------------------------
% STEP 2: Setup data & design
%---------------------------------------------------------------------------------------

% Specify data of interest
DCM.options.trials   = [1 2 3]; % Index of ERPs within ERP/ERF file
DCM.options.Tdcm     = [0 500]; % Peri-stimulus time to be modelled

% Specify between-condition trial effects
contrasts = [1 1 0]'; % Face Perception: Faces (Famous + Unfamiliar) - Scrambled
DCM.xU.X = contrasts; % Orientation is N_trials x N_contrasts
DCM.xU.name = {'Face Perception'};

%--------------------------------------------------------------------------
% STEP 3: Setup observation model
%--------------------------------------------------------------------------

% Location priors for dipoles
locs  = {
    [-38, -86, -14], 'lOFA';
    [+36, -86, -10], 'rOFA';
        
    [-42, -56, -20], 'lFFA';
    [+42, -52, -14], 'rFFA';   
};

DCM.Lpos  = cat(1, locs{:,1})';
DCM.Sname = locs(:,2)';
Nareas    = length(locs);

%--------------------------------------------------------------------------
% STEP 4: Setup neuronal model
%--------------------------------------------------------------------------

% A Matrix: Forward connections
DCM.A{1} = [
%    lOFA rOFA lFFA rFFA
    [  0    0    0    0  ];   % lOFA
    [  0    0    0    0  ];   % rOFA
    [  1    0    0    0  ];   % lFFA
    [  0    1    0    0  ];   % rFFA    
];

% A Matrix: Backward connections
DCM.A{2} = [
%    lOFA rOFA lFFA rFFA
    [  0    0    1    0  ];   % lOFA
    [  0    0    0    1  ];   % rOFA
    [  0    0    0    0  ];   % lFFA
    [  0    0    0    0  ];   % rFFA
];

% A Matrix: Lateral connections
DCM.A{3} = [
%    lOFA rOFA lFFA rFFA
    [  0    1    0    0  ];   % lOFA
    [  1    0    0    0  ];   % rOFA
    [  0    0    0    1  ];   % lFFA
    [  0    0    1    0  ];   % rFFA
];

% B Matrix: Modulation of connections
self_connections = eye(Nareas);
DCM.B{1} = double(DCM.A{1} | DCM.A{2} | DCM.A{3} | self_connections);

% C Matrix: Driving inputs
DCM.C = [1 1 0 0]';

% Save full model as a template
dcm_full_file = fullfile(fits_dir, 'templates', 'DCMs', strcat(DCM.name, '.mat'));
save(dcm_full_file, 'DCM')
DCM_Full = DCM;

%--------------------------------------------------------------------------
% STEP 5: Specify reduced models, if any
%--------------------------------------------------------------------------

% Reduced model with modulation of only self connections
DCM.name = 'DCM_Self';
DCM.B{1} = self_connections;

% Save reduced model as a template
dcm_self_file = fullfile(fits_dir, 'templates', 'DCMs', strcat(DCM.name, '.mat'));
save(dcm_self_file, 'DCM')

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
% Estimate specified DCMs for all subjects using the batch interface

%---------------------------------------------------------------------------------------
% STEP 1: Prepare inputs for batch job 
%---------------------------------------------------------------------------------------

% Cell array of paths to DCM files
dcm_files = {dcm_full_file}; % We fit only the full model here

% String identifier for estimated GCM & PEB files
name_tag = 'Full';

% Populate list of processed files as a column-order cell (N-files × 1)
% These files should contain forward models (with or without gain matrices)
files = dir(fullfile(base_dir, 'data', '**', 'maM*.mat'));
input_files = arrayfun(@(x) fullfile(x.folder, x.name), files, 'UniformOutput', false);

% Specify output directory for GCM specification
out_dir_spec = {fullfile(fits_dir, 'templates', 'GCMs', 'Full')};

% Specify output directory for estimated GCM 
out_dir_estim = {fits_dir};

%---------------------------------------------------------------------------------------
% STEP 2: Setup jobfile and inputs 
%---------------------------------------------------------------------------------------

% Specify jobfile
jobfile = {fullfile(code_dir, 'batch_dcm_peb_estimation_job.m')};

% Initialize inputs to jobfile as a columnar cell array
% There are 6 inputs needed by this jobfile
inputs = cell(6, 1);

% Populate inputs in order indicated in the jobfile (generated by the batch manager)
inputs{1, 1} = input_files; % DCM for M/EEG: M/EEG datasets - cfg_files (Should be cell array, row-order)
inputs{2, 1} = dcm_files; % DCM for M/EEG: DCM files - cfg_files (Should be cell array, not string/char)
inputs{3, 1} = out_dir_spec; % DCM for M/EEG: Directory - cfg_files (Should be cell array, not string/char)
inputs{4, 1} = out_dir_estim; % DCM estimation: Directory - cfg_files (Should be cell array, not string/char)
inputs{5, 1} = name_tag; % DCM estimation: Name - cfg_entry ('GCM_' suffix will be appended)
inputs{6, 1} = name_tag; % Specify / Estimate PEB: Name - cfg_entry ('PEB_' suffix will be appended)

%---------------------------------------------------------------------------------------
% STEP 3: Execute batch job using parallel compute backend
%---------------------------------------------------------------------------------------

% Initialize Parallel Compute Pool (Example Instructions for CBU Cluster)
delete(gcp('nocreate')) % Shut down any existing pool
n_workers = length(input_files);
P=cbupool(n_workers, '--mem-per-cpu=4G --time=12:00:00 --nodelist=node-j10');
parpool(P, P.NumWorkers);

% Call jobman with the jobfile and necessary inputs
spm_jobman('run', jobfile, inputs{:});

%---------------------------------------------------------------------------------------
% OUTPUTS
%---------------------------------------------------------------------------------------
% Running this batch job will produce the following outputs in the folder 'fits_dir'
% 1. GCM specification file called 'GCM_DCM_Full.mat' under fits_dir/templates/GCMs/Full
%       This is the GCM array with full DCM models which has not yet been fitted
% 2. Estimated GCM file called 'GCM_Full.mat' under fits_dir
%       This is the GCM array consisting of fitted DCM models
% 3. Estimated PEB file called 'PEB_Full.mat' under fits_dir
%       This is the group PEB estimate from the last step in the batch pipeline
% 4. Additionally, individual DCM fits are also stored in fits_dir/templates/GCMs/Full
%       These are the same as the estimated GCM file in (2) but are one per subject (16
%       total) instead of an array like GCM. Useful for quick inspection in the DCM GUI. 

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
% STEP 1: Prepare inputs for batch job 
%---------------------------------------------------------------------------------------

% Paths to estimated PEB and GCM structures
PEB_file = cellstr(fullfile(fits_dir, sprintf('PEB_%s.mat', name_tag)));
GCM_file = cellstr(fullfile(fits_dir, sprintf('GCM_%s.mat', name_tag)));

%---------------------------------------------------------------------------------------
% STEP 2: Setup jobfile and inputs 
%---------------------------------------------------------------------------------------

% Specify jobfile
jobfile = {fullfile(code_dir, 'batch_dcm_peb_search_job.m')}; 

% Initialize inputs to jobfile as a columnar cell array
% There are 2 inputs needed by this jobfile
inputs = cell(2, 1);

% Populate inputs in order indicated in the jobfile (generated by the batch manager)
inputs{1} = PEB_file;
inputs{2} = GCM_file;

%---------------------------------------------------------------------------------------
% STEP 3: Execute batch job
%---------------------------------------------------------------------------------------

% Call jobman with the jobfile and necessary inputs
spm_jobman('run', jobfile, inputs{:});

%---------------------------------------------------------------------------------------
% OUTPUTS
%---------------------------------------------------------------------------------------
% Running this batch job will produce the following output in the folder 'fits_dir'
% 1. BMA file called 'BMA_search_PEB_Full.mat' under fits_dir
%       This is the BMA obtained after averaging reduced models that contribute
%       significantly to model evidence. Can be reviewed by calling spm_dcm_peb_review.

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
% STEP 1: Prepare inputs for batch job 
%---------------------------------------------------------------------------------------

% Paths to DCM models and estimated PEB
PEB_file = cellstr(fullfile(fits_dir, sprintf('PEB_%s.mat', name_tag)));
DCM_files = {dcm_full_file; dcm_self_file};

% Path to GCM specification
out_dir_spec = {fullfile(fits_dir, 'templates', 'GCMs', 'Full_vs_Self')};

%---------------------------------------------------------------------------------------
% STEP 2: Setup jobfile and inputs 
%---------------------------------------------------------------------------------------

% Specify jobfile
jobfile = {fullfile(code_dir, 'batch_dcm_peb_bmc_modelspace_job.m')}; 

% Initialize inputs to jobfile as a columnar cell array
% There are 4 inputs needed by this jobfile
inputs = cell(4, 1);

% Populate inputs in order indicated in the jobfile (generated by the batch manager)
inputs{1} = input_files; % Path to subject data
inputs{2} = DCM_files; % Paths to files with template DCMs
inputs{3} = out_dir_spec; % Path to GCM specification with model space
inputs{4} = PEB_file; % Path to estimated PEB file

%---------------------------------------------------------------------------------------
% STEP 3: Execute batch job
%---------------------------------------------------------------------------------------

% Call jobman with the jobfile and necessary inputs
spm_jobman('run', jobfile, inputs{:});

%---------------------------------------------------------------------------------------
% OUTPUTS
%---------------------------------------------------------------------------------------
% Running this batch job will produce the following output in the folder 'fits_dir'
% 1. BMA file called 'BMA_PEB_Full.mat' under fits_dir
%       This is the BMA obtained after taking a weighted average of the full and reduced
%       self models. Can be reviewed by calling spm_dcm_peb_review.
% NOTE: Running this job with a different reduced model for BMC will overwrite this BMA
% file. Please rename this file appropriately to prevent overwriting and to reflect the
% comparison carried out. Example: 'BMA_PEB_Full_vs_Self.mat'. This can be done manually
% or via the batch interface by selecting the relevant file operations module.

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
DCM.B{1} =  [
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
DCM.B{1} =  [
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
DCM.B{1} =  [
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
DCM.B{1} =  [
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
DCM.B{1} =  [
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
DCM.B{1} =  [
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
DCM.B{1} =  [
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
    imagesc(GCM{1, k}.B{1} + 0.75*fliplr(eye(Nareas)));
    colormap(gray)
    caxis([0, 1])
    title(sprintf('Model %02d', k))
    axis square
end

%---------------------------------------------------------------------------------------
% STEP 2: Load estimated PEB and perform BMR of model space
%---------------------------------------------------------------------------------------

% Load estimated PEB from file
load(fullfile(fits_dir, sprintf('PEB_%s.mat', name_tag)))

% Bayesian Model Reduction (BMR) and comparison of models
[BMA, BMR] = spm_dcm_peb_bmc(PEB, GCM);

% Save BMA and BMR
outfile = fullfile(fits_dir, 'BMA_BMR_Families.mat');
save(outfile, 'BMA', 'BMR')

%---------------------------------------------------------------------------------------
% STEP 3: Group models into families and compare
%---------------------------------------------------------------------------------------
% Now we partition the model space into families and perform comparisons at the level of
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
% Family 2 has greater evidence (~0.8)

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
% Family 2 has moderately greater evidence (~0.6)

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
% We demonstrate here the inclusion of covariates for 2nd-level/group inference with PEB.
% The dataset includes ages of participants, and while we do not anticipate any effect
% of age on modulation of connections due to faces, we illustrate specification of age
% as a covariate in the PEB design matrix for inference.

%---------------------------------------------------------------------------------------
% STEP 1: Prepare inputs for batch job 
%---------------------------------------------------------------------------------------

% Define covariates, and assign appropriate labels
PEB_name = 'Age';
covariate_name = 'Age';
covariate_values = [31, 25, 30, 26, 23, 26, 31, 26, 29, 23, 24, 24, 25, 24, 30, 25]';
%covariate_name = 'Sex';
%covariate_values = [0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0]'; % 1=Female

% Mean-center the covariate (Optional)
covariate_values = covariate_values - mean(covariate_values);

%---------------------------------------------------------------------------------------
% STEP 2: Setup jobfile and inputs 
%---------------------------------------------------------------------------------------

% Specify jobfile
jobfile = {fullfile(code_dir, 'batch_dcm_peb_covariate_job.m')}; 

% Initialize inputs to jobfile as a columnar cell array
% There are 4 inputs needed by this jobfile
inputs = cell(4, 1);

% Populate inputs in order indicated in the jobfile (generated by the batch manager)
inputs{1, 1} = PEB_name; % Specify / Estimate PEB: Name - cfg_entry
inputs{2, 1} = GCM_file; % Specify / Estimate PEB: DCMs - cfg_files
inputs{3, 1} = covariate_name; % Specify / Estimate PEB: Name - cfg_entry
inputs{4, 1} = covariate_values; % Specify / Estimate PEB: Value - cfg_entry

%---------------------------------------------------------------------------------------
% STEP 3: Execute batch job
%---------------------------------------------------------------------------------------

spm_jobman('run', jobfile, inputs{:});

% Greedy search for nested models on this PEB can be done to perform inference and
% identify which connections modulated by faces are affected by age.

%---------------------------------------------------------------------------------------
% OUTPUTS
%---------------------------------------------------------------------------------------
% Running this batch job will produce the following outputs in the folder 'fits_dir'
% 1. Estimated PEB file called 'PEB_Age.mat' under fits_dir
%       This is the group PEB estimated with age as a covariate.
