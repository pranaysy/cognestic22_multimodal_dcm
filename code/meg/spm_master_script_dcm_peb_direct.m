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
base_dir = '/imaging/henson/Wakeman/cognestic22_multimodal_dcm'; % Change this to yours
addpath(genpath(fullfile(base_dir, 'code'))) % Add scripts & functions to workspace

% All fits go in this directory
fits_dir = fullfile(base_dir, 'fits', 'direct_script', 'meg');

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

% Save full model
DCM_Full = DCM; % Keep DCM in memory
dcm_full_file = fullfile(fits_dir, 'templates', 'DCMs', strcat(DCM.name, '.mat'));
save(dcm_full_file, 'DCM_Full')

%--------------------------------------------------------------------------
% STEP 5: Specify reduced models, if any
%--------------------------------------------------------------------------

% Reduced model with modulation of only self connections
DCM.name = 'DCM_Self';
DCM.B{1} = self_connections;

% Save reduced model
dcm_self_file = fullfile(fits_dir, 'templates', 'DCMs', strcat(DCM.name, '.mat'));
save(dcm_self_file, 'DCM')
DCM_Self = DCM; % Keep DCM in memory

%---------------------------------------------------------------------------------------
% STEP 6: Replicate DCM specification across subjects
%---------------------------------------------------------------------------------------

% Populate list of processed files as a column-order cell (N-files × 1)
% These files should contain forward models (with or without gain matrices)
files = dir(fullfile(base_dir, 'data', '**', 'maM*.mat'));
input_files = arrayfun(@(x) fullfile(x.folder, x.name), files, 'UniformOutput', false);

% Generate GCM with rows corresponding to DCMs for all subjects, columns representing
% models per subject (we only have one column here for the full model)
GCM = {};
for f=1:length(input_files)    
    GCM{f, 1} = DCM_Full; % Full model
    GCM{f, 1}.xY.Dfile = input_files{f}; % Add path to subject f's data file
    GCM{f, 1}.name = sprintf('%s_sub-%02d', name, f) % Add subject identifier
end   

% Save GCM specification
save(fullfile(fits_dir, 'templates', 'GCMs', 'Full', 'GCM_DCM_Full.mat'), 'GCM')

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

%---------------------------------------------------------------------------------------
% STEP 1: Manage environment
%---------------------------------------------------------------------------------------

% Initialize Parallel Compute Pool (Example Instructions for CBU Cluster)
delete(gcp('nocreate')) % Shut down any existing pool
n_workers = length(input_files);
P=cbupool(n_workers, '--mem-per-cpu=4G --time=12:00:00 --nodelist=node-j10');
parpool(P, P.NumWorkers);

% During fitting, individual subject-level DCMs will be estimated in the working folder
% We'll navigate to the templates GCMs folder to be consistent with the 'batch' script
cd(fullfile(fits_dir, 'templates', 'GCMs', 'Full'))

%---------------------------------------------------------------------------------------
% STEP 2: Fit DCMs in parallel
%---------------------------------------------------------------------------------------

% Fit GCM
GCM = spm_dcm_fit(GCM(:, 1), true);

% Switch back to base directory
cd(base_dir)

% Save fitted GCM
save(fullfile(fits_dir, 'GCM_Full'), 'GCM')

%---------------------------------------------------------------------------------------
% STEP 3: Estimate 2nd-level PEB model
%---------------------------------------------------------------------------------------
M = []; M.Q = 'all'; % Random effects over all parameters
PEB = spm_dcm_peb(GCM, M, {'B'});

% Save estimated PEB
save(fullfile(fits_dir, 'PEB_Full'), 'PEB')

%---------------------------------------------------------------------------------------
% STEP 4: Review estimated PEB
%---------------------------------------------------------------------------------------
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

% Construct model space
modelspace = {};
for k=1:size(GCM, 2)
    model = GCM{1, k};
    if isfield(model, 'M')
        model = rmfield(model, 'M'); % Fitted model creates errors
    end
    modelspace{k} = model;
end

% Save model space
save(fullfile(fits_dir,'templates', 'GCMs', 'Full_vs_Self', 'GCM_Full_vs_Self.mat'), 'modelspace')

[BMA, BMR] = spm_dcm_peb_bmc(PEB, modelspace);

% Write to disk
save('BMC_PEB_Full_vs_Self', 'BMA', 'BMR')

%---------------------------------------------------------------------------------------
% STEP 2: Review estimated BMA
%---------------------------------------------------------------------------------------
spm_dcm_peb_review(BMA, GCMp)

%---------------------------------------------------------------------------------------
%
%     888b     d888                               888888b.   888b     d888  .d8888b.  
%     8888b   d8888                               888  "88b  8888b   d8888 d88P  Y88b 
%     88888b.d88888                               888  .88P  88888b.d88888 888    888 
%     888Y88888P888  .d88b.  888d888 .d88b.       8888888K.  888Y88888P888 888        
%     888 Y888P 888 d88""88b 888P"  d8P  Y8b      888  "Y88b 888 Y888P 888 888        
%     888  Y8P  888 888  888 888    88888888      888    888 888  Y8P  888 888    888 
%     888   "   888 Y88..88P 888    Y8b.          888   d88P 888   "   888 Y88b  d88P 
%     888       888  "Y88P"  888     "Y8888       8888888P"  888       888  "Y8888P"
%
%---------------------------------------------------------------------------------------
% Since we observed significant modulation of between-region connections, we can zoom
% in and test for modulation of selective groups of connections like forward, backward
% and lateral connections. We demo this sequential hypothesis testing here using BMC.

%---------------------------------------------------------------------------------------
%
%     8888888888                                                  888 
%     888                                                         888 
%     888                                                         888 
%     8888888  .d88b.  888d888 888  888  888  8888b.  888d888 .d88888 
%     888     d88""88b 888P"   888  888  888     "88b 888P"  d88" 888 
%     888     888  888 888     888  888  888 .d888888 888    888  888 
%     888     Y88..88P 888     Y88b 888 d88P 888  888 888    Y88b 888 
%     888      "Y88P"  888      "Y8888888P"  "Y888888 888     "Y88888
%
%---------------------------------------------------------------------------------------
% Test for Modulation of Forward Connections

%---------------------------------------------------------------------------------------
% STEP 1: Define model space
%---------------------------------------------------------------------------------------

% Get full DCM specification
if isfield(DCM_Full, 'M')
    DCM_Full = rmfield(DCM_Full, 'M');
end
DCM = DCM_Full;

% Switch off Forward connections in B-matrix
DCM.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  1    1    1    0  ];   % lOFA
    [  1    1    0    1  ];   % rOFA
    [  0    0    1    1  ];   % lFFA
    [  0    0    1    1  ];   % rFFA
];

% Model Space
models = {DCM_Full, DCM};

%---------------------------------------------------------------------------------------
% STEP 2: BMR and BMC
%---------------------------------------------------------------------------------------
[BMA, BMR] = spm_dcm_peb_bmc(PEB, models);

% Write to disk
save('BMC_PEB_Forward', 'BMA', 'BMR')

%---------------------------------------------------------------------------------------
% STEP 3: Review
%---------------------------------------------------------------------------------------
spm_dcm_peb_review(BMA, GCMp)

% No forward connections from OFA to FFA are modulated by Faces

%---------------------------------------------------------------------------------------
%
%     888888b.                     888                                         888 
%     888  "88b                    888                                         888 
%     888  .88P                    888                                         888 
%     8888888K.   8888b.   .d8888b 888  888 888  888  888  8888b.  888d888 .d88888 
%     888  "Y88b     "88b d88P"    888 .88P 888  888  888     "88b 888P"  d88" 888 
%     888    888 .d888888 888      888888K  888  888  888 .d888888 888    888  888 
%     888   d88P 888  888 Y88b.    888 "88b Y88b 888 d88P 888  888 888    Y88b 888 
%     8888888P"  "Y888888  "Y8888P 888  888  "Y8888888P"  "Y888888 888     "Y88888
%
%---------------------------------------------------------------------------------------
% Test for Modulation of Backward Connections

%---------------------------------------------------------------------------------------
% STEP 1: Define model space
%---------------------------------------------------------------------------------------

% Get full DCM specification
if isfield(DCM_Full, 'M')
    DCM_Full = rmfield(DCM_Full, 'M');
end
DCM = DCM_Full;

% Switch off Backward connections in B-matrix
DCM.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  1    1    0    0  ];   % lOFA
    [  1    1    0    0  ];   % rOFA
    [  1    0    1    1  ];   % lFFA
    [  0    1    1    1  ];   % rFFA
];

% Model Space
models = {DCM_Full, DCM};

%---------------------------------------------------------------------------------------
% STEP 2: BMR and BMC
%---------------------------------------------------------------------------------------
[BMA, BMR] = spm_dcm_peb_bmc(PEB, models);

% Write to disk
save('BMC_PEB_Backward', 'BMA', 'BMR')

%---------------------------------------------------------------------------------------
% STEP 3: Review
%---------------------------------------------------------------------------------------
spm_dcm_peb_review(BMA, GCMp)

% Both backward connections from FFA to OFA are modulated by Faces

%---------------------------------------------------------------------------------------
%
%     888               888                            888 
%     888               888                            888 
%     888               888                            888 
%     888       8888b.  888888 .d88b.  888d888 8888b.  888 
%     888          "88b 888   d8P  Y8b 888P"      "88b 888 
%     888      .d888888 888   88888888 888    .d888888 888 
%     888      888  888 Y88b. Y8b.     888    888  888 888 
%     88888888 "Y888888  "Y888 "Y8888  888    "Y888888 888
%
%---------------------------------------------------------------------------------------
% Test for Modulation of Lateral Connections

%---------------------------------------------------------------------------------------
% STEP 1: Define model space
%---------------------------------------------------------------------------------------

% Get full DCM specification
if isfield(DCM_Full, 'M')
    DCM_Full = rmfield(DCM_Full, 'M');
end
DCM = DCM_Full;

% Switch off Lateral connections in B-matrix
DCM.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  1    0    1    0  ];   % lOFA
    [  0    1    0    1  ];   % rOFA
    [  1    0    1    0  ];   % lFFA
    [  0    1    0    1  ];   % rFFA
];

% Model Space
models = {DCM_Full, DCM};

%---------------------------------------------------------------------------------------
% STEP 2: BMR and BMC
%---------------------------------------------------------------------------------------
[BMA, BMR] = spm_dcm_peb_bmc(PEB, models);

% Write to disk
save('BMC_PEB_Lateral', 'BMA', 'BMR')

%---------------------------------------------------------------------------------------
% STEP 3: Review
%---------------------------------------------------------------------------------------
spm_dcm_peb_review(BMA, GCMp)

% Bidirectional lateral connections between both OFA and FFA are modulated by Faces
% Follow-up: Test for OFA and FFA separately?

%---------------------------------------------------------------------------------------
%
%      .d8888b.           888  .d888 
%     d88P  Y88b          888 d88P"  
%     Y88b.               888 888    
%      "Y888b.    .d88b.  888 888888 
%         "Y88b. d8P  Y8b 888 888    
%           "888 88888888 888 888    
%     Y88b  d88P Y8b.     888 888    
%      "Y8888P"   "Y8888  888 888
%
%---------------------------------------------------------------------------------------
% After testing which between-region connections are modulated, we test which self
% connections are modulated by faces.

%---------------------------------------------------------------------------------------
% STEP 1: Define model space
%---------------------------------------------------------------------------------------

% Get full DCM specification
if isfield(DCM_Full, 'M')
    DCM_Full = rmfield(DCM_Full, 'M');
end
DCM = DCM_Full;

% Switch off Self connections in B-matrix
DCM.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  0    1    1    0  ];   % lOFA
    [  1    0    0    1  ];   % rOFA
    [  1    0    0    1  ];   % lFFA
    [  0    1    1    0  ];   % rFFA
];

% Model Space
models = {DCM_Full, DCM};

%---------------------------------------------------------------------------------------
% STEP 2: BMR and BMC
%---------------------------------------------------------------------------------------
[BMA, BMR] = spm_dcm_peb_bmc(PEB, models);

% Write to disk
save('BMC_PEB_Self', 'BMA', 'BMR')

%---------------------------------------------------------------------------------------
% STEP 3: Review
%---------------------------------------------------------------------------------------
spm_dcm_peb_review(BMA, GCMp)

% No self connections in OFA and FFA are modulated by Faces

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
% Instead of testing for modulation of all self-connections versus none, we could ask
% whether any self-connection is being modulated at all. This involves multiple models
% with different combinations of self-connections being modulated. The model space
% therefore expands from just 2 models, like the previous BMC, to a much larger number
% of models. As an example here, we consider bilateral pairs of OFA and FFA, and test
% whether at least one self-connection was modulated by faces. These leads to a model
% space consisting of these four models:
%       1. Modulation of both OFA and FFA self-connections
%       2. Modulation of only OFA self-connections
%       3. Modulation of only FFA self-connections
%       4. Modulation of no self-connections
% The first three models encapsulate the hypothesis: are any self-connections modulated
% by faces? While the fourth model corresponds to the alternate hypothesis that no
% self-connections are modulated. Accordingly, we group these models into two families,
% and perform inference at the level of model families, rather than at the level of the
% individual models.

%---------------------------------------------------------------------------------------
% STEP 1: Setup model space and families
%---------------------------------------------------------------------------------------

% >> Family 1: Atleast one self-connection
% > Model 1 is the full model with both OFA and FFA self-connections
if isfield(DCM_Full, 'M')
    DCM_Full = rmfield(DCM_Full, 'M');
end
DCM_m1 = DCM_Full;

% > Model 2 has OFA self-connections only
% Get full DCM specification
DCM_m2 = DCM_m1;

% Switch off FFA self-connections in B-matrix
DCM_m2.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  1    1    1    0  ];   % lOFA
    [  1    1    0    1  ];   % rOFA
    [  1    0    0    1  ];   % lFFA
    [  0    1    1    0  ];   % rFFA
];

% > Model 3 has FFA self-connections only
% Get full DCM specification
DCM_m3 = DCM_m1;

% Switch off OFA self-connections in B-matrix
DCM_m3.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  0    1    1    0  ];   % lOFA
    [  1    0    0    1  ];   % rOFA
    [  1    0    1    1  ];   % lFFA
    [  0    1    1    1  ];   % rFFA
];

% >> Family 2: No self-connections
% > Model 4 has no self-connections
% Get full DCM specification
DCM_m4 = DCM_m1;

% Switch off all self-connections in B-matrix
DCM_m4.B{1} =  [
%    lOFA rOFA lFFA rFFA
    [  0    1    1    0  ];   % lOFA
    [  1    0    0    1  ];   % rOFA
    [  1    0    0    1  ];   % lFFA
    [  0    1    1    0  ];   % rFFA
];

% Define family-wise model space
models = {DCM_m1, DCM_m2, DCM_m3, DCM_m4};
families = [1, 1, 1, 2];

%---------------------------------------------------------------------------------------
% STEP 2: Perform BMR of model space nested under PEB
%---------------------------------------------------------------------------------------

% Bayesian Model Reduction (BMR) and comparison of models
[BMA, BMR] = spm_dcm_peb_bmc(PEB, models);
% Best model appears to be the one with modulation of only FFA self-connections
% This contradicts the previous BMC, where no self connections were modulated by Faces

% Write to disk
save('BMC_PEB_Self_ExpandedModelSpace', 'BMA', 'BMR')

% We now group models under families & consider each family as equally likely. This lets
% us pool model evidence within each family, and instead of picking a winning model, we
% perform inference about characteristics that define models grouped under that family.
% The characteristic here being modulation of self-connections, so models with any self
% connections i.e. models 1, 2 and 3, belong to the same family.

%---------------------------------------------------------------------------------------
% STEP 3: Compare families of models
%---------------------------------------------------------------------------------------

% Bayesian Model Selection (BMS) and averaging over families of models
[BMAf, fam] = spm_dcm_peb_bmc_fam(BMA, BMR, families, 'NONE');

% Write to disk
save('BMC_PEB_Self_ExpandedModelSpace_Families', 'BMAf', 'fam')

% The family of models with modulation of at least one self-connection has marginally
% higher posterior probability than the family with no modulation of self-connections.
% There is not enough evidence in favor of modulation of self-connections by faces.

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
covariate_name = 'Age';
covariate_values = [31, 25, 30, 26, 23, 26, 31, 26, 29, 23, 24, 24, 25, 24, 30, 25]';
%covariate_name = 'Sex';
%covariate_values = [0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0]'; % 1=Female

% Mean-center the covariate (Optional)
covariate_values = covariate_values - mean(covariate_values);

% Design Matrix
M.X = [ones([length(input_files), 1]), covariate_values]; % First covariate is group mean
M.Xnames = {'Commonalities', 'Age'};
M.Q = 'all'; % Random effects over all parameters

%---------------------------------------------------------------------------------------
% STEP 2: Fit 2nd-level PEB
%---------------------------------------------------------------------------------------
[PEB, GCMp] = spm_dcm_peb(GCMf, M, {'B'});

% Write to disk
save('PEB_Covariate_Fit', 'PEB', 'GCMp')
%---------------------------------------------------------------------------------------
% STEP 3: Review
%---------------------------------------------------------------------------------------
spm_dcm_peb_review(PEB, GCMp)

% Age has no effect on modulation of connections due to faces except for a small effect
% on rFFA self-connection.
