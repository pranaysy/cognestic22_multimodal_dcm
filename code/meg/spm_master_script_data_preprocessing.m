% SPM12 script to analyse the multi-subject, multi-modal human neuroimaging
% dataset described in Henson et al. (2018) Special Issue in Frontiers
%
% Note that you will need to have latest version of SPM12 on your MATLAB
% path, which you can download from here:
%       http://www.fil.ion.ucl.ac.uk/spm/
%
% plus the data, available from the OpenNeuro database in BIDS format: 
%      https://openneuro.org/datasets/ds000117.
%
% (A non-BIDS version is available here: 
%      ftp://ftp.mrc-cbu.cam.ac.uk/personal/rik.henson/wakemandg_hensonrn/
% but the BIDS-specific parts of code below will need changing)
%
% rik.henson@mrc-cbu.cam.ac.uk                              May 2018
% with help from Guillaume Flandin and Vladimir Litvak

clear;

% Add SPM12 to the MATLAB path
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
spm('asciiwelcome');

%% Input arguments
%==========================================================================
rawpth = '/imaging/henson/Wakeman/ds000117'; % Directory containing the raw data
outdir = '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset';
scrpth = fullfile(outdir,'code');                     % Directory containing the SPM analysis scripts
outpth = fullfile(outdir, 'data', 'derivatives',spm('Ver'));   % Output directory

keepdata   = false; % If false, intermediate files will be deleted to save disk space

% numworkers = 16; % Number of workers for distributed computing
% if numworkers, parpool(numworkers); end

delete(gcp('nocreate'))
numworkers = 16;
P=cbupool(numworkers, '--mem-per-cpu=4G --time=72:00:00 --nodelist=node-j10');
parpool(P, P.NumWorkers);

timewin    = [-100 500];

%% Parse BIDS-formatted dataset
%==========================================================================
BIDS   = spm_BIDS(rawpth);

subs   = spm_BIDS(BIDS,'subjects', 'task','facerecognition');
nsub   = numel(subs);
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);

%% Prepare output directory
%==========================================================================

fprintf('%-40s: %30s', 'Copy files in derivatives','...');              %-#

%-Create output directory tree if necessary
%--------------------------------------------------------------------------
spm_mkdir(outpth,subdir,{'meg','anat','fmri'});

%-Pipeline description
%--------------------------------------------------------------------------
spm_jsonwrite(fullfile(outpth,'pipeline_description.json'),struct(...
    'Name',spm('Ver'),...
    'Version',spm('Version'),...
    'CodeURL','http://www.fil.ion.ucl.ac.uk/spm/',...
    'License','Creative Commons Attribution 4.0 International Public License'),...
    struct('indent','  '));

%-Copy FIF files
%--------------------------------------------------------------------------
for s = 1:nsub
    runs = spm_BIDS(BIDS,'runs', 'sub',subs{s}, 'modality','meg', 'type','meg');
    for r = 1:numel(runs)
        f = fullfile(rawpth,'derivatives','meg_derivatives',subdir{s},'ses-meg','meg',[subdir{s} '_ses-meg_task-facerecognition_run-' runs{r} '_proc-sss_meg.fif']);
        spm_copy(f, fullfile(outpth,subdir{s},'meg'));
    end
end

%-Copy and gunzip T1 MPRAGE images
%--------------------------------------------------------------------------
for s = 1:nsub
    f = spm_BIDS(BIDS,'data','sub',subs{s},'modality','anat','type','T1w','acq','mprage');
    spm_copy(f, fullfile(outpth,subdir{s},'anat'), 'gunzip',true);
end

fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#



%% Preprocessing M/EEG data 
%==========================================================================
spm_jobman('initcfg');
spm('defaults','EEG');
spm_get_defaults('cmdline',true);

tic
parfor (s = 1:nsub, numworkers)
    if numworkers > 0
        spm_jobman('initcfg');
        spm('defaults', 'EEG');
        spm_get_defaults('cmdline',true);
    end
    
    runs = spm_BIDS(BIDS,'runs', 'sub',subs{s}, 'modality','meg', 'type','meg');
    
    % Change to subject's directory
    cd(fullfile(outpth,subdir{s},'meg'));
    
    Dc = {};
    
    %% Convert, update channels and event times, downsample
    for r = 1:numel(runs)
         
        S = [];
        %S.dataset = char(spm_BIDS(BIDS,'data','sub',subs{s},'type','meg','run', runs{r}, 'proc', 'sss'));
        % The sss maxfiltered .fif data files are in derivatives/meg_derivatives
        S.dataset = fullfile(outpth,subdir{s},'meg',[subdir{s} '_ses-meg_task-facerecognition_run-' runs{r} '_proc-sss_meg.fif']);
        S.mode = 'continuous';
        S.channels = {'EEG', 'MEGMAG', 'MEGPLANAR'};
        S.checkboundary = 1;
        D = spm_eeg_convert(S);
        
        % Set channel types and bad channels
        S = [];
        S.D    = D;
        S.task = 'bidschantype';
        S.save = 1;
        S.filename = fullfile(rawpth,subdir{s},'ses-meg',[subdir{s} '_ses-meg_task-facerecognition_channels.tsv']);
        D = spm_eeg_prep(S);
             
        S.D    = D;
        S.task = 'bidschanstatus';
        D = spm_eeg_prep(S);
        
        % Load events stored in BIDS
        S = [];
        S.D        = D;
        S.task     = 'loadbidsevents';
        S.replace  = 1;
        S.filename = char(spm_BIDS(BIDS,'data','ses','meg','sub',subs{s},'run', runs{r},'type','events'));
        D = spm_eeg_prep(S);
        
        % Do NOT downsample the data for DCM
%         S = [];
%         S.D = D;
%         S.method = 'resample';
%         S.fsample_new = 200;
%         D = spm_eeg_downsample(S);
%         
%         if ~keepdata, delete(S.D); end
        
        Dc{r} = D.save();
    end

    
    De = {};
    for r = 1:numel(runs)
        
% Baseline correction added during epoching should handle highpass
%         S = [];
%         S.D = Dc{r};
%         S.type = 'butterworth';
%         S.band = 'high';
%         S.freq = 1;
%         S.dir = 'twopass';
%         S.order = 5;
%         S.prefix = 'f';
%         D = spm_eeg_filter(S);
        
        S = [];
        S.D = Dc{r};
        S.type = 'butterworth';
        S.band = 'low';
        S.freq = 40;
        S.dir = 'twopass';
        S.order = 5;
        S.prefix = 'f';
        D = spm_eeg_filter(S);
        
%        if ~keepdata, delete(S.D); end
               
        S = [];
        S.D = D;
        S.timewin = timewin;
        S.trialdef(1).conditionlabel = 'Famous';
        S.trialdef(1).eventtype  = 'BIDS';
        S.trialdef(1).eventvalue = 'Famous';
        S.trialdef(2).conditionlabel = 'Unfamiliar';
        S.trialdef(2).eventtype  = 'BIDS';
        S.trialdef(2).eventvalue = 'Unfamiliar';
        S.trialdef(3).conditionlabel = 'Scrambled';
        S.trialdef(3).eventtype  = 'BIDS';
        S.trialdef(3).eventvalue = 'Scrambled';
        S.bc = 1; % Original script does not have baseline correction
        S.prefix = 'e';
        S.eventpadding = 0;
        D = spm_eeg_epochs(S);
        
        if ~keepdata, delete(S.D); end
        
        De{r} = D;
    end
    
    % Merge across runs
    S = [];
    S.D = De;
    S.prefix = 'c';
    D = spm_eeg_merge(S);
    
    if ~keepdata
        for i = 1:numel(De)
            delete(De{i});
        end
    end
    
    % Define average reference montage excluding bad channels
    eegchan = D.indchantype('EEG');
    goodind = D.indchantype('EEG', 'GOOD');
    
    goodind = find(ismember(eegchan, goodind));
    
    tra              =  eye(length(eegchan));
    tra(: ,goodind)  =  tra(:, goodind) - 1/length(goodind);
    
    montage          = [];
    montage.labelorg = D.chanlabels(eegchan);
    montage.labelnew = D.chanlabels(eegchan);
    montage.tra      = tra;
    
    % Apply montage
    S = [];
    S.D = D;
    S.mode = 'write';
    S.prefix = 'M';
    S.montage = montage;
    S.keepothers = 1;
    S.keepsensors = 1;
    D = spm_eeg_montage(S);
    
    if ~keepdata, delete(S.D); end
    
    S = [];
    S.D = D;
    S.mode = 'reject';
    S.badchanthresh = 0.2;
    S.methods.channels = {'EOG'};
    S.methods.fun = 'threshchan';
    S.methods.settings.threshold = 200;
    S.methods.settings.excwin = 1000;
    S.append = true;
    S.prefix = 'a';
    D = spm_eeg_artefact(S); 
    
    if ~keepdata, delete(S.D); end
    
    % Set condition order 
    D = conditions(D, [1 2 3], {'Famous',  'Unfamiliar', 'Scrambled'});
    D = condlist(D, {'Famous',  'Unfamiliar', 'Scrambled'});
    D.save;
    
    % Average
    S = [];
    S.D = D;
    S.robust.bycondition = 1;
    S.circularise = false;
    S.prefix = 'm';
    D = spm_eeg_average(S);
    
    if ~keepdata, delete(S.D); end
    
    % Set condition order 
    D = conditions(D, [1 2 3], {'Famous',  'Unfamiliar', 'Scrambled'});
    D = condlist(D, {'Famous',  'Unfamiliar', 'Scrambled'});
    D.save;

    % Contrast conditions
    S = [];
    S.D = D;
    S.c = [0 0 1; 0.5 0.5 0];
    S.label = {'Scrambled', 'Faces'}';
    S.weighted = 0;
    S.prefix = 'w';
    D = spm_eeg_contrast(S);   
    
   
    %% Source analysis (create forward model)
    mrifidfile = fullfile(rawpth,subdir{s},'ses-mri','anat', ['sub-' subs{s} '_ses-mri_acq-mprage_T1w.json']);

    jobfile = {fullfile(scrpth, 'meg', 'batch_forward_model_job.m')}; 
    inputs = cell(3,1);
    inputs{1} = {fullfile(D)};
    inputs{2} = {spm_select('FPList',fullfile(outpth,subdir{s},'anat'),'^sub-.*_T1w\.nii$')};
    inputs{3} = {mrifidfile};

    spm_jobman('run', jobfile, inputs{:});
    
    D = reload(D);
    [L, Dl] = spm_eeg_lgainmat(D); % Compute gain matrix and store
    
end
checkp1=toc;
