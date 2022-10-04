% SPM12 script to run DCM for fMRI on the multi-subject, multi-modal human neuroimaging
% dataset 
%
% Based on Henson et al. (2018) Special Issue in Frontiers and Zeidman ***
%
% Note that you will need to have latest version of SPM12 on your MATLAB
% path, which you can download from here:
%       http://www.fil.ion.ucl.ac.uk/spm/
%
% plus the data, available from the OpenNeuro database in BIDS format: 
%      https://openneuro.org/datasets/ds000117
%
% You will also need to have run the preprocessing steps described in
% Appendix 2 of https://www.frontiersin.org/articles/10.3389/fnins.2019.00300/full#supplementary-material
%
% rik.henson@mrc-cbu.cam.ac.uk                              Aug 2022
% with help from Pranay Yadav


% TO DO: drop "stat_concat"


clear;

% SPM12PATH = '/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'
% addpath(SPM12PATH);
spm fmri



%% Input arguments
%==========================================================================
rawpth = '/imaging/henson/users/rh01/Methods/DanData/OpenNeuro/ds000117' % where put raw BIDS data
derpth = fullfile(rawpth,'derivatives','SPM12') % will exist if you have run preprocessing scripts for Henson et al, 2019, or else download from ****

outpth = fullfile(derpth,'DCM');
try mkdir(outpth); end; 

% Create sub-directory for scripts,
scrpth = fullfile(outpth,'code');                 
try mkdir(scrpth); end; 
% ... and download files to here from ****

BIDS   = spm_BIDS(rawpth);

subs   = spm_BIDS(BIDS,'subjects', 'task','facerecognition');
matched_subs = {};
for s = [1:9 11:16]; % sub-10 only has 170 volumes in last run so cannot be combined in GCM below
    matched_subs{end+1} = subs{s};
end
subs = matched_subs;
nsub   = numel(subs)
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);

trialtypes = {'Famous','Unfamiliar','Scrambled'}; 

numworkers = nsub; % Number of workers for distributed computing (depends on system)
if numworkers > 0
    delete(gcp('nocreate')) % Shut down any existing pool
%     P=cbupool(numworkers, '--mem-per-cpu=4G --time=12:00:00 --nodelist=node-j10');
%     parpool(P, P.NumWorkers);
    parpool(numworkers);
end

runs = spm_BIDS(BIDS,'runs', 'modality','func', 'type','bold', 'task','facerecognition'); % across subjects
nrun = numel(runs);

nscan = repmat(208,1,nrun); % sub-15, run-01 has 209 scans, so ignore last

TR = 2;

%% 1. Combine conditions and the concatenate runs for DCM
%==========================================================================

conds = [];
conds.names{1} = 'All'; 
conds.names{2} = 'Faces';    % Famous and Nonfamous

% Create concatenated trial definition and movement parameter files
%----------------------------------------------------------------------
cd(outpth)

for s = 1:nsub    
    try mkdir(subdir{s}); end; 
    outdir = fullfile(outpth,subdir{s},'stat_concat')
    try mkdir(outdir); end;
    
    for t = 1:2
        conds.durations{t} = 0;
        conds.onsets{t} = [];
    end
    
    time_so_far = 0; volfiles = {}; movepar = [];
    for r = 1:length(runs)       
%         volfiles{r} = spm_select('ExtFPList',fullfile(derpth,subdir{s},'func'),sprintf('^swsub-.*run-%s_bold\\.nii$',runs{r}),[1 Inf]);
%         nscan(s,r)    = size(volfiles{r},1);
        volfiles{r} = spm_select('ExtFPList',fullfile(derpth,subdir{s},'func'),sprintf('^swsub-.*run-%s_bold\\.nii$',runs{r}),[1:nscan(r)]);
   
        trlfile = fullfile(derpth,subdir{s},'func',sprintf('sub-%s_run-%s_spmdef.mat',subs{s},runs{r}));
        d = load(trlfile);
%         % Trial definition files should exist if run preprocessing above, but in case not...
%         if ~exist(trlfile,'file')
%             d = spm_load(char(spm_BIDS(BIDS,'data','modality','func','type','events','sub',subs{s},'run',runs{r})));
%             clear conds
%             for t = 1:numel(trialtypes)
%                 conds.names{t} = trialtypes{t};
%                 conds.durations{t} = 0;
%                 conds.onsets{t} = d.onset(strcmpi(d.stim_type,trialtypes{t}));
%             end
%             save(fullfile(derpth,subdir{s},'func',sprintf('sub-%s_run-%s_spmdef.mat',subs{s},runs{r})),'-struct','conds');           
%         end       
%         for t = 1:numel(trialtypes)
%             conds.onsets{t} = [conds.onsets{t}; d.onsets{t} + time_so_far];
%         end
        conds.onsets{1} = [conds.onsets{1}; sort([d.onsets{1}; d.onsets{2}; d.onsets{3}]) + time_so_far];
        conds.onsets{2} = [conds.onsets{2}; sort([d.onsets{1}; d.onsets{2}]) + time_so_far];
 
%        time_so_far = time_so_far + nscan(s,r)*TR;
        time_so_far = time_so_far + nscan(r)*TR;
        
        d = load(spm_select('FPList',fullfile(derpth,subdir{s},'func'),sprintf('^rp.*run-%s.*\\.txt$',runs{r})));
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
    jobfile = {fullfile(scrpth,'batch_stats_fmri_concatenated_specify_job.m')};

    outdir = fullfile(outpth,subdir{s},'stat_concat')
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
%    spm_fmri_concatenate(spmfile, nscan(s,:)); % Assumes nscan already in memory from above (could save instead)
    spm_fmri_concatenate(spmfile, nscan); % Assumes nscan already in memory from above (could save instead)
    delete('SPM_backup.mat');
    
    % Estimate new model
    %----------------------------------------------------------------------    
    jobfile = {fullfile(scrpth,'batch_stats_fmri_concatenated_estimate_job.m')};
    inputs  = {};
    inputs{1} = cellstr(spmfile);
    spm_jobman('run', jobfile, inputs{:});
end

%% 2. Create VOI files (ie timeseries for ROIs)
%==========================================================================

ROI_names = {'lOFA','rOFA','lFFA','rFFA'};
ROI_coord = {[-38, -86, -14],[+36, -86, -10],[-42, -56, -20],[+42, -52, -14]};

parfor (s = 1:nsub, numworkers)
    if numworkers > 0
        spm_jobman('initcfg');
        spm('defaults', 'fmri');
        spm_get_defaults('cmdline',true);
    end
    
    outdir = fullfile(outpth,subdir{s},'stat_concat')
    cd(outdir)
    spmfile = fullfile(outdir,'SPM.mat');
    
    jobfile = {fullfile(scrpth,'batch_VOI_job.m')};  
    for r = 1:numel(ROI_names)
        inputs  = {};
        inputs{1} = cellstr(spmfile);
        inputs{2} = ROI_names{r};
        inputs{3} = ROI_coord{r};
        spm_jobman('run', jobfile, inputs{:});
    end
end

%% 3. Specify DCM for sub-01
%==========================================================================

ref_sub = 1; % eg first subject
outdir = fullfile(outpth,subdir{ref_sub},'stat_concat');
load(fullfile(outdir,'SPM.mat'));
        
DCM = [];

for r = 1:numel(ROI_names)
    load(fullfile(outdir,sprintf('VOI_%s_1.mat',ROI_names{r})),'xY');
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

DCM.options.nonlinear  = 0;
DCM.options.two_state  = 0;
DCM.options.stochastic = 0;
DCM.options.centre     = 0; % centre regressors (default = 0)
DCM.options.nograph    = 1;
DCM.options.maxnodes   = 8;
DCM.options.maxit      = 128;
DCM.options.hidden     = [];
DCM.options.induced    = 0;
        
DCM.M.options = struct(); % needed else crashes in estimation

% Full model
%----------------------------------------------------------------------

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
DCM.b(:,:,1) = zeros(DCM.n,DCM.n);
DCM.b(:,:,2) = [
     1     1     1     0
     1     1     0     1
     1     0     1     1
     0     1     1     1
];
DCM.d        = zeros(DCM.n,DCM.n,0); % needed else crashes in estimation

outfile = fullfile(outdir,'DCM_full.mat')
save(outfile,'DCM');

%DCM = spm_dcm_estimate(DCM)
%save(fullfile(outdir,'DCM_full_estimated.mat'),'DCM');
          
% Self-only (modulation) model
%----------------------------------------------------------------------

DCM.b(:,:,2) = [
     1     0     0     0
     0     1     0     0
     0     0     1     0
     0     0     0     1
];

outfile = fullfile(outdir,'DCM_self.mat')
save(outfile,'DCM');

% % Between-only (modulation) model
% %----------------------------------------------------------------------
% 
% DCM.b(:,:,2) = [
%      0     1     1     0
%      1     0     0     1
%      1     0     0     1
%      0     1     1     0
% ];
% 
% outfile = fullfile(outdir,'DCM_betw.mat')
% save(outfile,'DCM');
% 
% % Null (modulation) model
% %----------------------------------------------------------------------
% 
% DCM.b(:,:,2) = [
%      0     0     0     0
%      0     0     0     0
%      0     0     0     0
%      0     0     0     0
% ];
% 
% outfile = fullfile(outdir,'DCM_null.mat')
% save(outfile,'DCM');

%% 4. Specify GCM (DCM for each subject and model)
%==========================================================================
% (value of "ref_sub" must match value that generated DCMs above, eg 1)
GCMname = 'full_self';

jobfile = {fullfile(scrpth,'batch_create_gcm_job.m')}; 
inputs  = {};
inputs{1} = cellstr(outpth);
inputs{2} = GCMname;
inputs{3} = cellstr(fullfile(outpth,subdir{ref_sub},'stat_concat','DCM_full.mat'));
inputs{4} = {fullfile(outpth,subdir{ref_sub},'stat_concat','DCM_self.mat')};
% inputs{4} = {fullfile(outpth,subdir{ref_sub},'stat_concat','DCM_self.mat');...
%              fullfile(outpth,subdir{ref_sub},'stat_concat','DCM_betw.mat');...
%              fullfile(outpth,subdir{ref_sub},'stat_concat','DCM_null.mat')};
for s = 1:nsub
    inputs{5}{s,1} = fullfile(outpth,subdir{s},'stat_concat','SPM.mat');
end
for r = 1:numel(ROI_names)
    files = {};
    for s = 1:nsub
        inputs{5+r}{s,1} = fullfile(outpth,subdir{s},'stat_concat',sprintf('VOI_%s_1.mat',ROI_names{r}));
    end
end
spm_jobman('run', jobfile, inputs{:});

%% 5. Estimate GCM (may take time)
%==========================================================================
% addpath to get updated version of spm_dcm_fit where parfor==1
addpath /imaging/henson/users/rh01/Methods/DanData/OpenNeuro/ds000117/derivatives/SPM12/DCM/code 

jobfile = {fullfile(scrpth,'batch_fit_gcm_job.m')}; 
inputs  = {};
inputs{1} = cellstr(fullfile(outpth,GCMname));
inputs{2} = cellstr(outpth);
inputs{3} = sprintf(sprintf('GCM_%s_fit',GCMname));
inputs{4} = 2; % 2 for iterative PEB fit
spm_jobman('run', jobfile, inputs{:});


%% 6. Run PEB
%==========================================================================
jobfile = {fullfile(scrpth,'batch_fit_peb_job.m')}; 
inputs  = {};
inputs{1} = sprintf('%s_fit',GCMname);
inputs{2} = cellstr(fullfile(outpth,sprintf('GCM_%s_fit.mat',GCMname)));
inputs{3} = 1;
spm_jobman('run', jobfile, inputs{:});


%% 7. Run BMC across all possible models (BMR)
%==========================================================================
jobfile = {fullfile(scrpth,'batch_peb_bmr_search_job.m')}; 
inputs  = {};
inputs{1} = cellstr(fullfile(outpth,sprintf('PEB_%s_fit.mat',GCMname)));
inputs{2} = cellstr(fullfile(outpth,sprintf('GCM_%s.mat',GCMname))); % This only needed for plotting - doesn't restrict to two models
spm_jobman('run', jobfile, inputs{:});


%% 8. Run BMC on just two models (with/without between-region connections)
%==========================================================================
jobfile = {fullfile(scrpth,'batch_peb_bmc_job.m')}; 
inputs  = {};
inputs{1} = cellstr(fullfile(outpth,sprintf('PEB_%s_fit.mat',GCMname)));
inputs{2} = cellstr(fullfile(outpth,sprintf('GCM_%s_fit.mat',GCMname)));
spm_jobman('run', jobfile, inputs{:});


%% Check DCM reproduces activations

load(fullfile(outpth,'GCM_full_self_fit.mat'));
m = 1; % Cannot estimate model 2 because batch does not update GCM
cols = [1:2]; 
betas = []; dcm_betas = []; pvar = [];
for s = 1:nsub
    outdir = fullfile(outpth,subdir{s},'stat_concat')
    for r = 1:numel(ROI_names)        
        %% Get GLM betas 
        load(fullfile(outdir,'SPM.mat'));      
        VOI =  load(fullfile(outdir,sprintf('VOI_%s_1.mat',ROI_names{r})),'Y');
        tmp = pinv(SPM.xX.xKXs.X)*VOI.Y;
        betas(s,r,:) = tmp(cols);
        
        %% Get DCM betas
        load(GCM{s,m})
        Y2 = DCM.y(:,r); % fitted data
        Y2 = Y2/DCM.Y.scale;
        tmp = pinv(SPM.xX.xKXs.X)*Y2;
        dcm_betas(s,r,:) = tmp(cols);
        
        PSS   = sum(DCM.y(:,r).^2);
        RSS   = sum(DCM.R(:,r).^2);
        pvar(s,r) = PSS/(PSS + RSS);        
    end
    allF(s) = DCM.F;
end
figure,boxplot(pvar)
%figure,boxplot(allF)
mean(pvar)
%mean(allF)

cw = [0 1]; % Faces only
T=[]; p=[]; figure; clf
for r = 1:numel(ROI_names)   
    for c=1:size(cw,1)
         dc(:,c) = squeeze(betas(:,r,:))*cw(c,:)';
         T(c) = mean(dc(:,c))/(std(dc(:,c))/sqrt(size(dc,1)));
    end  
    p = t2p(T,size(dc,1)-1,2);
    
    fprintf('%s\n',ROI_names{r})
    disp([cw T' p'])
    
    for c=1:size(cw,1)
         dc(:,c) = squeeze(dcm_betas(:,r,:))*cw(c,:)';
         T(c) = mean(dc(:,c))/(std(dc(:,c))/sqrt(size(dc,1)));
    end  
    p = t2p(T,size(dc,1)-1,2);
    
    fprintf('DCM: %s\n',ROI_names{r})
    disp([cw T' p'])
    
    subplot(2,2,r),
    mbeta = mean(squeeze(betas(:,r,:)))';
    mbeta = mbeta/mean(mbeta)
    mdcm_beta = mean(squeeze(dcm_betas(:,r,:)))';  
    mdcm_beta = mdcm_beta/mean(mdcm_beta)
    bar([mbeta mdcm_beta])
    legend('real data','DCM fitted');
    Title = ROI_names;
    
    title(Title{r},'Interpreter','none')
end

b=2;
for r = 1:numel(ROI_names)   
    corr(betas(:,r,b),dcm_betas(:,r,b))
end
