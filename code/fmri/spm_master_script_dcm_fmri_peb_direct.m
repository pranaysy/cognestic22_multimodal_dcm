% SPM12 script to run DCM for fMRI on the multi-subject, multi-modal human neuroimaging dataset 
%
% Preprocessing described in Henson et al. (2019) doi.org/10.3389/fnins.2019.00300
%
% Note that you will need to have latest version of SPM12 on your MATLAB
% path, which you can download from here: 
%       https://www.fil.ion.ucl.ac.uk/spm/software/download/
%
% You can either:
% 
% 1. downlaod the raw data, available from the OpenNeuro database in BIDS format: 
%       https://openneuro.org/datasets/ds000117
% and run the preprocessing steps described in Appendix 2 of:
%       https://www.frontiersin.org/articles/10.3389/fnins.2019.00300/full#supplementary-material

% ...or 2. download preprocessed data from Figshare link provided in paper
%
% rik.henson@mrc-cbu.cam.ac.uk                              Aug 2022
% with help from Pranay Yadav

clear

SPM12PATH = '/imaging/local/software/spm_cbu_svn/releases/spm12_latest/' % insert path to your local install of SPM12
addpath(SPM12PATH);
spm fmri

%% Input arguments
%==========================================================================
rawpth = '/imaging/henson/users/rh01/Methods/DanData/OpenNeuro/ds000117' % where put raw BIDS data
derpth = fullfile(rawpth,'derivatives','SPM12') 
% above will exist if you have run preprocessing scripts for Henson et al, 2019, 
% or else download and extract fmri_data.tar.gz from Figshare link

outpth = fullfile(derpth,'DCM');
try mkdir(outpth); end; 

% Create sub-directory for scripts,
scrpth = fullfile(outpth,'code');                 
try mkdir(scrpth); end; 
% ... and download and extract code.tar.gz to here from Figshare link
addpath(scrpth) % Needed for parfor version of spm_dcm_peb_fit and spm_fmri_concatenate to return SPM

% If you have all raw data
% BIDS   = spm_BIDS(rawpth);
% subs   = spm_BIDS(BIDS,'subjects', 'task','facerecognition');
% runs = spm_BIDS(BIDS,'runs', 'modality','func', 'type','bold', 'task','facerecognition'); 

% ...else just re-specify
subs = {}; for s = [1:9 11:16]; subs{end+1} = sprintf('%02d',s); end % subject 10 had fewer scans in last run
nsub   = numel(subs)
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);

runs = {}; for r = [1:9]; runs{end+1} = sprintf('%02d',r); end
nrun = numel(runs)


nscan = repmat(208,1,nrun); % sub-15, run-01 has 209 scans, so ignore last
TR = 2;

numworkers = nsub; % Number of workers for distributed computing (depends on system)
if numworkers > 0
    delete(gcp('nocreate')) % Shut down any existing pool
%     P=cbupool(numworkers, '--mem-per-cpu=4G --time=12:00:00 --nodelist=node-j10');
%     parpool(P, P.NumWorkers);
    parpool(numworkers);
end


%% 1. Combine conditions and the concatenate runs for DCM
%==========================================================================

% Could just load separate-run SPM.mat files, concatenate and re-save, but
% recreate below in case not been estimated before

% Care: if run in parfor, will hang if SPM.mat files already exist, so make
% sure delete them before re-running

parfor (s = 1:nsub, numworkers)
    if numworkers > 0
        spm_jobman('initcfg');
        spm('defaults', 'fmri');
        spm_get_defaults('cmdline',true);
    end
        
    SPM = [];
    SPM.xY.RT  = 2;
    SPM.nscan = sum(nscan);
    SPM.xBF.T = 16;
    SPM.xBF.T0 = 8;
    SPM.xBF.name = 'hrf';
    SPM.xBF.UNITS = 'secs';
    SPM.xBF.Volterra = 1;
    
    conds = [];
    conds.names{1} = 'All';
    conds.names{2} = 'Faces';    % Famous and Nonfamous

    outdir = fullfile(outpth,subdir{s})
    try mkdir(outdir); end
    cd(outdir)
    for t = 1:2; conds.onsets{t} = []; end
      
    time_so_far = 0; volfiles = {}; movepar = [];
    for r = 1:length(runs)       
        volfiles{r} = spm_select('ExtFPList',fullfile(derpth,subdir{s},'func'),sprintf('^swsub-.*run-%s_bold\\.nii$',runs{r}),[1:nscan(r)]);
   
        trlfile = fullfile(derpth,subdir{s},'func',sprintf('sub-%s_run-%s_spmdef.mat',subs{s},runs{r}));
        d = load(trlfile);

        conds.onsets{1} = [conds.onsets{1}; sort([d.onsets{1}; d.onsets{2}; d.onsets{3}]) + time_so_far];
        conds.onsets{2} = [conds.onsets{2}; sort([d.onsets{1}; d.onsets{2}]) + time_so_far];
 
        time_so_far = time_so_far + nscan(r)*TR;
        
        d = load(spm_select('FPList',fullfile(derpth,subdir{s},'func'),sprintf('^rp.*run-%s.*\\.txt$',runs{r})));
        d = d(1:nscan(r),:);
        movepar = [movepar; d];
    end
    
    % Specify concatenated model
    %----------------------------------------------------------------------    
    
    for c = 1:length(conds.names)
        SPM.Sess(1).U(c).name{1} = conds.names{c};
        SPM.Sess(1).U(c).ons = conds.onsets{c};
        SPM.Sess(1).U(c).dur = 0;
        SPM.Sess(1).U(c).P.name = 'none';
    end
    
    SPM.Sess(1).C.C = movepar;
    SPM.Sess(1).C.name = {'x';'y';'z';'pitch';'roll';'yaw'};
     
    SPM.xY.P = strvcat(volfiles{:});
    SPM.xX.K.HParam = 128;  
    SPM.xVi.form = 'AR(0.2)';
    SPM.xGX.iGXcalc = 'none';
    
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

%% 2. Create VOI files (ie timeseries for ROIs)
%==========================================================================

ROI_names = {'lOFA','rOFA','lFFA','rFFA'};
ROI_coord = {[-38, -86, -14],[+36, -86, -10],[-42, -56, -20],[+42, -52, -14]};
rad = 10; % radius

% Preload SPMs just so parfor can work without loading
allSPMs = {};
for s = 1:nsub
    tmp = load(fullfile(outpth,subdir{s},'SPM.mat'));
    allSPMs{s} = tmp.SPM;
end
    
parfor (s = 1:nsub, numworkers) % can't parallelise with needing to load SPM
    if numworkers > 0
        spm_jobman('initcfg');
        spm('defaults', 'fmri');
        spm_get_defaults('cmdline',true);
    end
    
    outdir = fullfile(outpth,subdir{s})
    SPM = allSPMs{s};
%     spmfile = fullfile(outdir,'SPM.mat');
%     load(spmfile)
    
    [x,y,z] = ndgrid(1:SPM.xVol.DIM(1),1:SPM.xVol.DIM(2),1:SPM.xVol.DIM(3));
    XYZ     = [x(:),y(:),z(:)]'; 
    XYZmm   = SPM.xVol.M(1:3,:) * [XYZ;ones(1,size(XYZ,2))];
    Q       = ones(1,size(XYZmm,2));   
    
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
    
    for r = 1:length(ROI_names)
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

%% 3. Specify DCM for sub-01
%==========================================================================

ref_sub = 1; % eg first subject
outdir = fullfile(outpth,subdir{ref_sub});
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

DCM.b(:,:,2) = [
     1     0     0     0
     0     1     0     0
     0     0     1     0
     0     0     0     1
];
outfile = fullfile(outdir,'DCM_self.mat')
save(outfile,'DCM');


%% Estimate full model for single subject and review
load(fullfile(outdir,'DCM_full.mat'));
DCM = spm_dcm_fit(DCM,(numworkers > 0)); DCM = DCM{1};
save(fullfile(outdir,'DCM_full_estimated.mat'),'DCM');
spm_dcm_review(DCM);
spm_dcm_fmri_check(DCM);

%% 4. Specify GCM (DCM for each subject and model)
%==========================================================================
% (value of "ref_sub" must match value that generated DCMs above, eg 1)

models = {'full','self'};
altpar = struct([]); % no other alternative parameterisations for moment (cf Step 8 below)
altpar(2).B = [
     1     0     0     0
     0     1     0     0
     0     0     1     0
     0     0     0     1
];

GCM = {};
for s = 1:nsub
    outdir = fullfile(outpth,subdir{s})
    for m = 1:length(models)
        DCMfile = fullfile(outdir,['DCM_' models{m} '.mat']);
        
        if ~exist(DCMfile)
 
            load(fullfile(outpth,subdir{ref_sub},['DCM_' models{m} '.mat']));
            
            % Update VOI data
            for r = 1:numel(ROI_names)
                load(fullfile(outdir,sprintf('VOI_%s_1.mat',ROI_names{r})),'xY');
                DCM.xY(r) = xY;
            end 
           
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
            
            save(DCMfile,'DCM');
        end
        
        GCM{s,m} = DCMfile;
    end
end

cd(outpth)
save(['GCM_' models{1}],'GCM');    


%% 5. Estimate GCM 
%==========================================================================
GCM = spm_dcm_fit(GCM(:,1),(numworkers > 0));

% save(['GCM_' models{1}],'GCM'); % This saves big file with all estimated
% values, but to mimic batch interface, update the DCM files in the subject
% directories instead (note files not appended with "_m0001" etc like in
% batched version)

GCM_fit = GCM;
load(['GCM_' models{1}])
for s=1:size(GCM,1) 
    DCM = GCM_fit{s,1};
    save(GCM{s,1},'DCM')
end
clear GCM_fit
spm_dcm_fmri_check(GCM);

% Or iterative fit, but may take time
% M = []; M.X = ones(nsub,1);
% M.Q = 'single'; % default in batch for PEB, so assume also default in batch for peb_fit
% GCM = spm_dcm_peb_fit(GCM,M,{'A','B'},(numworkers > 0)); % 'A' and 'B' because batch only allows this
% save(['GCM_' models{1} 'peb_fit'],'GCM');


%% 6. Run PEB (on B matrices)
%==========================================================================
M = []; M.X = ones(nsub,1); M.Q = 'all'; 
[PEB, RCM] = spm_dcm_peb(GCM(:,1),M,{'B'});
save(['PEB_' models{1}]','PEB') %,'RCM');

spm_dcm_peb_review(PEB,GCM);


%% 7. Run BMC across all possible models (BMR)
%==========================================================================
[BMA,BMR] = spm_dcm_peb_bmc(PEB);
save(['BMA_search_PEB_' models{1}],'BMA','BMR')

spm_dcm_peb_review(BMA,GCM);


%% 8. Run BMC on just two models (with/without between-region connections)
%==========================================================================
[BMA, BMR] = spm_dcm_peb_bmc(PEB, GCM(1,:));
save(['BMA_PEB_' models{1}],'BMA','BMR')
spm_dcm_peb_review(BMA,GCM);




%% Check DCM reproduces activations
% 
% load(fullfile(outpth,'GCM_full_self_fit.mat'));
% m = 1; % Cannot estimate model 2 because batch does not update GCM
% cols = [1:2]; 
% betas = []; dcm_betas = []; pvar = [];
% for s = 1:nsub
%     outdir = fullfile(outpth,subdir{s},'stat_concat')
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
