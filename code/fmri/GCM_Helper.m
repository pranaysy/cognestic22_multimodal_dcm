% Specify root working directory 
base_dir = '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset'; % Change this to yours

% Assign operational directories to variables
derpth = fullfile(base_dir, 'data','derivatives','SPM12');
fits_dir = fullfile(base_dir, 'fits');

% Set up variables
subs = compose('%02g', [1:16]); % subject 10 had fewer scans in last run
nsub   = numel(subs);
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);

%% Load template DCM and repeat over all subjects
%---------------------------------------------------------------------------------------
% STEP 1: Load template DCM and prepare output variables
%---------------------------------------------------------------------------------------

outfile_full = fullfile(derpth, 'sub-15', 'fmri', 'CatGLM', 'DCM_Full.mat');

% Name this GCM (will be used as name of GCM file, with 'GCM_' appended)
GCM_name = 'Full';

% Path to template folder: This is where specified GCM will be saved
templatedir = fullfile(fits_dir, 'templates', 'GCMs', GCM_name);

%---------------------------------------------------------------------------------------
% STEP 2: Specify each subject's DCM from template DCM
%---------------------------------------------------------------------------------------

% Load template model
model = load(outfile_full); % This will load the 'Full' model
models = {GCM_name};        % Can be multiple
DCM_Full = model.DCM;
ROI_names = {DCM_Full.xY.name};

GCM = cell(nsub, 1);
for s = 1:nsub
    
    % Path to this subject's folder under derivatives
    outdir = fullfile(derpth,subdir{s}, 'fmri', 'CatGLM');
    
    % Iterate over models (we only have one here)
    for m = 1:length(models)
        
        % Specify output file for DCM specification
        DCMfile = fullfile(outdir,['DCM_' models{m} '.mat']);
        
        % Make copy of the full model template
        DCM = DCM_Full;
        
        % Update VOI data
        for r = 1:numel(ROI_names)
            load(fullfile(outdir,sprintf('VOI_%s_1.mat',ROI_names{r})), 'xY');
            DCM.xY(r) = xY;
        end 
        
        % Set up DCM.Y
        DCM.v   = length(DCM.xY(1).u); % number of time points
        DCM.Y.Q = spm_Ce(ones(1,DCM.n)*DCM.v);
        DCM.Y.X0  = DCM.xY(1).X0;
        DCM.Y.y = []; % clear to allow different nscans
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end

        % Update onsets    
        load(fullfile(outdir,'SPM.mat'));
        DCM.U.u = []; % clear to allow different nscans
        for u = 1:2
            DCM.U.u(:,u)  = SPM.Sess(1).U(u).u((32+1):end); % DCM allows for 2 TRs before first stimulus
            DCM.U.name{u} = SPM.Sess(1).U(u).name{1};
        end
        
        % Save this DCM specified for this subject
        save(DCMfile,'DCM');
        
        % Append to GCM array
        GCM{s,m} = DCMfile;
        
    end
end

% Save GCM
GCM_file = fullfile(templatedir, 'GCM_Full.mat');
save(GCM_file,'GCM');  

%% Make GCM with modelspace for full vs self-only models
%  Load GCM with 'Full' model
load(GCM_file)
GCM = spm_dcm_load(GCM);

% Remove priors (legacy, needed to prevent errors downstream)
GCM = cellfun(@(x) rmfield(x, 'M'), GCM, 'UniformOutput', false);

% Create a second column for 'Self' model
GCM(:, 2) = GCM;

% Loop over each row (subject) and set 'b' matrix to diagonal
for k=1:length(GCM)
    GCM{k,2}.b(:,:,2) = eye(3);
end

% Save GCM 
templatedir = fullfile(fits_dir, 'templates', 'GCMs', 'Full_vs_Self');
save(fullfile(templatedir, 'GCM_Full_vs_Self.mat'), 'GCM')
