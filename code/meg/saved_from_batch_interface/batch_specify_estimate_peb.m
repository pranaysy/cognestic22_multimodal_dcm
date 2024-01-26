% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/code/meg/saved_from_batch_interface/batch_specify_estimate_peb_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
