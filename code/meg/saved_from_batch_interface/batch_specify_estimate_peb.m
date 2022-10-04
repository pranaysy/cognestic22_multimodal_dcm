% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'/imaging/henson/Wakeman/cognestic2022_dcm_meeg/code/saved_from_batch_interface/batch_specify_estimate_peb_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
