%-----------------------------------------------------------------------
% Job saved on 19-Sep-2022 11:08:38 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%%
matlabbatch{1}.spm.dcm.spec.meeg.D = {
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-01/meg/maMceffdspmeeg_sub-01_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-02/meg/maMceffdspmeeg_sub-02_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-03/meg/maMceffdspmeeg_sub-03_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-04/meg/maMceffdspmeeg_sub-04_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-05/meg/maMceffdspmeeg_sub-05_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-06/meg/maMceffdspmeeg_sub-06_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-07/meg/maMceffdspmeeg_sub-07_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-08/meg/maMceffdspmeeg_sub-08_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-09/meg/maMceffdspmeeg_sub-09_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-10/meg/maMceffdspmeeg_sub-10_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-11/meg/maMceffdspmeeg_sub-11_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-12/meg/maMceffdspmeeg_sub-12_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-13/meg/maMceffdspmeeg_sub-13_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-14/meg/maMceffdspmeeg_sub-14_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-15/meg/maMceffdspmeeg_sub-15_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/data/derivatives/SPM12/sub-16/meg/maMceffdspmeeg_sub-16_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      };
%%
matlabbatch{1}.spm.dcm.spec.meeg.dcmmat = {
                                           '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/templates/DCMs/DCM_Full.mat'
                                           '/imaging/henson/Wakeman/cognestic2022_dcm_meeg/templates/DCMs/DCM_Self.mat'
                                           };
matlabbatch{1}.spm.dcm.spec.meeg.pE = {''};
matlabbatch{1}.spm.dcm.spec.meeg.P = {''};
matlabbatch{1}.spm.dcm.spec.meeg.feedback = 0;
matlabbatch{1}.spm.dcm.spec.meeg.output = 'GCM';
matlabbatch{1}.spm.dcm.spec.meeg.dir = {'/imaging/henson/Wakeman/cognestic2022_dcm_meeg/templates/GCMs/Full_vs_Self'};
matlabbatch{2}.spm.dcm.peb.compare.peb_mat = {'/imaging/henson/Wakeman/cognestic2022_dcm_meeg/fits/PEB_Full.mat'};
matlabbatch{2}.spm.dcm.peb.compare.model_space_mat(1) = cfg_dep('DCM for M/EEG: GCM mat File(s)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','gcmmat'));
matlabbatch{2}.spm.dcm.peb.compare.show_review = 1;
