%-----------------------------------------------------------------------
% Job saved on 25-Jan-2024 11:45:57 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%%
matlabbatch{1}.spm.dcm.spec.meeg.D = {
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-01/meg/wmaMcefspmeeg_sub-01_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-02/meg/wmaMcefspmeeg_sub-02_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-03/meg/wmaMcefspmeeg_sub-03_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-04/meg/wmaMcefspmeeg_sub-04_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-05/meg/wmaMcefspmeeg_sub-05_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-06/meg/wmaMcefspmeeg_sub-06_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-07/meg/wmaMcefspmeeg_sub-07_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-08/meg/wmaMcefspmeeg_sub-08_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-09/meg/wmaMcefspmeeg_sub-09_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-10/meg/wmaMcefspmeeg_sub-10_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-11/meg/wmaMcefspmeeg_sub-11_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-12/meg/wmaMcefspmeeg_sub-12_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-13/meg/wmaMcefspmeeg_sub-13_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-14/meg/wmaMcefspmeeg_sub-14_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-15/meg/wmaMcefspmeeg_sub-15_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data/derivatives/SPM12/sub-16/meg/wmaMcefspmeeg_sub-16_ses-meg_task-facerecognition_run-01_proc-sss_meg.mat'
                                      };
%%
matlabbatch{1}.spm.dcm.spec.meeg.dcmmat = {
                                           '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/fits/batch_gui/meg/templates/DCMs/DCM_Full.mat'
                                           '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/fits/batch_gui/meg/templates/DCMs/DCM_Self.mat'
                                           };
matlabbatch{1}.spm.dcm.spec.meeg.pE = {''};
matlabbatch{1}.spm.dcm.spec.meeg.P = {''};
matlabbatch{1}.spm.dcm.spec.meeg.feedback = 0;
matlabbatch{1}.spm.dcm.spec.meeg.output = 'GCM';
matlabbatch{1}.spm.dcm.spec.meeg.dir = {'/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/fits/batch_gui/meg/templates/GCMs/Full_vs_Self'};
matlabbatch{2}.spm.dcm.peb.compare.peb_mat = {'/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/fits/batch_gui/meg/PEB_Full.mat'};
matlabbatch{2}.spm.dcm.peb.compare.model_space_mat(1) = cfg_dep('DCM for M/EEG: GCM mat File(s)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','gcmmat'));
matlabbatch{2}.spm.dcm.peb.compare.show_review = 1;
