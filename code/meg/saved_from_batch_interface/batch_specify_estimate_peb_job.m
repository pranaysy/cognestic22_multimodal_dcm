%-----------------------------------------------------------------------
% Job saved on 24-Jan-2024 11:51:31 by cfg_util (rev $Rev: 7345 $)
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
matlabbatch{1}.spm.dcm.spec.meeg.dcmmat = {'/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/fits/batch_gui/meg/templates/DCMs/DCM_Full.mat'};
matlabbatch{1}.spm.dcm.spec.meeg.pE = {''};
matlabbatch{1}.spm.dcm.spec.meeg.P = {''};
matlabbatch{1}.spm.dcm.spec.meeg.feedback = 0;
matlabbatch{1}.spm.dcm.spec.meeg.output = 'GCM';
matlabbatch{1}.spm.dcm.spec.meeg.dir = {'/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/fits/batch_gui/meg/templates/GCMs/Full'};
matlabbatch{2}.spm.dcm.estimate.dcms.gcmmat(1) = cfg_dep('DCM for M/EEG: GCM mat File(s)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','gcmmat'));
matlabbatch{2}.spm.dcm.estimate.output.single.dir = {'/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/fits/batch_gui/meg'};
matlabbatch{2}.spm.dcm.estimate.output.single.name = 'Full';
matlabbatch{2}.spm.dcm.estimate.est_type = 1;
matlabbatch{2}.spm.dcm.estimate.fmri.analysis = 'time';
matlabbatch{3}.spm.dcm.peb.specify.name = 'Full';
matlabbatch{3}.spm.dcm.peb.specify.model_space_mat(1) = cfg_dep('DCM estimation: GCM mat File(s)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','gcmmat'));
matlabbatch{3}.spm.dcm.peb.specify.dcm.index = 1;
matlabbatch{3}.spm.dcm.peb.specify.cov.none = struct([]);
matlabbatch{3}.spm.dcm.peb.specify.fields.custom = {'B'};
matlabbatch{3}.spm.dcm.peb.specify.priors_between.components = 'All';
matlabbatch{3}.spm.dcm.peb.specify.priors_between.ratio = 16;
matlabbatch{3}.spm.dcm.peb.specify.priors_between.expectation = 0;
matlabbatch{3}.spm.dcm.peb.specify.priors_between.var = 0.0625;
matlabbatch{3}.spm.dcm.peb.specify.priors_glm.group_ratio = 1;
matlabbatch{3}.spm.dcm.peb.specify.estimation.maxit = 64;
matlabbatch{3}.spm.dcm.peb.specify.show_review = 1;
