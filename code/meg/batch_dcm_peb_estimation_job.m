%-----------------------------------------------------------------------
% Job saved on 04-Aug-2022 11:26:07 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.dcm.spec.meeg.D = '<UNDEFINED>';
matlabbatch{1}.spm.dcm.spec.meeg.dcmmat = '<UNDEFINED>';
matlabbatch{1}.spm.dcm.spec.meeg.pE = {''};
matlabbatch{1}.spm.dcm.spec.meeg.P = {''};
matlabbatch{1}.spm.dcm.spec.meeg.feedback = 0;
matlabbatch{1}.spm.dcm.spec.meeg.output = 'GCM';
matlabbatch{1}.spm.dcm.spec.meeg.dir = '<UNDEFINED>';
matlabbatch{2}.spm.dcm.estimate.dcms.gcmmat(1) = cfg_dep('DCM for M/EEG: GCM mat File(s)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','gcmmat'));
matlabbatch{2}.spm.dcm.estimate.output.single.dir = '<UNDEFINED>';
matlabbatch{2}.spm.dcm.estimate.output.single.name = '<UNDEFINED>';
matlabbatch{2}.spm.dcm.estimate.est_type = 1;
matlabbatch{2}.spm.dcm.estimate.fmri.analysis = 'time';
matlabbatch{3}.spm.dcm.peb.specify.name = '<UNDEFINED>';
matlabbatch{3}.spm.dcm.peb.specify.model_space_mat(1) = cfg_dep('DCM estimation: GCM mat File(s)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','gcmmat'));
matlabbatch{3}.spm.dcm.peb.specify.dcm.index = 1;
matlabbatch{3}.spm.dcm.peb.specify.cov.none = struct([]);
%matlabbatch{3}.spm.dcm.peb.specify.fields.default = {
%                                                     'A'
%                                                     'B'
%                                                     }';
matlabbatch{3}.spm.dcm.peb.specify.fields.custom = {'B'};
matlabbatch{3}.spm.dcm.peb.specify.priors_between.components = 'All';
matlabbatch{3}.spm.dcm.peb.specify.priors_between.ratio = 16;
matlabbatch{3}.spm.dcm.peb.specify.priors_between.expectation = 0;
matlabbatch{3}.spm.dcm.peb.specify.priors_between.var = 0.0625;
matlabbatch{3}.spm.dcm.peb.specify.priors_glm.group_ratio = 1;
matlabbatch{3}.spm.dcm.peb.specify.estimation.maxit = 64;
matlabbatch{3}.spm.dcm.peb.specify.show_review = 1;
