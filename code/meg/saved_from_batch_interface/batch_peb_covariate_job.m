%-----------------------------------------------------------------------
% Job saved on 10-Nov-2022 18:57:13 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.dcm.peb.specify.name = 'Age';
matlabbatch{1}.spm.dcm.peb.specify.model_space_mat = {'/imaging/henson/Wakeman/cognestic22_multimodal_dcm/fits/batch_gui/meg/GCM_Full.mat'};
matlabbatch{1}.spm.dcm.peb.specify.dcm.index = 1;
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor.name = 'Age';
%%
matlabbatch{1}.spm.dcm.peb.specify.cov.regressor.value = [4.625
                                                          -1.375
                                                          3.625
                                                          -0.375
                                                          -3.375
                                                          -0.375
                                                          4.625
                                                          -0.375
                                                          2.625
                                                          -3.375
                                                          -2.375
                                                          -2.375
                                                          -1.375
                                                          -2.375
                                                          3.625
                                                          -1.375];
%%
matlabbatch{1}.spm.dcm.peb.specify.fields.custom = {'B'};
matlabbatch{1}.spm.dcm.peb.specify.priors_between.components = 'All';
matlabbatch{1}.spm.dcm.peb.specify.priors_between.ratio = 16;
matlabbatch{1}.spm.dcm.peb.specify.priors_between.expectation = 0;
matlabbatch{1}.spm.dcm.peb.specify.priors_between.var = 0.0625;
matlabbatch{1}.spm.dcm.peb.specify.priors_glm.group_ratio = 1;
matlabbatch{1}.spm.dcm.peb.specify.estimation.maxit = 64;
matlabbatch{1}.spm.dcm.peb.specify.show_review = 1;
