%-----------------------------------------------------------------------
% Job saved on 16-Aug-2022 17:19:05 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.stats.fmri_est.spmmat = '<UNDEFINED>';
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{2}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.con.consess{1}.fcon.name = 'Effects of Interest';
matlabbatch{2}.spm.stats.con.consess{1}.fcon.weights = [1 0; 0 1];
matlabbatch{2}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{2}.spm.stats.con.consess{2}.tcon.name = 'All > Baseline';
matlabbatch{2}.spm.stats.con.consess{2}.tcon.weights = [1 0];
matlabbatch{2}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{2}.spm.stats.con.delete = 0;
