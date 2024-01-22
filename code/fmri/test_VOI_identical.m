% Sample script check if VOI masks are identical
sub = 13;

%% Load VOI files and show total diff
rootdir = '/imaging/henson/Wakeman/pranay_does_things/CBU_Neuroimaging_2024_Test_Reset/data';
subdirs = {'derivatives', 'SPM12', sprintf('sub-%02d', sub), 'fmri', 'CatGLM'};
fnames = {'VOI_bVC_mask.nii', 'VOI_lFFA_mask.nii', 'VOI_rFFA_mask.nii'};

% VOIs extracted using batch script are located under: 
%       fullfile(rootdir, 'derivatives', 'SPM12')
% VOIs extracted using direct script are located under: 
%       fullfile(rootdir, 'testderiv', 'derivatives', 'SPM12')

fbatch = fullfile(rootdir, subdirprefix{:}, fnames);
fdirect = fullfile(rootdir, 'testderiv', subdirprefix{:}, fnames);

nfiles = length(fnames);
vbatch = cell(1,nfiles); vdirect = cell(1,nfiles);
fprintf('Diff between batch and direct for sub-%02d\n', sub)
for k=1:nfiles
    vbatch{1,k} = spm_read_vols(spm_vol(fbatch{k}));
    vdirect{1,k} = spm_read_vols(spm_vol(fdirect{k}));
    
    % Simple sum of diff
    vdiff = vbatch{k}(:) - vdirect{k}(:);
    svdiff = sum(vdiff);
    fprintf('\t%s: %d\t(%d voxels differ)\n', fnames{k}, svdiff, nnz(vdiff))
end

