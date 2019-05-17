function out = segment_with_lesion(inv_wtVSD,wT1w_with_tumour,nbGaussian,affectedtissue)

% do the segmentation
spm_jobman('initcfg')

%% using the batch

matlabbatch{1}.spm.tools.USwLtools.uswl.imgMsk = {inv_wtVSD};
matlabbatch{1}.spm.tools.USwLtools.uswl.imgRef = {wT1w_with_tumour};
matlabbatch{1}.spm.tools.USwLtools.uswl.imgStruc = '';
matlabbatch{1}.spm.tools.USwLtools.uswl.imgOth = '';
matlabbatch{1}.spm.tools.USwLtools.uswl.options.imgTpm = {[fileparts(which('spm')) filesep 'tpm' filesep 'TPM.nii']};
matlabbatch{1}.spm.tools.USwLtools.uswl.options.NbGaussian = [2 2 nbGaussian 2 2 0 1]; 
matlabbatch{1}.spm.tools.USwLtools.uswl.options.tpm4lesion = affectedtissue;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biasreg = 0.001;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biasfwhm = 60;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biaswr = [1 1];
matlabbatch{1}.spm.tools.USwLtools.uswl.options.ICVmsk = 1;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.mrf = 2;
out = spm_jobman('run', matlabbatch);


