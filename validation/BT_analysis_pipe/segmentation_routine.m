function out = segmentation_routine(T1name,FLAIRname,mask,nbGaussian,affectedtissue)

spm_jobman('initcfg')

% do the segmentation
matlabbatch{1}.spm.tools.USwLtools.uswl.imgMsk = {mask};
matlabbatch{1}.spm.tools.USwLtools.uswl.imgRef = {T1name};
imgs{1} = T1name; imgs{2} = FLAIRname;
matlabbatch{1}.spm.tools.USwLtools.uswl.imgStruc = imgs';
matlabbatch{1}.spm.tools.USwLtools.uswl.imgOth = '';
matlabbatch{1}.spm.tools.USwLtools.uswl.options.imgTpm =  {[fileparts(which('spm')) filesep 'tpm' filesep 'TPM.nii']}; % {'C:\Users\s1835343\mri_stuff\spm12\toolbox\USwLesion\eTPM\eTPM_wBG.nii'}; % {[fileparts(which('spm')) filesep 'tpm' filesep 'TPM.nii']};
matlabbatch{1}.spm.tools.USwLtools.uswl.options.NbGaussian = [2 2 nbGaussian 2 2 0 2]; 
matlabbatch{1}.spm.tools.USwLtools.uswl.options.tpm4lesion = affectedtissue;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biasreg = 0.001;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biasfwhm = 60;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biaswr = [1 1];
matlabbatch{1}.spm.tools.USwLtools.uswl.options.ICVmsk = 1;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.mrf = 2;

out = spm_jobman('run', matlabbatch);