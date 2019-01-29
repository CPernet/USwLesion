function out = segmentation_routine(T1name,FLAIRname,mask,nbGaussian,affectedtissue)

spm_jobman('initcfg')

% do the segmentation
matlabbatch{1}.spm.tools.USwLtools.uswl.imgMsk = {mask};
matlabbatch{1}.spm.tools.USwLtools.uswl.imgRef = {T1name};
matlabbatch{1}.spm.tools.USwLtools.uswl.imgStruc = {FLAIRname};
matlabbatch{1}.spm.tools.USwLtools.uswl.imgOth = '';
matlabbatch{1}.spm.tools.USwLtools.uswl.options.imgTpm = {'C:\Users\s1835343\mri_stuff\spm12\toolbox\USwLesion\eTPM\eTPM_wBG.nii'}; % {[fileparts(which('spm')) filesep 'tpm' filesep 'TPM.nii']};
matlabbatch{1}.spm.tools.USwLtools.uswl.options.NbGaussian = [2 2 nbGaussian 2 2 0 2]; 
matlabbatch{1}.spm.tools.USwLtools.uswl.options.tpm4lesion = affectedtissue;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biasreg = 0.001;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biasfwhm = 60;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biaswr = [1 1];
matlabbatch{1}.spm.tools.USwLtools.uswl.options.ICVmsk = 1;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.mrf = 2;

% clean up
matlabbatch{2}.spm.tools.USwLtools.USwLutils.FxLesMsk.fnMsk(1) = cfg_dep('US with lesion: c3 image', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c3'));
matlabbatch{2}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.minVol = 10;
matlabbatch{2}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.fnOth = '';

out = spm_jobman('run', matlabbatch);