%-----------------------------------------------------------------------
% Job saved on 24-Jan-2019 12:37:29 by cfg_util (rev $Rev: 6942 $)
% spm SPM - SPM12 (7219)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.USwLtools.uswl.imgMsk = {'C:\Users\s1835343\mri_stuff\BRAT\BRATS2015_Training\HGG\brats_2013_pat0001_1\VSD.Brain.XX.O.MR_Flair.54512\voi1.nii,1'};
matlabbatch{1}.spm.tools.USwLtools.uswl.imgRef = {'C:\Users\s1835343\mri_stuff\BRAT\BRATS2015_Training\HGG\brats_2013_pat0001_1\VSD.Brain.XX.O.MR_T1.54513\VSD.nii,1'};
matlabbatch{1}.spm.tools.USwLtools.uswl.imgStruc = {'C:\Users\s1835343\mri_stuff\BRAT\BRATS2015_Training\HGG\brats_2013_pat0001_1\VSD.Brain.XX.O.MR_Flair.54512\VSD.nii,1'};
matlabbatch{1}.spm.tools.USwLtools.uswl.imgOth = '';
matlabbatch{1}.spm.tools.USwLtools.uswl.options.imgTpm = {'C:\Users\s1835343\mri_stuff\spm12\toolbox\USwLesion\eTPM\eTPM_wBG.nii'};
matlabbatch{1}.spm.tools.USwLtools.uswl.options.NbGaussian = [2 2 3 2 3 4 2 1];
matlabbatch{1}.spm.tools.USwLtools.uswl.options.tpm4lesion = 3;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biasreg = 0.001;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biasfwhm = 60;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biaswr = [1 1];
matlabbatch{1}.spm.tools.USwLtools.uswl.options.ICVmsk = 1;
matlabbatch{1}.spm.tools.USwLtools.uswl.options.mrf = 2;
matlabbatch{2}.spm.tools.USwLtools.USwLutils.FxLesMsk.fnMsk(1) = cfg_dep('US with lesion: c3 image', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c3'));
matlabbatch{2}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.minVol = Inf;
matlabbatch{2}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.fnOth = '';
