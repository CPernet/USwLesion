function test_combine = test_combine_job()


%-----------------------------------------------------------------------
% Job saved on 10-May-2019 15:48:19 by cfg_util (rev $Rev: 6942 $)
% spm SPM - SPM12 (7219)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
spm_jobman('initcfg')

matlabbatch{1}.spm.util.imcalc.input = {
                                        'C:\Users\s1835343\mri_stuff\BRAT\F1000\0025337\T1w_skull_stripped.nii,1'
                                        'C:\Users\s1835343\mri_stuff\BRAT\BRATS2015_Training\LGG\brats_2013_pat0002_1\VSD.Brain.XX.O.MR_T1.54639\normalization_segmentation_nbG3_tissue3\rc3kVSD.nii,1'
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'test_combine_images';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'i2.*(i1>100)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

test_combine = spm_jobman('run', matlabbatch);

