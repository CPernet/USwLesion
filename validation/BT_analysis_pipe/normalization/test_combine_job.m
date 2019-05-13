function combine = test_combine_job()


%-----------------------------------------------------------------------
% Job saved on 13-May-2019 14:40:52 by cfg_util (rev $Rev: 6942 $)
% spm SPM - SPM12 (7219)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

spm_jobman('initcfg')


matlabbatch{1}.spm.util.imcalc.input = {
                                        'C:\Users\s1835343\mri_stuff\BRAT\F1000\0025336\T1w_skull_stripped.nii,1'
                                        'C:\Users\s1835343\mri_stuff\BRAT\BRATS2015_Training\LGG\brats_2013_pat0001_1\VSD.Brain.XX.O.MR_T1.54633\normalization_segmentation_nbG2_tissue2\rc3kVSD_thresholded.nii,1'    
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'output';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = '(i1.*((i1>0)-(i2>0))) ';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 1;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

combine = spm_jobman('run', matlabbatch);
