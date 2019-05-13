function subtract = test_subtraction_job(healthy_brain,c3_warped_mask)

%-----------------------------------------------------------------------
% Job saved on 13-May-2019 10:33:49 by cfg_util (rev $Rev: 6942 $)
% spm SPM - SPM12 (7219)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
spm_jobman('initcfg')

matlabbatch{1}.spm.util.imcalc.input = {
                                        [healthy_brain]
                                        [c3_warped_mask]
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'output';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'i1-i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

subtract = spm_jobman('run', matlabbatch);
