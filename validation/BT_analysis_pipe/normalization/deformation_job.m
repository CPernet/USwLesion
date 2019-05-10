function apply_deform = deformation_job(deform_field,mask_to_deform)

%-----------------------------------------------------------------------
% Job saved on 10-May-2019 14:19:08 by cfg_util (rev $Rev: 6942 $)
% spm SPM - SPM12 (7219)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

spm_jobman('initcfg')

matlabbatch{1}.spm.spatial.realign.write.data = {
                                                 [deform_field]
                                                 [mask_to_deform]
                                                 };
matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 1 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';

apply_deform = spm_jobman('run', matlabbatch);
