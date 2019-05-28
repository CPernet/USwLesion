function standard_segment =standard_spm_segment_job(wT1w)

%-----------------------------------------------------------------------
% Job saved on 22-May-2019 11:16:21 by cfg_util (rev $Rev: 6942 $)
% spm SPM - SPM12 (7219)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
spm_jobman('initcfg')
matlabbatch{1}.spm.spatial.preproc.channel.vols = {wT1w};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'C:\Users\s1835343\mri_stuff\spm12\tpm\TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'C:\Users\s1835343\mri_stuff\spm12\tpm\TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'C:\Users\s1835343\mri_stuff\spm12\tpm\TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'C:\Users\s1835343\mri_stuff\spm12\tpm\TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 0;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'C:\Users\s1835343\mri_stuff\spm12\tpm\TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'C:\Users\s1835343\mri_stuff\spm12\tpm\TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [1 1];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];

standard_segment = spm_jobman('run', matlabbatch);
