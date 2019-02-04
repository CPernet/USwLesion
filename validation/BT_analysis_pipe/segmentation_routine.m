function out = segmentation_routine(T1name,FLAIRname,mask,nbGaussian,affectedtissue)


% do the segmentation
spm_jobman('initcfg')

% Mask image
fn_in{1} = mask;
% Structural reference for affine transform
fn_in{2} = T1name;
% All Structurals used for segmentation
imgs{1} = T1name; 
imgs{2} = FLAIRname;
fn_in{3} = char(imgs');
% others
fn_in{4} = char('');

% option sttucture
options = struct('imgTpm',[],'NbGaussian',[2 2 nbGaussian 2 2 2], 'tpm4lesion', affectedtissue, ...
    'ICVmsk', 1, 'mrf', 2, 'biasreg', 1.0000e-03, 'biasfwhm', 60, 'biaswr', [1 1], 'cleanup', 0);

TPM = {[fileparts(which('spm')) filesep 'toolbox' filesep 'USwLesion' filesep 'validation' filesep 'BT_analysis_pipe' filesep 'TPM_noskull.nii']};
options.imgTpm = TPM;

% run it
out = crc_USwL(fn_in,options);


% same using the batch
% matlabbatch{1}.spm.tools.USwLtools.uswl.imgMsk = 
% matlabbatch{1}.spm.tools.USwLtools.uswl.imgRef = 
% matlabbatch{1}.spm.tools.USwLtools.uswl.imgStruc = imgs';
% matlabbatch{1}.spm.tools.USwLtools.uswl.imgOth = '';
% matlabbatch{1}.spm.tools.USwLtools.uswl.options.imgTpm =  {[fileparts(which('spm')) filesep 'tpm' filesep 'TPM.nii']};
% matlabbatch{1}.spm.tools.USwLtools.uswl.options.NbGaussian = [2 2 nbGaussian 2 2 0 2]; 
% matlabbatch{1}.spm.tools.USwLtools.uswl.options.tpm4lesion = affectedtissue;
% matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biasreg = 0.001;
% matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biasfwhm = 60;
% matlabbatch{1}.spm.tools.USwLtools.uswl.options.bias.bias_yes.biaswr = [1 1];
% matlabbatch{1}.spm.tools.USwLtools.uswl.options.ICVmsk = 1;
% matlabbatch{1}.spm.tools.USwLtools.uswl.options.mrf = 2;
% out = spm_jobman('run', matlabbatch);