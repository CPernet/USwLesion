 %% Application

%% set the directory and defaults

current = pwd;
global defaults
defaults = spm('defaults','FMRI');

% data dir
patient_dir = uigetdir(pwd,'Select VBM-BT directory');

%% co-register the Flair/T2 images for each patient to the T1

Healthy_dir = uigetdir(pwd,'Select F1000 directory');

for patient = 1:54
    
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {'C:\Users\s1835343\mri_stuff\VBM-BT\016bcf55-37c5-45e7-b95f-0a3ff0176354\sub-016bcf55-37c5-45e7-b95f-0a3ff0176354_T1w.nii,1'};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {'C:\Users\s1835343\mri_stuff\VBM-BT\016bcf55-37c5-45e7-b95f-0a3ff0176354\sub-016bcf55-37c5-45e7-b95f-0a3ff0176354_T2w.nii,1'};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    
end






