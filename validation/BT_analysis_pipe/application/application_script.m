%% Application

%% set the directory and defaults

current = pwd;
global defaults
defaults = spm('defaults','FMRI');

% data dir
patient_dir = uigetdir(pwd,'Select VBM-BT directory');
cd(patient_dir); local = dir;

%% co-register the Flair/T2 images for each patient to the T1

for patient = 1:54
    cd(local(patient+2).name)
    
    folderinfo = dir(pwd);
    if strfind (folderinfo(3).name,'T1w')
        Ref =folderinfo(3).name;
    else
        Coreg =folderinfo(3).name; 
    end
    
    if strfind (folderinfo(4).name,'T1w')
        Ref =folderinfo(4).name;
    else
        Coreg =folderinfo(4).name;
    end
    
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {Ref};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {Coreg};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    spm_jobman('run', matlabbatch);
    
    cd ..
    
end

%% set origin to AC

for patient = 1:54
    cd(local(patient+2).name)
    
    folderinfo = dir(pwd);
    alt_image = folderinfo(3).name;
    
    if strfind (folderinfo(4).name,'T1w')
        ref_image =folderinfo(4).name;
    elseif strfind (folderinfo(5).name,'T1w')
        ref_image =folderinfo(5).name;
    end
    

%     spm_select(XXX)
%     spmup_auto_reorient(XXX,1);

end





