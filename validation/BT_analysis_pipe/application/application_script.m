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
    
    [t,sts] = spm_select(2,'image','Select files...',{},pwd,'.*','1');
    V = cellstr(t);
    spmup_auto_reorient(V,1);
    
    cd ..
end

%% flip images so the tumour is on the same side (left)

for patient = [3 6:8 14:18 21:22 26:27 29 36 41 45 50 53:54]
    cd(local(patient+2).name)
    
    folderinfo = dir(pwd);
    image_1 = folderinfo(3).name;
    
    if strfind (folderinfo(4).name,'T1w')
        image_2 =folderinfo(4).name;
    elseif strfind (folderinfo(5).name,'T1w')
        image_2 =folderinfo(5).name;
    end
    
    image_1_vol = spm_vol(image_1);
    image_1_read = spm_read_vols(image_1_vol);
    image_2_vol = spm_vol(image_2);
    image_2_read = spm_read_vols(image_2_vol);
    A = flipud(image_1_read);
    B = flipud(image_2_read);
    
    image_1 = spm_write_vol(image_1_vol,A);
    image_2 = spm_write_vol(image_2_vol,B);
    
    cd ..
end

%% USwL normalise/segment
 
for patient = 2:6
    cd(local(patient+2).name)
    
    folderinfo = dir(pwd);
    lesion_mask = folderinfo(3).name;
    image_1 = folderinfo(4).name;
    
    if strfind (folderinfo(5).name,'T1w')
        image_2 =folderinfo(5).name;
    elseif strfind (folderinfo(6).name,'T1w')
        image_2 =folderinfo(6).name;
    end
    
    for nbGaussian = 1:2
        for affectedtissue = 1:2 % add +1 for GM+WM or GM+WM+CSF
            out = USwL_segmentation(image_2,image_1,lesion_mask,nbGaussian,affectedtissue);
            % clean up
            destination = [pwd filesep 'nbG' num2str(nbGaussian) '_tissue' num2str(affectedtissue+1)];
            mkdir(destination)
            movefile(char(out{1}.segmImg.c1),[destination filesep 'c1.nii']);
            movefile(char(out{1}.segmImg.c2),[destination filesep 'c2.nii']);
            movefile(char(out{1}.segmImg.c3),[destination filesep 'c3.nii']);
            movefile(char(out{1}.segmImg.c4),[destination filesep 'c4.nii']);
            movefile(char(out{1}.segmImg.wc1),[destination filesep 'wc1.nii']);
            movefile(char(out{1}.segmImg.wc2),[destination filesep 'wc2.nii']);
            movefile(char(out{1}.segmImg.wc3),[destination filesep 'wc3.nii']);
            movefile(char(out{1}.segmImg.wc4),[destination filesep 'wc4.nii']);
            movefile(char(out{1}.segmImg.mwc1),[destination filesep 'mwc1.nii']);
            movefile(char(out{1}.segmImg.mwc2),[destination filesep 'mwc2.nii']);
            movefile(char(out{1}.segmImg.mwc3),[destination filesep 'mwc3.nii']);
            movefile(char(out{1}.segmImg.mwc4),[destination filesep 'mwc4.nii']);
            movefile(char(out{1}.segmImg.rc1),[destination filesep 'rc1.nii']);
            movefile(char(out{1}.segmImg.rc2),[destination filesep 'rc2.nii']);
            movefile(char(out{1}.segmImg.rc3),[destination filesep 'rc3.nii']);
            movefile(char(out{1}.ICVmsk),[destination filesep 'ICVmsk.nii']);
            movefile(char(out{1}.wICVmsk),[destination filesep 'wICVmsk.nii']);
            movefile(char(out{1}.Struc_1),[destination filesep 'Struc_1.nii']);
            movefile(char(out{1}.Struc_2),[destination filesep 'Struc_2.nii']);
            movefile(char(out{1}.wStruc_1),[destination filesep 'wStruc_1.nii']);
            movefile(char(out{1}.wStruc_2),[destination filesep 'wStruc_2.nii']);
            movefile(char(out{1}.TPMl),[destination filesep 'TPMl.nii']);
            movefile([pwd filesep 'wlesion mask.nii'],[destination filesep 'wlesion mask.nii']);
            movefile([pwd filesep 'swlesion mask.nii'],[destination filesep 'swlesion mask.nii']);
            movefile([pwd filesep 'dtlesion mask.nii'],[destination filesep 'dtlesion mask.nii']);
            movefile([pwd filesep 'BiasField_krsub*.nii'],[destination filesep 'BiasField_alt']);
            movefile([pwd filesep 'BiasField_ksub*.nii'],[destination filesep 'BiasField_T1']);
            movefile([pwd filesep 'iy*.nii'],[destination filesep 'iy']);
            movefile([pwd filesep 'y*.nii'],[destination filesep 'y']);
            movefile([pwd filesep 'krsub*.nii'],[destination filesep 'k_alt']);
            movefile([pwd filesep 'ksub*.nii'],[destination filesep 'k_T1']);
            movefile([pwd filesep 'mkrsub*.nii'],[destination filesep 'mk_alt']);
            movefile([pwd filesep 'mksub*.nii'],[destination filesep 'mk_T1']);
            movefile([pwd filesep 'swicv*.nii'],[destination filesep 'swicv']);
            movefile([pwd filesep 'wmkr*.nii'],[destination filesep 'wmk_alt']);
            movefile([pwd filesep 'wmk*.nii'],[destination filesep 'wmk_T1']);
            delete([pwd filesep 'c5*.nii']);
            delete([pwd filesep 'c6*.nii']);
            delete([pwd filesep 'ksub*.mat']);

        end
    end
    
    cd ..
end

