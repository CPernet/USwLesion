%% Normalisation

%% set the directory and defaults

current = pwd;
global defaults
defaults = spm('defaults','FMRI');

% data dir
BRAT_dir = uigetdir(pwd,'Select BRAT directory');
if BRAT_dir == 0
    return
else
    cd(BRAT_dir);
    try
        cd('HGG'); BRAT_dir1 = pwd;
        cd ..
        cd('LGG'); BRAT_dir2 = pwd;
    catch
        error('HGG directory not found');
    end
end

Healthy_dir = uigetdir(pwd,'Select F1000 directory');

%% SPM segment healthy brains into standard space

cd(Healthy_dir); local = dir;
for subject = 1:30
    cd(local(subject+2).name)
    normalise = normalisation_batch_job();
    
    %inverse normalise skull-stripped brain back into subject space
    
    matlabbatch{1}.spm.util.bbox.image = {[pwd filesep 'anat.nii']};
    matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
    out = spm_jobman('run', matlabbatch);
    
    clear matlabbatch
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[pwd filesep 'iy_anat.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[pwd filesep 'T1w_skull_stripped.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = out{1}.bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    spm_jobman('run', matlabbatch);
    
    clear matlabbatch
    
    cd ..
end


%% intensity normalisation of healthy and patient brains (1st healthy brain matched to 1st patient brain, 2nd to 2nd, etc)

patient_to_segment = cell(30,1);
control_to_use     = cell(30,1);

for subject = 1:30
    cd(Healthy_dir); local = dir;
    cd(local(subject+2).name)
    HV = spm_vol([pwd filesep 'T1w_skull_stripped.nii']);
    healthy_V = spm_read_vols(HV);
    % start values at 0
    healthy_min = min(healthy_V(:));
    healthy_V   = healthy_V-healthy_min;
    healthy_max = max(healthy_V(:));
    
    if subject < 11
        tumour_type = 2;
    else
        tumour_type = 1;
    end
        % loop in High Grade Glioma or Low Grade Glioma
        % --------------------------------------------
        BRAT_dir = eval(['BRAT_dir' num2str(tumour_type)]);
        cd(BRAT_dir); folders = dir('brats*');
        if tumour_type == 1
            folders = folders(1:20);
            patient_dir = [BRAT_dir filesep folders(subject-10).name];
        else
            folders = folders(1:10);
            patient_dir = [BRAT_dir filesep folders(subject).name];
        end
                
        tmp = dir([patient_dir filesep 'VSD.Brain.XX.O.MR_T1.*']);
        patient_VSD = [patient_dir filesep tmp.name filesep 'VSD.nii'];
        PVSD = spm_vol(patient_VSD);
        patient_VSD = spm_read_vols(PVSD);
        patient_min = min(patient_VSD(:));
        patient_VSD = patient_VSD - patient_min;
        patient_max = max(patient_VSD(:));
        
        if patient_max < healthy_max
            scaling     = healthy_max/patient_max;
            healthy_V   = healthy_V./scaling;
        elseif patient_max > healthy_max
            scaling     = patient_max/healthy_max;
            patient_VSD = patient_VSD./scaling;
        end
        HV.fname    = [HV.fname(1:end-22) 'Scaled_T1w_skull_stripped.nii'];
        spm_write_vol(HV,healthy_V);
        control_to_use{subject}= HV.fname;
        PVSD.fname  = [PVSD.fname(1:end-7) 'Scaled_VSD.nii'];
        spm_write_vol(PVSD,patient_VSD);
        patient_to_segment{subject}= PVSD.fname;
end


%% tumour segmentation to get lesion masks

% segmentation
spm_jobman('initcfg')
 for subject = 1:30
    
    if subject < 11
        tumour_type = 2;
    else
        tumour_type = 1;
    end
    
    % loop in High Grade Glioma or Low Grade Glioma
    % --------------------------------------------
    BRAT_dir = eval(['BRAT_dir' num2str(tumour_type)]);
    cd(BRAT_dir); folders = dir('brats*');
    if tumour_type == 1
        folders = folders(1:20);
        patient_dir = [BRAT_dir filesep folders(subject-10).name];
        T1name = patient_to_segment{subject};
    else
        folders = folders(1:10);
        patient_dir = [BRAT_dir filesep folders(subject).name];
        T1name = patient_to_segment{subject};
    end

    tmp = dir([patient_dir filesep '*Flair*']);
    FLAIRname = [patient_dir filesep tmp.name filesep 'VSD.nii'];
    tmp = dir([patient_dir filesep 'VSD.Brain_*more*']);
    mask = [patient_dir filesep tmp.name filesep 'VSD.nii'];
    out = segmentation_routine(T1name,FLAIRname,mask,3,3);
    
    % normalized scaled_VSD
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = out;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[patient_to_segment{subject}]}; 
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72
        90 90 108];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1.5 1.5 1.5];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    MNI_patient = spm_jobman('run', matlabbatch);
    
    clear matlabbatch

        % normalized original tumour 
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = out;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {mask};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72
        90 90 108];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1.5 1.5 1.5];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    MNI_tumour = spm_jobman('run', matlabbatch);
    
    clear matlabbatch
    
end

%% set a tumour inside healthy brain
 
for subject = 1:30
    
    cd(Healthy_dir); local = dir;
    cd(local(subject+2).name)   
    T1w = [pwd filesep 'Scaled_T1w_skull_stripped.nii'];
    
    if subject < 11
        tumour_type = 2;
    else
        tumour_type = 1;
    end
    
    % loop in High Grade Glioma or Low Grade Glioma
    % --------------------------------------------
    BRAT_dir = eval(['BRAT_dir' num2str(tumour_type)]);
    cd(BRAT_dir); folders = dir('brats*');
    if tumour_type == 1
        folders = folders(1:20);
        patient_dir = [BRAT_dir filesep folders(subject-10).name];
    else
        folders = folders(1:10);
        patient_dir = [BRAT_dir filesep folders(subject).name];
    end
 
    tmp = dir([patient_dir filesep 'VSD.Brain_*more*']);
    Tumour = [patient_dir filesep tmp.name filesep 'wtVSD.nii'];

    tmp = dir([patient_dir filesep 'VSD.Brain.XX.O.MR_T1.*']);
    Patient = [patient_dir filesep tmp.name filesep 'wScaled_VSD.nii'];
    
    V1 = spm_vol(T1w);    T1w = spm_read_vols(V1);
    V2 = spm_vol(Tumour); Tumour_Mask = spm_read_vols(V2);
    V3 = spm_vol(Patient);Patient = spm_read_vols(V3);
    mask = Tumour_Mask > 0.1;
    T1w(mask) = Patient(mask);
    cd(Healthy_dir); local = dir;
    cd(local(subject+2).name)  
    V1.fname = [pwd filesep 'T1w_with_tumour.nii'];
    spm_write_vol(V1,T1w);
    
    %get the correct bounding box
    matlabbatch{1}.spm.util.bbox.image = {[pwd filesep 'anat.nii']};
    matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
    out = spm_jobman('run', matlabbatch);
    clear matlabbatch

    %inverse normalise the healthy brain+tumour AND tumour mask back into subject space
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[pwd filesep 'iy_anat.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[pwd filesep 'T1w_with_tumour.nii']
                                                                [Tumour]
                                                                };
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = out{1}.bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    inv = spm_jobman('run', matlabbatch);
    clear matlabbatch
    movefile(cell2mat(inv{1}.files(2)),[pwd filesep 'inv_wtVSD.nii']);
    
    %get the correct bounding box
    matlabbatch{1}.spm.util.bbox.image = {[pwd filesep 'wT1w_with_tumour.nii']};
    matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
    out = spm_jobman('run', matlabbatch);
    clear matlabbatch
    
    %inverse normalise the brain mask into healthy patient subject space
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[pwd filesep 'iy_anat.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[pwd filesep 'brain_mask.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = out{1}.bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    spm_jobman('run', matlabbatch);
    clear matlabbatch
       
end

%% segment and compute similarity 

% standard segmentation/USwL segmentation 

for subject = 1:30
    cd(Healthy_dir); local = dir;
    cd(local(subject+2).name)
    
    wT1w_skull_stripped = [pwd filesep 'wT1w_skull_stripped.nii'];
    wT1w_with_tumour = [pwd filesep 'wT1w_with_tumour.nii'];
    lesion_mask = [pwd filesep 'inv_wtVSD.nii']; lesion_mask_V = spm_read_vols(spm_vol(lesion_mask));
    
    
    %standard SPM segmentation of wT1w_skull_stripped.nii from healthy subject
    standard_H =standard_spm_segment_job(wT1w_skull_stripped);
    destination = [pwd filesep 'healthy_standard_segmentation'];
    mkdir(destination)
    movefile(char(standard_H{1}.tiss(1).c),[destination filesep 'c1wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(2).c),[destination filesep 'c2wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(3).c),[destination filesep 'c3wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(4).c),[destination filesep 'c4wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(5).c),[destination filesep 'c5wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(1).wc),[destination filesep 'wc1wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(2).wc),[destination filesep 'wc2wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(3).wc),[destination filesep 'wc3wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(5).wc),[destination filesep 'wc5wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(6).wc),[destination filesep 'wc6wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(1).mwc),[destination filesep 'mwc1wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(2).mwc),[destination filesep 'mwc2wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(3).mwc),[destination filesep 'mwc3wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(5).mwc),[destination filesep 'mwc5wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.tiss(6).mwc),[destination filesep 'mwc6wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.fordef),[destination filesep 'y_wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.invdef),[destination filesep 'iy_wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.param),[destination filesep 'wT1w_skull_stripped_seg8.mat']);
    movefile(char(standard_H{1}.channel.biasfield),[destination filesep 'BiasField_wT1w_skull_stripped.nii']);
    movefile(char(standard_H{1}.channel.biascorr),[destination filesep 'mwT1w_skull_stripped.nii']);
    %get the correct bounding box
    clear matlabbatch
    matlabbatch{1}.spm.util.bbox.image = {[pwd filesep 'healthy_tumour_USwL' filesep 'wmkwT1w_with_tumour.nii']};
    matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
    out = spm_jobman('run', matlabbatch);
    % apply y deformation field to mwT1w_with_tumour
    clear matlabbatch
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[destination filesep 'y_wT1w_skull_stripped.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[destination filesep 'mwT1w_skull_stripped.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = out{1}.bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    spm_jobman('run', matlabbatch);
    clear matlabbatch   
    
    %standard SPM segmentation of wT1w_with_tumour.nii
    standard_HT =standard_spm_segment_job(wT1w_with_tumour);
    destination = [pwd filesep 'healthy_tumour_standard_segmentation'];
    mkdir(destination)
    movefile(char(standard_HT{1}.tiss(1).c),[destination filesep 'c1wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(2).c),[destination filesep 'c2wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(3).c),[destination filesep 'c3wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(4).c),[destination filesep 'c4wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(5).c),[destination filesep 'c5wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(1).wc),[destination filesep 'wc1wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(2).wc),[destination filesep 'wc2wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(3).wc),[destination filesep 'wc3wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(5).wc),[destination filesep 'wc5wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(6).wc),[destination filesep 'wc6wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(1).mwc),[destination filesep 'mwc1wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(2).mwc),[destination filesep 'mwc2wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(3).mwc),[destination filesep 'mwc3wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(5).mwc),[destination filesep 'mwc5wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.tiss(6).mwc),[destination filesep 'mwc6wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.fordef),[destination filesep 'y_wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.invdef),[destination filesep 'iy_wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.param),[destination filesep 'wT1w_with_tumour_seg8.mat']);
    movefile(char(standard_HT{1}.channel.biasfield),[destination filesep 'BiasField_wT1w_with_tumour.nii']);
    movefile(char(standard_HT{1}.channel.biascorr),[destination filesep 'mwT1w_with_tumour.nii']);
    %get the correct bounding box
    clear matlabbatch
    matlabbatch{1}.spm.util.bbox.image = {[pwd filesep 'healthy_tumour_USwL' filesep 'wmkwT1w_with_tumour.nii']};
    matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
    out = spm_jobman('run', matlabbatch);
    % apply y deformation field to mwT1w_with_tumour
    clear matlabbatch
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[destination filesep 'y_wT1w_with_tumour.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[destination filesep 'mwT1w_with_tumour.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = out{1}.bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    spm_jobman('run', matlabbatch);
    clear matlabbatch    
    
    %USwL segmentation of wT1w_with_tumour.nii
    for nbGaussian = 1:2
        for affectedtissue = 1:2 % add +1 for GM+WM or GM+WM+CSF 
            seg_with_lesion_HT = segment_with_lesion(lesion_mask,wT1w_with_tumour,nbGaussian,affectedtissue);
            destination = [pwd filesep 'healthy_tumour_USwL'];
            %cleanup
            mkdir(destination)
            movefile(char(seg_with_lesion_HT{1}.segmImg.c1),[destination filesep 'c1kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.c2),[destination filesep 'c2kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.c3),[destination filesep 'c3kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.c4),[destination filesep 'c4kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.wc1),[destination filesep 'wc1kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.wc2),[destination filesep 'wc2kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.wc3),[destination filesep 'wc3kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.wc4),[destination filesep 'wc4kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.mwc1),[destination filesep 'mwc1kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.mwc2),[destination filesep 'mwc2kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.mwc3),[destination filesep 'mwc3kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.mwc4),[destination filesep 'mwc4kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.rc1),[destination filesep 'rc1kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.rc2),[destination filesep 'rc2kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.segmImg.rc3),[destination filesep 'rc3kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.ICVmsk),[destination filesep 'icv_kwT1w_with_tumour.nii.nii']);
            movefile(char(seg_with_lesion_HT{1}.wICVmsk),[destination filesep 'wicv_kwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.Struc_1),[destination filesep 'kmkwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.wStruc_1),[destination filesep 'kwmkwT1w_with_tumour.nii']);
            movefile(char(seg_with_lesion_HT{1}.TPMl),[destination filesep 'TPM_les.nii']);
        end
    end

cd ..

end

% test similarity 

SSIM = NaN(30,4);
rms = NaN(30,4);

for subject = 1:30
    
    cd(Healthy_dir); local = dir;
    cd(local(subject+2).name)
    Tumour = [pwd filesep 'inv_wtVSD.nii'];
    
    %generate brain mask without tumour
    lesion_mask = [Tumour]; lesion_mask_V = spm_read_vols(spm_vol(lesion_mask));
    brain_mask = [pwd filesep 'wbrain_mask.nii']; %brain mask from healthy subject in subject space
    X = spm_vol(brain_mask); brain_mask_V = spm_read_vols(X);
    Y = spm_vol(lesion_mask); lesion_mask = logical(spm_read_vols(Y));
    brain_mask_V(lesion_mask) = 0;
    X.fname = [pwd filesep 'brain_mask_minus_tumour.nii'];
    spm_write_vol(X,brain_mask_V);
    brain_mask_minus_tumour = [pwd filesep 'brain_mask_minus_tumour.nii'];
    
    %get the correct bounding box
    matlabbatch{1}.spm.util.bbox.image = {[pwd filesep 'healthy_tumour_USwL' filesep 'wmkwT1w_with_tumour.nii']};
    matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
    out = spm_jobman('run', matlabbatch);
    clear matlabbatch
    
    %normalise brain mask (minus lesion mask) and lesion mask so they are in same dimensions as other images
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[pwd filesep 'y_anat.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[pwd filesep 'brain_mask_minus_tumour.nii']
                                                                [Tumour]
                                                               };
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = out{1}.bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    out = spm_jobman('run', matlabbatch);
    SSIM_brain_mask_minus_tumour = out{1}.files{1};
    SSIM_lesion_mask = out{1}.files{2};
    SSIM_brain_mask_minus_tumour = spm_read_vols(spm_vol(SSIM_brain_mask_minus_tumour));
    SSIM_lesion_mask = spm_read_vols(spm_vol(SSIM_lesion_mask));
    clear matlabbatch out
   
    %get the images for measuring similarity of whole brain
    standard_healthy_brain = [pwd filesep 'healthy_standard_segmentation' filesep 'wmwT1w_skull_stripped.nii'];
    standard_healthy_brain = spm_read_vols(spm_vol(standard_healthy_brain));
    standard_healthy_tumour_brain = [pwd filesep 'healthy_tumour_standard_segmentation' filesep 'wmwT1w_with_tumour.nii'];
    standard_healthy_tumour_brain = spm_read_vols(spm_vol(standard_healthy_tumour_brain));
    USwL_healthy_tumour_brain = [pwd filesep 'healthy_tumour_USwL' filesep 'wmkwT1w_with_tumour.nii'];
    USwL_healthy_tumour_brain = spm_read_vols(spm_vol(USwL_healthy_tumour_brain));
    
    %measure similarity of whole brain
    SSIM(subject,1) = SSI(standard_healthy_brain,standard_healthy_tumour_brain,SSIM_lesion_mask,1);
    SSIM(subject,2) = SSI(standard_healthy_brain,standard_healthy_tumour_brain,SSIM_brain_mask_minus_tumour,1);
    SSIM(subject,3) = SSI(standard_healthy_brain,USwL_healthy_tumour_brain,SSIM_lesion_mask,1);
    SSIM(subject,4) = SSI(standard_healthy_brain,USwL_healthy_tumour_brain,SSIM_brain_mask_minus_tumour,1);

    % get dimensions of rms mwc1 (has different dimensions)
    matlabbatch{1}.spm.util.bbox.image = {[pwd filesep 'healthy_tumour_standard_segmentation' filesep 'mwc1wT1w_with_tumour.nii']};
    matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
    out = spm_jobman('run', matlabbatch);
    clear matlabbatch 
    
    %normalise brain mask (minus lesion mask) and lesion mask so they are in same dimensions as rms mwc1
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[pwd filesep 'y_anat.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[pwd filesep 'brain_mask_minus_tumour.nii']
                                                                [Tumour]
                                                               };
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = out{1}.bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1.5 1.5 1.5];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'rms';
    spm_jobman('run', matlabbatch);
    rms_brain_mask_minus_tumour = [pwd filesep 'rmsbrain_mask_minus_tumour.nii']; 
    rms_lesion_mask = [pwd filesep 'rmsinv_wtVSD.nii'];
    rms_brain_mask_minus_tumour = spm_read_vols(spm_vol(rms_brain_mask_minus_tumour));
    rms_lesion_mask = spm_read_vols(spm_vol(rms_lesion_mask));
    clear matlabbatch 
   
    %get the images for measuring rms
    c1_image_1 = [pwd filesep 'healthy_standard_segmentation' filesep 'mwc1wT1w_skull_stripped.nii'];
    c1_image_1V = spm_vol(c1_image_1); standard_healthy_c1 = spm_read_vols(c1_image_1V);
    c1_image_2 = [pwd filesep 'healthy_tumour_standard_segmentation' filesep 'mwc1wT1w_with_tumour.nii'];
    c1_image_2V = spm_vol(c1_image_2); standard_healthy_tumour_c1 = spm_read_vols(c1_image_2V);
    c1_image_3 = [pwd filesep 'healthy_tumour_USwL' filesep 'mwc1kwT1w_with_tumour.nii'];
    c1_image_3V = spm_vol(c1_image_3); USwL_healthy_tumour_c1 = spm_read_vols(c1_image_3V);
    
    %do we get the same values per voxel?
    img1 = standard_healthy_c1(logical(rms_lesion_mask)); 
    img2 = standard_healthy_tumour_c1(logical(rms_lesion_mask));
    rms(subject,1) = sqrt(mean((img1(:)-img2(:)).^2));
    img1 = standard_healthy_c1(logical(logical(rms_brain_mask_minus_tumour))); 
    img2 = standard_healthy_tumour_c1(logical(rms_brain_mask_minus_tumour));
    rms(subject,2) = sqrt(mean((img1(:)-img2(:)).^2));
    img1 = standard_healthy_c1(logical(rms_lesion_mask)); 
    img2 = USwL_healthy_tumour_c1(logical(rms_lesion_mask));
    rms(subject,3) = sqrt(mean((img1(:)-img2(:)).^2));
    img1 = standard_healthy_c1(logical(logical(rms_brain_mask_minus_tumour))); 
    img2 = USwL_healthy_tumour_c1(logical(rms_brain_mask_minus_tumour));
    rms(subject,4) = sqrt(mean((img1(:)-img2(:)).^2));

    cd ..
end

%write results as table
cd(fileparts(which('normalisation_script.m')))
title = {'lesion_standardhealthy_VS_standardtumour';'brainmask_standardhealthy_VS_standardtumour';'lesion_standardhealthy_VS_USwL';'brainmask_standardhealthy_VS_USwL'};
SSIM_results = array2table(SSIM,'VariableNames',title);
writetable(SSIM_results,[pwd filesep 'SSIM_results.csv']);
rms_results = array2table(rms,'VariableNames',title);
writetable(rms_results,[pwd filesep 'rms_results.csv']);

IMP = importdata('SSIM_results.csv'); SSIM = IMP.data;
IMP = importdata('rms_results.csv'); rms = IMP.data;

% test whether there is a statistical difference in similarity
%SSIM

SSIM_lesion = SSIM(:,[1 3]);
SSIM_brain = SSIM(:,[2 4]);

s = 1:30; s(19) = [];
SSIM_lesion_norm = SSIM_lesion(s,:);
%plot for SSIM within tumour (non-outlier results)
figure; 
for sub = 1:29
    t = [1:29];
    x = repmat(t,29);
    scatter(x(:,1),SSIM(s,1),30,'filled');
    hold on;
    scatter(x(:,2),SSIM(s,3),30,'filled');
    hold on;
    plot([1 2],SSIM_lesion_norm(sub,[1 2]),'LineWidth',2);
    hold on;
end
box on; grid on
axis([0.9 2.1 -1*10^5 7*10^5])
title('SSIM of normalized tumour area compared to healthy area','FontSize',16)
ylabel('SSIM','FontSize',14)

o = 19; SSIM_lesion_outlier = SSIM_lesion(o,:);
%plot for SSIM within tumour (outlier result)
figure;
scatter(1,SSIM(o,1),30,'filled');
hold on;
scatter(2,SSIM(o,3),30,'filled');
hold on;
plot([1 2],SSIM_lesion_outlier(1,[1 2]),'LineWidth',2);
hold on;
box on; grid on
axis([0.9 2.1 1*10^10 2*10^10])
title('SSIM of normalized tumour area compared to healthy area','FontSize',16)
ylabel('SSIM','FontSize',14)

s = 1:30; s([2 19 27 28 29]) = []; SSIM_brain_norm = SSIM_brain(s,:);
%plot for SSIM of brain minus tumour
figure;
for sub = 1:25
    t = [1:25];
    x = repmat(t,25);
    scatter(x(:,1),SSIM(s,2),30,'filled');
    hold on;
    scatter(x(:,2),SSIM(s,4),30,'filled');
    hold on;
    plot([1 2],SSIM_brain_norm(sub,[1 2]),'LineWidth',2);
    hold on;
end
box on; grid on
axis([0.9 2.1 0*10^5 18*10^6])
title('SSIM of normalized brain (minus tumour) compared to healthy brain (minus tumour)','FontSize',16)
ylabel('SSIM','FontSize',14)

o = [2 19 27 28 29]; SSIM_brain_outlier = SSIM_brain(o,:);
%plot for SSIM of brain minus tumour (outliers)
figure;
for sub = 1:5
    t = [1:5];
    x = repmat(t,5);
    scatter(x(:,1),SSIM(o,2),30,'filled');
    hold on;
    scatter(x(:,2),SSIM(o,4),30,'filled');
    hold on;
    plot([1 2],SSIM_brain_outlier(sub,[1 2]),'LineWidth',2);
    hold on;
end
box on; grid on
axis([0.9 2.1 3*10^9 5*10^9])
title('SSIM of normalized brain (minus tumour) compared to healthy brain (minus tumour)','FontSize',16)
ylabel('SSIM','FontSize',14)

%RMS
s = 1:30; 
%plot for rms within tumour 
figure;
for sub = 1:30
    t = [1:30];
    x = repmat(t,30);
    scatter(x(:,1),rms(s,1),30,'filled');
    hold on;
    scatter(x(:,2),rms(s,3),30,'filled');
    hold on;
    plot([1 2],rms_lesion_norm(sub,[1 2]),'LineWidth',2);
    hold on;
end
box on; grid on
title('RMS of normalized tumour area compared to healthy area','FontSize',16)
ylabel('SSIM','FontSize',14)

%plot for rms of brain minus tumour
figure;
for sub = 1:30
    t = [1:30];
    x = repmat(t,30);
    scatter(x(:,1),rms(s,2),30,'filled');
    hold on;
    scatter(x(:,2),rms(s,4),30,'filled');
    hold on;
    plot([1 2],rms_brain_norm(sub,[1 2]),'LineWidth',2);
    hold on;
end
box on; grid on
title('RMS of normalized brain (minus tumour) compared to healthy brain (minus tumour)','FontSize',16)
ylabel('SSIM','FontSize',14)

% [medianssim,CIssim] = rst_data_plot(SSIM(:,[1 3 2 4]),'estimator','median','newfig','yes');
[diffssim,CIdssim,pssim,alphav,hssim]= rst_multicompare(SSIM,[1 3;2 4],'alphav',0.05,'estimator','median','newfig','yes');
% [medianrms,CIrms] = rst_data_plot(rms(:,[1 3 2 4]),'estimator','median','newfig','yes');
[diffrms,CIdrms,prms,alphav,hrms]= rst_multicompare(rms,[1 3;2 4],'alphav',0.05,'estimator','median','newfig','yes');


      