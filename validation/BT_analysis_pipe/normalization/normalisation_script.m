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
    Tumour = [patient_dir filesep tmp.name filesep 'wtVSD.nii'];%MNI_tumour{1}.xxx;

    tmp = dir([patient_dir filesep 'VSD.Brain.XX.O.MR_T1.*']);
    Patient = [patient_dir filesep tmp.name filesep 'wmkScaled_VSD.nii'];%MNI_patient{1}.xxx;
    
    V1 = spm_vol(T1w);    T1w = spm_read_vols(V1);
    V2 = spm_vol(Tumour); Tumour = spm_read_vols(V2);
    V3 = spm_vol(Patient);Patient = spm_read_vols(V3);
    mask = Tumour > 0.1;
    T1w(mask) = Patient(mask);
    cd(Healthy_dir); local = dir;
    cd(local(subject+2).name)  
    V1.fname = [pwd filesep 'T1w_with_tumour_TEST.nii'];
    spm_write_vol(V1,T1w);

    %get the correct bounding box
    matlabbatch{1}.spm.util.bbox.image = {'C:\Users\s1835343\mri_stuff\BRAT\F1000\0025337\anat.nii,1'};
    matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
    out = spm_jobman('run', matlabbatch);
    clear matlabbatch

    %inverse normalise the healthy brain+tumour back into subject space
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[pwd filesep 'iy_anat.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[pwd filesep 'T1w_with_tumour.nii']};%{'C:\Users\s1835343\mri_stuff\BRAT\BRATS2015_Training\LGG\brats_2013_pat0002_1\VSD.Brain.XX.O.MR_T1.54633\Scaled_VSD.nii,1'};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = out{1}.bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    spm_jobman('run', matlabbatch);
    clear matlabbatch
    
    %inverse normalise the tumour mask into healthy patient subject space
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[pwd filesep 'iy_anat.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {Tumour};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = out{1}.bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w_HV_';
    spm_jobman('run', matlabbatch);
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

% standard SPM segmentation of wT1w_skull_stripped.nii from healthy subject

for subject = 1:30
    cd(Healthy_dir); local = dir;
    cd(local(subject+2).name)   
    wT1w_skull_stripped = [pwd filesep 'wT1w_skull_stripped.nii']; 
    brain_mask = [pwd filesep 'wbrain_mask.nii']; %brain mask from healthy subject
    wT1w_with_tumour = [pwd filesep 'wT1w_with_tumour.nii']; 
    
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
    inv_wtVSD = [patient_dir filesep tmp.name filesep 'ww_HV_wtVSD.nii'];
    
    cd(Healthy_dir); local = dir;
    cd(local(subject+2).name)
    
    %standard SPM segmentation of wT1w_skull_stripped.nii from healthy subject
    standard_H =standard_spm_segment_job(wT1w_skull_stripped);
    image_1 = [pwd filesep 'healthy_standard_segmentation' filesep 'c1wT1w_skull_stripped.nii'];
    %[image_1{1,1}.tiss.c]
    
    %standard SPM segmentation of wT1w_with_tumour.nii
    
    standard_HT =standard_spm_segment_job(wT1w_with_tumour);
    image_2 = [pwd filesep 'healthy_tumour_standard_segmentation' filesep 'c1wT1w_with_tumour.nii'];
    
    %USwL segmentation of wT1w_with_tumour.nii
    seg_with_lesion_HT = segment_with_lesion(inv_wtVSD,wT1w_with_tumour);

    %image 1 vs image 2
    %measure similarity of whole brain
    SSIM_1 = SSI(image_1,image_2,brain_mask,1);
    
    %do we get the same values per voxel?
    rms_1 = sqrt(mean((c1image1-c1image2).^2));
    
    %image 1 vs image 3
    %measure similarity of whole brain
    SSIM_2 = SSI(image1,image3,roi,1);
    
    %do we get the same values per voxel?
    rms_2 = sqrt(mean((c1image1-c1image2).^2));

end




image1 --> standard SPM segmentation of wT1w_skull_stripped.nii from healthy subject
roi --> brain_mask.nii from healthy subject
image2 --> standard SPM segmentation of wT1w_with_tumour.nii
SSIM = SSI(image1,image2,roi,1) % how similar the whole brain is
rms = sqrt(mean((c1image1-c1image2).^2)); % do we get the same values per voxel
4 times
image3 --> segmentation with lesions of wT1w_with_tumour.nii
SSIM = SSI(image1,image3,roi,1)
rms = sqrt(mean((c1image1-c1image3).^2));


      