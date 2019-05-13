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
            HV.fname    = [HV.fname(1:end-22) 'Scaled_T1w_skull_stripped.nii'];
            spm_write_vol(HV,healthy_V);
        elseif patient_max > healthy_max
            scaling     = patient_max/healthy_max;
            patient_VSD = patient_VSD./scaling;
            PVSD.fname  = [PVSD.fname(1:end-7) 'Scaled_VSD.nii'];
            spm_write_vol(PVSD,patient_VSD);
        end
        patient_to_segment{subject} = PVSD.fname;
        control_to_use{subject}     = HV.fname;
end


%% tumour segmentation to get lesion masks

% segmentation

Gaussian_param = [2 3];
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
     
    %T1name = patient_to_segment{subject};
    tmp = dir([patient_dir filesep '*Flair*']);
    FLAIRname = [patient_dir filesep tmp.name filesep 'VSD.nii'];
    tmp = dir([patient_dir filesep 'VSD.Brain_*more*']);
    mask = [patient_dir filesep tmp.name filesep 'VSD.nii'];
    for nbGaussian = 1:2
        for affectedtissue = 1:2 % add +1 for GM+WM or GM+WM+CSF
            out = segmentation_routine(T1name,FLAIRname,mask,Gaussian_param(nbGaussian),affectedtissue+1);
            
            % clean up
            cd(fileparts(T1name))
            destination = [pwd filesep 'normalization_segmentation_nbG' num2str(Gaussian_param(nbGaussian)) '_tissue' num2str(affectedtissue+1)];
            mkdir(destination)
            movefile(cell2mat(out{1}.segmImg.c1),[destination filesep 'c1kVSD.nii']);
            movefile(cell2mat(out{1}.segmImg.c2),[destination filesep 'c2kVSD.nii']);
            movefile(cell2mat(out{1}.segmImg.c3),[destination filesep 'c3kVSD.nii']);
            movefile(cell2mat(out{1}.segmImg.c4),[destination filesep 'c4kVSD.nii']);
            delete([pwd filesep 'c5kVSD.nii']);
            delete([pwd filesep 'c6kVSD.nii']);
            delete(cell2mat(out{1}.ICVmsk));
            delete(cell2mat(out{1}.wICVmsk));
            delete(cell2mat(out{1}.Struc_1));
            delete(cell2mat(out{1}.Struc_2));
            delete(cell2mat(out{1}.wStruc_1));
            delete(cell2mat(out{1}.wStruc_1));
            delete(cell2mat(out{1}.TPMl));
            delete(cell2mat(out{1}.segmImg.wc1));
            delete(cell2mat(out{1}.segmImg.wc2));
            delete(cell2mat(out{1}.segmImg.wc3));
            delete(cell2mat(out{1}.segmImg.wc4));
            delete(cell2mat(out{1}.segmImg.mwc1));
            delete(cell2mat(out{1}.segmImg.mwc2));
            delete(cell2mat(out{1}.segmImg.mwc3));
            delete(cell2mat(out{1}.segmImg.mwc4));
            delete(cell2mat(out{1}.segmImg.rc1));
            delete(cell2mat(out{1}.segmImg.rc2));
            delete(cell2mat(out{1}.segmImg.rc3));
            delete([pwd filesep 'BiasField_kVSD.nii']);
            delete([pwd filesep 'iy_kVSD.nii']);
            delete([pwd filesep 'kVSD.nii']);
            delete([pwd filesep 'kVSD_seg8.mat']);
            delete([pwd filesep 'mkVSD.nii']);
            delete([pwd filesep 'swicv_kVSD.nii']);
            delete([pwd filesep 'wmkVSD.nii']);
            delete([pwd filesep 'y_kVSD.nii']);
        end
    end
end

% apply thresholding method 3

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
    else
        folders = folders(1:10);
        patient_dir = [BRAT_dir filesep folders(subject).name];
    end
    
    tmp = dir([patient_dir filesep 'VSD.Brain_*more*']);
    ground_truth = [patient_dir filesep tmp.name filesep 'VSD.nii'];
    tmp = dir([patient_dir filesep 'VSD.Brain.XX.O.MR_T1.*']);
    
    % set indexing at 1 for each subject
    param_index = 1;
    for nbGaussian = 1:2
        for affectedtissue = 1:2 % add +1 for GM+WM or GM+WM+CSF
            % get the files
            root = [patient_dir filesep tmp.name filesep 'normalization_segmentation_nbG' num2str(Gaussian_param(nbGaussian)) '_tissue' num2str(affectedtissue+1)];
            c1 = [root filesep 'c1kVSD.nii'];
            c2 = [root filesep 'c2kVSD.nii'];
            c3 = [root filesep 'c3kVSD.nii'];
            c4 = [root filesep 'c4kVSD.nii'];
            
            % check dimensions of images match
            
            V1 = spm_vol(c1);
            V2 = spm_vol(c2);
            V3 = spm_vol(c3);
            V4 = spm_vol(c4);
            
            if isequal(V1.dim,V2.dim,V3.dim,V4.dim);
                disp('dimensions do match');
            else
                error('image dimensions do not match');
            end
            
            c1_vol = spm_read_vols(V1);
            c2_vol = spm_read_vols(V2);
            c3_vol = spm_read_vols(V3);
            c4_vol = spm_read_vols(V4);
            
            c3_threshold = (c3_vol > c1_vol) & (c3_vol > c2_vol) & (c3_vol > c4_vol);
            c3_threshold = extent_thresold(c3_threshold); %function to keep only the biggest clusters
            [root,filename,ext]=fileparts(V3.fname);
            V3.fname = [root filesep filename '_thresholded' ext];
            V3.descrip = 'c3 > c1 & c3 > c2 & c3 > c4';
            spm_write_vol(V3,c3_threshold);
            
            % update parameter indexing nbGaussian*affectedtissue
            param_index = param_index+1;
        end
    end
    
end

%% create new brains 

% (wT1 - mask) + c3(masked)

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
    else
        folders = folders(1:10);
        patient_dir = [BRAT_dir filesep folders(subject).name];
    end
    
    tmp = dir([patient_dir filesep 'VSD.Brain.XX.O.MR_T1.*']);
        
        param_index = 1;
        for nbGaussian = 1:2
            for affectedtissue = 1:2 % add +1 for GM+WM or GM+WM+CSF
                
                % get the patient mask
                root = [patient_dir filesep tmp.name filesep 'normalization_segmentation_nbG' num2str(Gaussian_param(nbGaussian)) '_tissue' num2str(affectedtissue+1)];
                mask_to_deform = [root filesep 'c3kVSD_thresholded.nii']; 
                
                %apply deformation field to patient mask
                cd(Healthy_dir); local = dir;
                cd(local(subject+2).name)
                deform_field = [pwd filesep 'y_anat.nii'];
                c3_warped = deformation_job(deform_field,mask_to_deform);
                
                %locate new (warped) mask
                c3_mask = spm_vol([root filesep 'rc3kVSD.nii']);
                c3_warped_mask = spm_read_vols(c3_mask);
             
                %get the healthy brain
                healthy_B = spm_vol(control_to_use{subject});
                healthy_brain = spm_read_vols(healthy_B);
                
                %combine healthy brain and patient mask
                test_combine = test_combine_job();
                
                % update parameter index
                param_index = param_index+1;
            end
        end    
end
