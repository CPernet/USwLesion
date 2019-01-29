%% Data analysis script the details the different step taken the explain tests


%% deal with datapath; set the directory and defaults

current = pwd;
global defaults
defaults = spm('defaults','FMRI');
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




%% Step 1: for each patient, runs the function image_overlap which computes the similarity between two 3D images based on a number of measures 
% mean Jaccard index, mean Hausdorrf distance, overlap based on voxel/cluster matching -- plus the Dice coefficient

index = 1;
for tumour_type = 1:2
    % loop in High Grade Glioma or Low Grade Glioma
    % --------------------------------------------
    BRAT_dir = eval(['BRAT_dir' num2str(tumour_type)]);
    cd(BRAT_dir); folders = dir;
    if tumour_type == 1
        folders = folders(1:20);
    else
        folders = folders(1:10);
    end
    
    for patient = 1:size(folders,1)-2
        % loop for each patient
        % ---------------------
        patient_dir = [BRAT_dir filesep folders(patient+2).name];
        tmp = dir([patient_dir filesep 'VSD.Brain_*more*']); 
        ground_truth = [patient_dir filesep tmp.name filesep 'VSD.nii'];
        tmp = dir([patient_dir filesep '*Flair*']);
        [mJ1(index),mHd1(index),overlap1(index)] = image_overlap([patient_dir filesep tmp.name filesep 'voi1.nii'],ground_truth);
        Dice1(index) = 2*overlap1(index).voxel.tp / (2*overlap1(index).voxel.tp+2*overlap1(index).voxel.fp+2*overlap1(index).voxel.fn);
        [mJ2(index),mHd2(index),overlap2(index)] = image_overlap([patient_dir filesep tmp.name filesep 'voi2.nii'],ground_truth);
        Dice2(index) = 2*overlap2(index).voxel.tp / (2*overlap2(index).voxel.tp+2*overlap2(index).voxel.fp+2*overlap2(index).voxel.fn);
        
        index = index+1;
    end
end

% figure
x = 1:index;
plot(x,Dice1,'g',x,Dice2,'o');


%% Step 2: call the segmentation routine
% parameters we change are T1name, FLAIRname, voi mask, nbGaussian and affected tissues

count = 1;
Gaussian_param = [2 3 Inf];
for tumour_type = 1:2
    % loop in High Grade Glioma or Low Grade Glioma
    % --------------------------------------------
    BRAT_dir = eval(['BRAT_dir' num2str(tumour_type)]);
    cd(BRAT_dir); folders = dir;
    if tumour_type == 1
        folders = folders(1:20);
    else
        folders = folders(1:10);
    end
       
    for patient = 1:size(folders,1)-2;
        % loop for each patient
        % ---------------------
        patient_dir = [BRAT_dir filesep folders(patient+2).name];
        tmp = dir([patient_dir filesep 'VSD.Brain.XX.O.MR_T1.*']);
        T1name = [patient_dir filesep tmp.name filesep 'VSD.nii'];
        
        tmp = dir([patient_dir filesep '*Flair*']);
        FLAIRname = [patient_dir filesep tmp.name filesep 'VSD.nii']; 
                            
        index = 1;
        for voi = 1:2;
            mask = [patient_dir filesep tmp.name filesep 'voi' num2str(voi) '.nii'];
            for nbGaussian = 1:3;
                for affectedtissue = 1:2 % add +1 for GM+WM or GM+WM+CSF
                    out{index} = segmentation_routine(T1name,FLAIRname,mask,Gaussian_param(nbGaussian),affectedtissue+1);
                    index = index+1;
                end
            end
        end
%         
%         % clean-up
%         mkdir()
%         movefile(out(..)) 
%         delete
%         
    count = count+1;
    end
end