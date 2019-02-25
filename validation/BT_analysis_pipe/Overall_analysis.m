%% Data analysis script the details the different step taken the explain tests


%% deal with datapath; set the directory and defaults

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

% result dir
toolbox_dir = fileparts(which('crc_USwL.m'));
save_in = [toolbox_dir filesep 'validation' filesep 'BT_analysis_pipe'];

cd(BRAT_dir); save temp_data

%% Step 1: Establish baseline measurement between the ground truth and manual masks
% for each patient, runs the function image_overlap which computes
% the similarity between two 3D images based on a number of measures mean Jaccard index,
% mean Hausdorrf distance, overlap based on voxel/cluster matching -- plus the Dice coefficient

index = 1;
for tumour_type = 1:2
    % loop in High Grade Glioma or Low Grade Glioma
    % --------------------------------------------
    BRAT_dir = eval(['BRAT_dir' num2str(tumour_type)]);
    cd(BRAT_dir); folders = dir;
    if tumour_type == 1
        folders = folders(1:22);
    else
        folders = folders(1:12);
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

% save as mJ1 mJ2 mHd1 mHd2 Dice1 Dice2 MCC1 MCC2 Kappa1 Kappa2

for subject=1:30
        MCC1(subject) = overlap1(subject).voxel.mcc;
        Kappa1(subject) = overlap1(subject).voxel.CK;
        MCC2(subject) = overlap1(subject).voxel.mcc;
        Kappa2(subject) = overlap1(subject).voxel.CK;
end

baseline = table(mJ1', mJ2', mHd1', mHd2', Dice1', Dice2', MCC1', MCC2', Kappa1', Kappa2', ...
    'VariableNames',{'Jaccard1', 'Jaccard2', 'Hausdorff1', 'Hausdorff2', 'Dice1', 'Dice2', 'MatthewCorrCoef1', 'MatthewCorrCoef2', 'KappaCoef1', 'KappaCoef2'});
writetable(baseline,[save_in filesep 'baseline_measures.csv'])

cd(BRAT_dir); save temp_data


%% Step 2: perform segmentation to obtain new tumour masks
% call the segmentation routine in parallel loop (+clean up)
% parameters we change are T1name, FLAIRname, voi mask, nbGaussian and affected tissues

Gaussian_param = [2 3 Inf];
for tumour_type = 1:2
%     % loop in High Grade Glioma or Low Grade Glioma
%     % --------------------------------------------
    BRAT_dir = eval(['BRAT_dir' num2str(tumour_type)]);
    cd(BRAT_dir); folders = dir;
    if tumour_type == 1
        folders = folders(1:22);
    else
        folders = folders(1:12);
    end
    
    parfor patient = 1:size(folders,1)-2
        % loop for each patient
        % ---------------------
        patient_dir = [BRAT_dir filesep folders(patient+2).name];
        tmp = dir([patient_dir filesep 'VSD.Brain.XX.O.MR_T1.*']);
        T1name = [patient_dir filesep tmp.name filesep 'VSD.nii'];
        
        tmp = dir([patient_dir filesep '*Flair*']);
        FLAIRname = [patient_dir filesep tmp.name filesep 'VSD.nii'];
        
        for voi = 1:2
            mask = [patient_dir filesep tmp.name filesep 'voi' num2str(voi) '.nii'];
            for nbGaussian = 1:3
                for affectedtissue = 1:2 % add +1 for GM+WM or GM+WM+CSF
                    out = segmentation_routine(T1name,FLAIRname,mask,Gaussian_param(nbGaussian),affectedtissue+1);
                    % clean up
                    cd(fileparts(T1name))
                    destination = [pwd filesep 'segmentation_voi' num2str(voi) '_nbG' num2str(Gaussian_param(nbGaussian)) '_tissue' num2str(affectedtissue+1)];
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
    end
end

cd(BRAT_dir); save temp_data


%% step 3: compute the similarity between ground truth and masks
mJ    = NaN(30,12,5);
mHd   = NaN(30,12,5);
Dice  = NaN(30,12,5);
mcc   = NaN(30,12,5);
kappa = NaN(30,12,5);

subj_index  = 1;
for tumour_type = 1:2
    BRAT_dir = eval(['BRAT_dir' num2str(tumour_type)]);
    cd(BRAT_dir); folders = dir;
    if tumour_type == 1
        folders = folders(1:22);
    else
        folders = folders(1:12);
    end
    
    for patient = 1:size(folders,1)-2
        patient_dir = [BRAT_dir filesep folders(patient+2).name];
        tmp = dir([patient_dir filesep 'VSD.Brain_*more*']);
        ground_truth = [patient_dir filesep tmp.name filesep 'VSD.nii'];
        tmp = dir([patient_dir filesep 'VSD.Brain.XX.O.MR_T1.*']);

        % set indexing at 1 for each subject
        param_index = 1;
        for voi = 1:2
            for nbGaussian = 1:3
                for affectedtissue = 1:2 % add +1 for GM+WM or GM+WM+CSF
                    % get the files
                    root = [patient_dir filesep tmp.name filesep 'segmentation_voi' num2str(voi) '_nbG' num2str(Gaussian_param(nbGaussian)) '_tissue' num2str(affectedtissue+1)];
                    c1 = [root filesep 'c1kVSD.nii'];
                    c2 = [root filesep 'c2kVSD.nii'];
                    c3 = [root filesep 'c3kVSD.nii'];
                    c4 = [root filesep 'c4kVSD.nii'];
                    
                    % compute similarity
                    [mJ(subj_index,param_index,:),mHd(subj_index,param_index,:),overlap,Dice(subj_index,param_index,:)] = compare2gdtruth(c1,c2,c3,c4,ground_truth,'on');
                    for m=1:5;
                        mcc(subj_index,param_index,m) = overlap(m).voxel.mcc;
                        kappa(subj_index,param_index,m) = overlap(m).voxel.CK;
                    end
                    
                   % update parameter indexing voi*nbGaussian*affectedtissue
                   param_index = param_index+1;
                end
            end
        end
        subj_index = subj_index+1;
    end
end

cd(BRAT_dir); save temp_data

param_index = 1;
Names = cell(1,60);
for th = 1:5
    for voi = 1:2
        for nbGaussian = 1:3
            for affectedtissue = 1:2 % add +1 for GM+WM or GM+WM+CSF
                Names{param_index} = ['voi' num2str(voi) '_nbG' num2str(nbGaussian) '_tissue' num2str(affectedtissue+1) '_threshold' num2str(th)];
                param_index = param_index+1;
            end
        end
    end
end

% save the datasets for each similarity measure as 30x60 tables 

for sim = 1:5
    if sim==1
        sim_name = Dice; sim_nname = 'Dice';
    elseif sim==2
        sim_name = mJ; sim_nname = 'mJ';
    elseif sim==3
        sim_name = mHd; sim_nname = 'mHd';
    elseif sim==4
        sim_name = mcc; sim_nname = 'overlap_mcc';
    elseif sim== 5
        sim_name = kappa; sim_nname = 'overlap_kappa';
    end
    
    M = reshape(sim_name,[30,60]);
    Results = array2table(M,'VariableNames',Names);
    writetable(Results,[save_in filesep sim_nname '_results.csv']);
    
end

% generate similarity results postsegmentation RELATIVE to pre-segmentation (i.e. subtract pre-segmentation from post-segmentation).

for t = 1:5
    if t==1
        tname = ['Dice_results.csv']; tnname = 'Dice';
    elseif t==2
        tname = ['mHd_results.csv']; tnname = 'mHd'; 
    elseif t==3
        tname = ['mJ_results.csv']; tnname = 'mJ'; 
    elseif t==4
        tname = ['overlap_kappa_results.csv']; tnname = 'kappa'; 
    elseif t==5
        tname = ['overlap_mcc_results.csv']; tnname = 'mcc'; 
    end
    
    IMP = importdata([pwd filesep tname]);
    data = IMP.data;
    clear IMP
       
    voi1_col = data(:,[1:6 13:18 25:30 37:42 49:54]); 
    voi2_vol = data(:,[7:12 19:24 31:36 43:48 55:60]);
    
end

%% step 4: statistically test which masks are the best and in which conditions

stats_analysis(save_in)

