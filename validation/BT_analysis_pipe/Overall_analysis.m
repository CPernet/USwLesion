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

MJ = reshape(mJ,[30,60]);
mJ_results = array2table(MJ,'VariableNames',Names);
writetable(mJ_results,[save_in filesep 'mJ_results.csv']);

MHD = reshape(mHd,[30,60]);
mHd_results = array2table(MHD,'VariableNames',Names);
writetable(mHd_results,[save_in filesep 'mHd_results.csv'])

DICE = reshape(Dice,[30,60]);
Dice_results = array2table(DICE,'VariableNames',Names);
writetable(Dice_results,[save_in filesep 'Dice_results.csv'])

MCC = reshape(mcc,[30,60]);
overlap_mcc_results = array2table(MCC,'VariableNames',Names);
writetable(overlap_mcc_results,[save_in filesep 'overlap_mcc_results.csv'])

KAPPA = reshape(kappa,[30,60]);
overlap_kappa_results = array2table(KAPPA,'VariableNames',Names);
writetable(overlap_kappa_results,[save_in filesep 'overlap_kappa_results.csv'])


%% step 4: statistically test which masks are the best and in which conditions


% IMP = importdata([save_in filesep 'mJ_results.csv']); 
% IMP2 = importdata([save_in filesep 'baseline_measures.csv']);
% x = [1:30];
% y1 = IMP.data(:,1); 
% y2 = IMP2.data(:,1);
% scatter(x,y1,'filled');
% hold on;
% scatter(x,y2,'filled');

IMP = importdata([save_in filesep 'baseline_measures.csv']);
x = [1:30];
y = IMP.data(:,1); 
g = gramm('x',x,'y',y);
g.geom_point();
g.stat_glm();
g.set_names('x','Patient','y','mJ');
g.set_title('Mean Jaccard');
g.draw();
