function [mJ,mHd,overlap,Dice] = compare2gdtruth(c1,c2,c3,c4,ground_truth,save_option);

% FORMAT: [mJ,mHd,overlap,Dice] = compare2gdtruth(c1,c2,c3,c4,ground_truth,save_option);
%
% INPUT c1, c2, c3, c4 
%          ground_truth
%        save_option 'on' or 'off' 
%

% OUTPUT mJ,mHd,overlap,Dice
 
% default
save_option = 'off';

%% 1st load images and check all dimension match

c1 = [patient_dir filesep tmp.name filesep ['segmentation_' num2str(index)] filesep 'c1kVSD.nii'];
c2 = [patient_dir filesep tmp.name filesep ['segmentation_' num2str(index)] filesep 'c2kVSD.nii'];
c3 = [patient_dir filesep tmp.name filesep ['segmentation_' num2str(index)] filesep 'c3kVSD.nii'];
c4 = [patient_dir filesep tmp.name filesep ['segmentation_' num2str(index)] filesep 'c4kVSD.nii'];

c1_vol = spm_read_vols(spm_vol(c1));
c2_vol = spm_read_vols(spm_vol(c2));
c3_vol = spm_read_vols(spm_vol(c3));
c4_vol = spm_read_vols(spm_vol(c4));

if isequal(size(c1_vol),size(c2_vol),size(c3_vol),size(c4_vol)); 
    disp('dimensions do match');
else disp('dimensions do not match');
end

%% create binary images of c3

% tumour prob > sum of all other classes 
% normal tissue is c1+c2+c4

normal_tissue = c1_vol + c2_vol + c4_vol;
c3a = c3_vol > normal_tissue;

if strcmp(save_option,'on')
    spm_write_vol
end

% tumour prob > any of other classes 

c3b = (c3_vol > c1_vol) | (c3_vol > c2_vol) | (c3_vol > c4_vol);

if strcmp(save_option,'on')
    spm_write_vol
end

% tumour > 99%

c3c = sum(c3_vol(:)>0.99) / sum(c3_vol(:)~=0);

if strcmp(save_option,'on')
    spm_write_vol
end


clean c1 c2 c3 c4

%% now compute overlap

GT = spm_read_vols(spm_vol(ground_truth));

[mJ(1),mHd(1),overlap(1)] = image_overlap(c3a,GT);
Dice(1) = 2*overlap1(1).voxel.tp / (2*overlap1(1).voxel.tp+2*overlap1(1).voxel.fp+2*overlap1(1).voxel.fn);

[mJ(2),mHd(3),overlap(2)] = image_overlap(c3b,GT);
Dice(2) = 2*overlap1(2).voxel.tp / (2*overlap1(2).voxel.tp+2*overlap1(2).voxel.fp+2*overlap1(2).voxel.fn);

[mJ(3),mHd(3),overlap(3)] = image_overlap(c3c,GT);
Dice(3) = 2*overlap1(3).voxel.tp / (2*overlap1(3).voxel.tp+2*overlap1(3).voxel.fp+2*overlap1(3).voxel.fn);

