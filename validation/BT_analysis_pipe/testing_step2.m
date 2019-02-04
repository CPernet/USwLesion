     
tmp = 'C:\Users\s1835343\mri_stuff\BRAT\BRATS2015_Training\HGG\brats_2013_pat0001_1\VSD.Brain.XX.O.MR_Flair.54512';
cd(tmp);

% c1 = 2;
% c2 = 3;
% x = [c1,c2];
% sum(x)

C1 = [tmp filesep ['segmentation_' num2str(1)] filesep 'c1kVSD.nii'];
C2 = [tmp filesep ['segmentation_' num2str(1)] filesep 'c2kVSD.nii'];
C3 = [tmp filesep ['segmentation_' num2str(1)] filesep 'c3kVSD.nii'];
C4 = [tmp filesep ['segmentation_' num2str(1)] filesep 'c4kVSD.nii'];

C1_vol = spm_read_vols(spm_vol(C1));
C2_vol = spm_read_vols(spm_vol(C2));
C3_vol = spm_read_vols(spm_vol(C3));
C4_vol = spm_read_vols(spm_vol(C4));

% C1_sum = sum(C1_vol(:)~=0);
% C2_sum = sum(C2_vol(:)~=0);
% C3_sum = sum(C3_vol(:)~=0);
% C4_sum = sum(C4_vol(:)~=0);


x = C3_vol > C1_vol;
y = C3_vol > C2_vol;
z = C3_vol > C4_vol;


normal_tissue = C1_vol + C2_vol + C4_vol;
c3a = (C3_vol > C1_vol) & (C3_vol > C2_vol) & (C3_vol > C4_vol);
c3b = (C3_vol > C1_vol) | (C3_vol > C2_vol) | (C3_vol > C4_vol);

% c3c = sum(C3_vol(:)>0.99) / sum(C3_vol(:)~=0);
c3c = C3_vol > 0.99; 