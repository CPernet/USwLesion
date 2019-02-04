

c1 = 'C:\Users\s1835343\mri_stuff\BRAT\BRATS2015_Training\HGG\brats_2013_pat0004_1\VSD.Brain.XX.O.MR_T1.54531\segmentation_voi1_nbG2_tissue2\c1kVSD.nii';
c2 = 'C:\Users\s1835343\mri_stuff\BRAT\BRATS2015_Training\HGG\brats_2013_pat0004_1\VSD.Brain.XX.O.MR_T1.54531\segmentation_voi1_nbG2_tissue2\c2kVSD.nii';
c3 = 'C:\Users\s1835343\mri_stuff\BRAT\BRATS2015_Training\HGG\brats_2013_pat0004_1\VSD.Brain.XX.O.MR_T1.54531\segmentation_voi1_nbG2_tissue2\c3kVSD.nii';
c4 = 'C:\Users\s1835343\mri_stuff\BRAT\BRATS2015_Training\HGG\brats_2013_pat0004_1\VSD.Brain.XX.O.MR_T1.54531\segmentation_voi1_nbG2_tissue2\c4kVSD.nii';



c1_vol = spm_read_vols(spm_vol(c1));
c2_vol = spm_read_vols(spm_vol(c2));
c3_vol = spm_read_vols(spm_vol(c3));
c4_vol = spm_read_vols(spm_vol(c4));

        if isequal(size(c1_vol),size(c2_vol),size(c3_vol),size(c4_vol));
            disp('dimensions do match');
        else disp('dimensions do not match');
        end
    