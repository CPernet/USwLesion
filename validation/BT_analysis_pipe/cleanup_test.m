index = 1;
tmp = 'C:\Users\s1835343\mri_stuff\BRAT\BRATS2015_Training\HGG\brats_2013_pat0001_1\VSD.Brain.XX.O.MR_Flair.54512';
cd(tmp);
new_folder = mkdir(['segmentation_' num2str(index)]);

files = dir;
for file = 3:size(files,1); 
        
    if strcmp(files(file).name,'c1kVSD.nii') == 1 || strcmp(files(file).name,'c2kVSD.nii') == 1 || ...
            strcmp(files(file).name,'c3kVSD.nii') == 1|| strcmp(files(file).name,'c4kVSD.nii') == 1;
        movefile(files(file).name,['segmentation_' num2str(index)]);
    else
        if strfind(files(file).date,['Feb-2019']);
            disp('delete')
            % delete files(file).name
        end
    end
end
