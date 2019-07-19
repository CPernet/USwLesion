%% Application

%% set the directory and defaults

current = pwd;
global defaults
defaults = spm('defaults','FMRI');

% data dir
patient_dir = uigetdir(pwd,'Select VBM-BT directory');
cd(patient_dir); local = dir;

%% co-register the Flair/T2 images for each patient to the T1

for patient = 1:54
    cd(local(patient+2).name)
    
    folderinfo = dir(pwd);
    if strfind (folderinfo(3).name,'T1w')
        Ref =folderinfo(3).name;
    else
        Coreg =folderinfo(3).name; 
    end
    
    if strfind (folderinfo(4).name,'T1w')
        Ref =folderinfo(4).name;
    else
        Coreg =folderinfo(4).name;
    end
    
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {Ref};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {Coreg};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    spm_jobman('run', matlabbatch);
    
    cd ..
    
end

%% set origin to AC

for patient = 1:54
    cd(local(patient+2).name)
    
    [t,sts] = spm_select(2,'image','Select files...',{},pwd,'.*','1');
    V = cellstr(t);
    spmup_auto_reorient(V,1);
    
    cd ..
end

%% flip images so the tumour is on the same side (left)

for patient = [3 6:8 14:18 21:22 26:27 29 36 41 45 50 53:54]
    cd(local(patient+2).name)
    
    folderinfo = dir(pwd);
    image_1 = folderinfo(3).name;
    
    if strfind (folderinfo(4).name,'T1w')
        image_2 =folderinfo(4).name;
    elseif strfind (folderinfo(5).name,'T1w')
        image_2 =folderinfo(5).name;
    end
    
    image_1_vol = spm_vol(image_1);
    image_1_read = spm_read_vols(image_1_vol);
    image_2_vol = spm_vol(image_2);
    image_2_read = spm_read_vols(image_2_vol);
    A = flipud(image_1_read);
    B = flipud(image_2_read);
    
    image_1 = spm_write_vol(image_1_vol,A);
    image_2 = spm_write_vol(image_2_vol,B);
    
    cd ..
end

%% USwL normalise/segment
 
s = 1:54; s(40) = [];
for patient = s
    cd(local(patient+2).name)
    
    folderinfo = dir(pwd);
    lesion_mask = folderinfo(3).name;
    image_1 = folderinfo(4).name;
    
    if strfind (folderinfo(5).name,'T1w')
        image_2 =folderinfo(5).name;
    elseif strfind (folderinfo(6).name,'T1w')
        image_2 =folderinfo(6).name;
    end
    
    for nbGaussian = 1:2
        for affectedtissue = 1:2 % add +1 for GM+WM or GM+WM+CSF
            out = USwL_segmentation(image_2,image_1,lesion_mask,nbGaussian,affectedtissue);
            % clean up
            destination = [pwd filesep 'nbG' num2str(nbGaussian) '_tissue' num2str(affectedtissue+1)];
            mkdir(destination)
            movefile(char(out{1}.segmImg.c1),[destination filesep 'c1.nii']);
            movefile(char(out{1}.segmImg.c2),[destination filesep 'c2.nii']);
            movefile(char(out{1}.segmImg.c3),[destination filesep 'c3.nii']);
            movefile(char(out{1}.segmImg.c4),[destination filesep 'c4.nii']);
            movefile(char(out{1}.segmImg.wc1),[destination filesep 'wc1.nii']);
            movefile(char(out{1}.segmImg.wc2),[destination filesep 'wc2.nii']);
            movefile(char(out{1}.segmImg.wc3),[destination filesep 'wc3.nii']);
            movefile(char(out{1}.segmImg.wc4),[destination filesep 'wc4.nii']);
            movefile(char(out{1}.segmImg.mwc1),[destination filesep 'mwc1.nii']);
            movefile(char(out{1}.segmImg.mwc2),[destination filesep 'mwc2.nii']);
            movefile(char(out{1}.segmImg.mwc3),[destination filesep 'mwc3.nii']);
            movefile(char(out{1}.segmImg.mwc4),[destination filesep 'mwc4.nii']);
            movefile(char(out{1}.segmImg.rc1),[destination filesep 'rc1.nii']);
            movefile(char(out{1}.segmImg.rc2),[destination filesep 'rc2.nii']);
            movefile(char(out{1}.segmImg.rc3),[destination filesep 'rc3.nii']);
            movefile(char(out{1}.ICVmsk),[destination filesep 'ICVmsk.nii']);
            movefile(char(out{1}.wICVmsk),[destination filesep 'wICVmsk.nii']);
            movefile(char(out{1}.Struc_1),[destination filesep 'Struc_1.nii']);
            movefile(char(out{1}.Struc_2),[destination filesep 'Struc_2.nii']);
            movefile(char(out{1}.wStruc_1),[destination filesep 'wStruc_1.nii']);
            movefile(char(out{1}.wStruc_2),[destination filesep 'wStruc_2.nii']);
            movefile(char(out{1}.TPMl),[destination filesep 'TPMl.nii']);
            movefile([pwd filesep 'wlesion mask.nii'],[destination filesep 'wlesion mask.nii']);
            movefile([pwd filesep 'swlesion mask.nii'],[destination filesep 'swlesion mask.nii']);
            movefile([pwd filesep 'dtlesion mask.nii'],[destination filesep 'dtlesion mask.nii']);
            movefile([pwd filesep 'BiasField_krsub*.nii'],[destination filesep 'BiasField_alt']);
            movefile([pwd filesep 'BiasField_ksub*.nii'],[destination filesep 'BiasField_T1']);
            movefile([pwd filesep 'iy*.nii'],[destination filesep 'iy']);
            movefile([pwd filesep 'y*.nii'],[destination filesep 'y']);
            movefile([pwd filesep 'krsub*.nii'],[destination filesep 'k_alt']);
            movefile([pwd filesep 'ksub*.nii'],[destination filesep 'k_T1']);
            movefile([pwd filesep 'mkrsub*.nii'],[destination filesep 'mk_alt']);
            movefile([pwd filesep 'mksub*.nii'],[destination filesep 'mk_T1']);
            movefile([pwd filesep 'swicv*.nii'],[destination filesep 'swicv']);
            movefile([pwd filesep 'wmkr*.nii'],[destination filesep 'wmk_alt']);
            movefile([pwd filesep 'wmk*.nii'],[destination filesep 'wmk_T1']);
            movefile([pwd filesep 'ksub*.mat'],[destination filesep 'ksub']);
            delete([pwd filesep 'c5*.nii']);
            delete([pwd filesep 'c6*.nii']);
            
        end
    end
    
    cd ..
end

%threshold and calculate STAPLE

s = 1:54; s(40) = []; s(31) = [];
for patient = s
    cd(local(patient+2).name)
    for nbGaussian = 1:2
        for affectedtissue = 1:2 % add +1 for GM+WM or GM+WM+CSF
            folder = [pwd filesep 'nbG' num2str(nbGaussian) '_tissue' num2str(affectedtissue+1)];
            wc1 = [folder filesep 'wc1.nii']; wc1 = cellstr(wc1);
            wc2 = [folder filesep 'wc2.nii']; wc2 = cellstr(wc2);
            wc3 = [folder filesep 'wc3.nii']; wc3 = cellstr(wc3);
            wc4 = [folder filesep 'wc4.nii']; wc4 = cellstr(wc4);
            masks = [wc1;wc2;wc3;wc4];
            
            matlabbatch{1}.spm.util.imcalc.input = [masks];
            matlabbatch{1}.spm.util.imcalc.output = 'thresholded_mask';
            matlabbatch{1}.spm.util.imcalc.outdir = {folder};
            matlabbatch{1}.spm.util.imcalc.expression = '(i3>i1) & (i3>i2)&(i3>i4)';
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 0;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
            
            matlabbatch{2}.spm.tools.USwLtools.USwLutils.FxLesMsk.fnMsk(1) = cfg_dep('Image Calculator: ImCalc Computed Image: thresholded_mask', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
            matlabbatch{2}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.minVol = Inf;
            matlabbatch{2}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.fnOth = '';
            spm_jobman('run', matlabbatch);
            clear matlabbatch
            
        end
    end
    
    mask_1 = [pwd filesep 'nbG1_tissue2' filesep 'tthresholded_mask.nii']; mask_1 = cellstr(mask_1);
    mask_2 = [pwd filesep 'nbG1_tissue3' filesep 'tthresholded_mask.nii']; mask_2 = cellstr(mask_2);
    mask_3 = [pwd filesep 'nbG2_tissue2' filesep 'tthresholded_mask.nii']; mask_3 = cellstr(mask_3);
    mask_4 = [pwd filesep 'nbG2_tissue3' filesep 'tthresholded_mask.nii']; mask_4 = cellstr(mask_4);
    P = [mask_1;mask_2;mask_3;mask_4];
    crc_STAPLE(P);
    
    cd ..
end

% create sum of (STAPLE) lesion masks for motor/non-motor

i = 1;
for patient = group_1
    cd(local(patient+2).name)
    STAPLE_mask_gr1{i} = [pwd filesep 'nbG1_tissue2' filesep 'tstaple_mask.nii'];
    i = i +1;
    cd ..
end
STAPLE_mask_gr1 = cellstr(STAPLE_mask_gr1');

i = 1;
for patient = group_2
    cd(local(patient+2).name)
    STAPLE_mask_gr2{i} = [pwd filesep 'nbG1_tissue2' filesep 'tstaple_mask.nii'];
    i = i +1;
    cd ..
end
STAPLE_mask_gr2 = cellstr(STAPLE_mask_gr2');

matlabbatch{1}.spm.util.imcalc.input = [STAPLE_mask_gr1];
matlabbatch{1}.spm.util.imcalc.output = 'group1_sum_of_masks';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'sum(X,1)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run', matlabbatch);

matlabbatch{1}.spm.util.imcalc.input = [STAPLE_mask_gr2];
matlabbatch{1}.spm.util.imcalc.output = 'group2_sum_of_masks';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'sum(X,1)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run', matlabbatch);


%% get images for group 1/2

%compute mean of wc1 images

s = 1:54; s(40) = []; s(31) = []; s(36) = [];
for patient = s
    cd(local(patient+2).name)
    wc1_1 = [pwd filesep 'nbG1_tissue2' filesep 'wc1.nii']; wc1_1_vol = spm_vol(wc1_1) ; wc1_1 = spm_read_vols(wc1_1_vol);
    wc1_2 = [pwd filesep 'nbG1_tissue3' filesep 'wc1.nii']; wc1_2_vol = spm_vol(wc1_2) ; wc1_2 = spm_read_vols(wc1_2_vol);
    wc1_3 = [pwd filesep 'nbG2_tissue2' filesep 'wc1.nii']; wc1_3_vol = spm_vol(wc1_3) ; wc1_3 = spm_read_vols(wc1_3_vol);
    wc1_4 = [pwd filesep 'nbG2_tissue3' filesep 'wc1.nii']; wc1_4_vol = spm_vol(wc1_4) ; wc1_4 = spm_read_vols(wc1_4_vol);
    
    wc1_mean = (wc1_1 + wc1_2 + wc1_3 + wc1_4)/4;
    wc1_1_vol.fname = [wc1_1_vol.fname(1:end-20) 'wc1_mean.nii'];
    spm_write_vol(wc1_1_vol,wc1_mean);
    
    %smooth the mean image
    matlabbatch{1}.spm.spatial.smooth.data = {[pwd filesep 'wc1_mean.nii']};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run', matlabbatch);
    clear matlabbatch
    
    cd ..
end

%% get the Total Intracranial Volume

% get volumes
s = 1:54; s(40) = []; 
for patient = s
    cd(local(patient+2).name)
    
    for nbGaussian = 1:2
        for affectedtissue = 1:2 % add +1 for GM+WM or GM+WM+CSF
            
            destination = [pwd filesep 'nbG' num2str(nbGaussian) '_tissue' num2str(affectedtissue+1)];
            movefile([destination filesep '*'],[pwd])
            
            matlabbatch{1}.spm.util.tvol.matfiles = {[pwd filesep 'seg8.mat']};
            matlabbatch{1}.spm.util.tvol.tmax = 4;
            matlabbatch{1}.spm.util.tvol.mask = {'C:\Users\s1835343\mri_stuff\spm12\tpm\mask_ICV.nii,1'};
            matlabbatch{1}.spm.util.tvol.outf = 'volumes';
            spm_jobman('run', matlabbatch);
            clear matlabbatch
            
            movefile([pwd filesep '*'],[destination]);
            movefile([destination filesep 'nbG*'],[pwd]);
            movefile([destination filesep 'lesion mask.nii'],[pwd]);
            movefile([destination filesep 'mwc1_mean.nii'],[pwd]);
            movefile([destination filesep 'smwc1_mean.nii'],[pwd]);
            movefile([destination filesep 'rsub*'],[pwd]);
            movefile([destination filesep 'sub*'],[pwd]);
    
        end
    end
    cd ..
end

% sum volumes to get TIV and take mean across 4 segmentations

TIV = NaN(54,1);
s = 1:54; s(40) = []; s(36) = []; s(31) = []; 
for patient = s
    cd(local(patient+2).name)

    IMP = importdata([pwd filesep 'nbG1_tissue2' filesep 'volumes.csv']); nbG1_tiss2_vols = IMP.data;
    nbG1_tiss2_sum = sum(nbG1_tiss2_vols(:));
    IMP = importdata([pwd filesep 'nbG1_tissue3' filesep 'volumes.csv']); nbG1_tiss3_vols = IMP.data;
    nbG1_tiss3_sum = sum(nbG1_tiss3_vols(:));
    IMP = importdata([pwd filesep 'nbG2_tissue2' filesep 'volumes.csv']); nbG2_tiss2_vols = IMP.data;
    nbG2_tiss2_sum = sum(nbG2_tiss2_vols(:));
    IMP = importdata([pwd filesep 'nbG2_tissue3' filesep 'volumes.csv']); nbG2_tiss3_vols = IMP.data;
    nbG2_tiss3_sum = sum(nbG2_tiss3_vols(:));
    
    TIV(patient,1) = (nbG1_tiss2_sum + nbG1_tiss3_sum + nbG2_tiss2_sum + nbG2_tiss3_sum)/4;
    %index = index+1;
    cd ..
end

cd('C:\Users\s1835343\mri_stuff\spm12\toolbox\USwLesion\validation\BT_analysis_pipe\application')
csvwrite('TIV.csv',TIV);

%% Use SPM summarize to extract data at specific coordinates 

group_1 = [1 4 5 6 7 8 9 10 11 12 14 15 16 17 18 21 22 23 24 25 26 27 29 30 32 34 38 39 43 44 46 49 50 51 52];
group_2 = [2 3 13 19 20 28 33 35 37 41 42 45 47 48 53 54];

group1_TIV = TIV(group_1); group2_TIV = TIV(group_2); 

% MOTOR THALAMUS

index = 1;
Lthal_gr1 = NaN(length(group_1),1);
for patient = group_1
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [Lthal_gr1(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[-19 -15 8]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
Lthal_gr2 = NaN(length(group_2),1);
for patient = group_2
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [Lthal_gr2(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[-19 -15 8]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
Rthal_gr1 = NaN(length(group_1),1);
for patient = group_1
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [Rthal_gr1(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[19 -15 8]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
Rthal_gr2 = NaN(length(group_2),1);
for patient = group_2
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [Rthal_gr2(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[19 -15 8]'),@mean);
    index = index+1;
    cd ..
end

ThData = NaN(length(group_1),4); 
ThData(:,1) = Lthal_gr1; ThData(1:length(group_2),2) = Lthal_gr2;
ThData(:,3) = Rthal_gr1; ThData(1:length(group_2),4) = Rthal_gr2;
[est_thal,HDI_thal] = rst_data_plot(ThData,'estimator','median','newfig','yes');

% MOTOR HANDS (4a/p) 

index = 1;
Lmotor_hands_gr1 = NaN(length(group_1),1);
for patient = group_1
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [Lmotor_hands_gr1(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[-33.2 -26.3 61.5]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
Lmotor_hands_gr2 = NaN(length(group_2),1);
for patient = group_2
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [Lmotor_hands_gr2(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[-33.2 -26.3 61.5]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
Rmotor_hands_gr1 = NaN(length(group_1),1);
for patient = group_1
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [Rmotor_hands_gr1(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[39.1 -22.6 57.2]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
Rmotor_hands_gr2 = NaN(length(group_2),1);
for patient = group_2
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [Rmotor_hands_gr2(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[39.1 -22.6 57.2]'),@mean);
    index = index+1;
    cd ..
end

motorhandsData = NaN(length(group_1),4); 
motorhandsData(:,1) = Lmotor_hands_gr1; motorhandsData(1:length(group_2),2) = Lmotor_hands_gr2;
motorhandsData(:,3) = Rmotor_hands_gr1; motorhandsData(1:length(group_2),4) = Rmotor_hands_gr2;
[est_motor_hands,HDI_motor_hands] = rst_data_plot(motorhandsData,'estimator','median','newfig','yes');

% MOTOR (FEET) - SMA

index = 1;
Lmotor_feet_gr1 = NaN(length(group_1),1);
for patient = group_1
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [Lmotor_feet_gr1(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[-6.6 -9.0 60.9]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
Lmotor_feet_gr2 = NaN(length(group_2),1);
for patient = group_2
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [Lmotor_feet_gr2(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[-6.6 -9.0 60.9]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
Rmotor_feet_gr1 = NaN(length(group_1),1);
for patient = group_1
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [Rmotor_feet_gr1(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[5.6 -8.5 60.9]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
Rmotor_feet_gr2 = NaN(length(group_2),1);
for patient = group_2
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [Rmotor_feet_gr2(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[5.6 -8.5 60.9]'),@mean);
    index = index+1;
    cd ..
end

motorfeetData = NaN(length(group_1),4); 
motorfeetData(:,1) = Lmotor_feet_gr1; motorfeetData(1:length(group_2),2) = Lmotor_feet_gr2;
motorfeetData(:,3) = Rmotor_feet_gr1; motorfeetData(1:length(group_2),4) = Rmotor_feet_gr2;
[est_motor_feet,HDI_motor_feet] = rst_data_plot(motorfeetData,'estimator','median','newfig','yes');

% VI CEREBELLUM

index = 1;
LVI_cerebellum_gr1 = NaN(length(group_1),1);
for patient = group_1
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [LVI_cerebellum_gr1(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[-27.0 -60.0 -28.0]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
LVI_cerebellum_gr2 = NaN(length(group_2),1);
for patient = group_2
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [LVI_cerebellum_gr2(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[-27.0 -60.0 -28.0]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
RVI_cerebellum_gr1 = NaN(length(group_1),1);
for patient = group_1
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [RVI_cerebellum_gr1(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[27.0 -60.0 -28.0]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
RVI_cerebellum_gr2 = NaN(length(group_2),1);
for patient = group_2
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [RVI_cerebellum_gr2(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[27.0 -60.0 -28.0]'),@mean);
    index = index+1;
    cd ..
end

VI_cerebellum_Data = NaN(length(group_1),4); 
VI_cerebellum_Data(:,1) = LVI_cerebellum_gr1; VI_cerebellum_Data(1:length(group_2),2) = LVI_cerebellum_gr2;
VI_cerebellum_Data(:,3) = RVI_cerebellum_gr1; VI_cerebellum_Data(1:length(group_2),4) = RVI_cerebellum_gr2;
[est_VI_cerebellum,HDI_VI_cerebellum] = rst_data_plot(VI_cerebellum_Data,'estimator','median','newfig','yes'); 

% Insula

index = 1;
L_insula_gr1 = NaN(length(group_1),1);
for patient = group_1
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [L_insula_gr1(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[-44.9 -0.5 9.2]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
L_insula_gr2 = NaN(length(group_2),1);
for patient = group_2
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [L_insula_gr2(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[-44.9 -0.5 9.2]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
R_insula_gr1 = NaN(length(group_1),1);
for patient = group_1
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [R_insula_gr1(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[44.9 -0.5 9.2]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
R_insula_gr2 = NaN(length(group_2),1);
for patient = group_2
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [R_insula_gr2(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[44.9 -0.5 9.2]'),@mean);
    index = index+1;
    cd ..
end

Insula_Data = NaN(length(group_1),4); 
Insula_Data(:,1) = L_insula_gr1; Insula_Data(1:length(group_2),2) = L_insula_gr2;
Insula_Data(:,3) = R_insula_gr1; Insula_Data(1:length(group_2),4) = R_insula_gr2;
[est_insula,HDI_insula] = rst_data_plot(Insula_Data,'estimator','median','newfig','yes'); 

% Parietal operculum

index = 1;
L_parietal_gr1 = NaN(length(group_1),1);
for patient = group_1
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [L_parietal_gr1(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[-23.1 -5.3 8.7]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
L_parietal_gr2 = NaN(length(group_2),1);
for patient = group_2
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [L_parietal_gr2(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[-23.1 -5.3 8.7]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
R_parietal_gr1 = NaN(length(group_1),1);
for patient = group_1
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [R_parietal_gr1(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[23.1 -5.3 8.7]'),@mean);
    index = index+1;
    cd ..
end

index = 1;
R_parietal_gr2 = NaN(length(group_2),1);
for patient = group_2
    cd(local(patient+2).name)
    wc1_img = [pwd filesep 'mwc1_mean.nii'];
    [R_parietal_gr2(index,1), xY] = spm_summarise(wc1_img,struct('def','sphere', 'spec',8, 'xyz',[23.1 -5.3 8.7]'),@mean);
    index = index+1;
    cd ..
end

Parietal_Data = NaN(length(group_1),4); 
Parietal_Data(:,1) = L_parietal_gr1; Parietal_Data(1:length(group_2),2) = L_parietal_gr2;
Parietal_Data(:,3) = R_parietal_gr1; Parietal_Data(1:length(group_2),4) = R_parietal_gr2;
[est_parietal,HDI_parietal] = rst_data_plot(Parietal_Data,'estimator','median','newfig','yes'); 


%correct for difference in total intracranial volume
adj_ThData = NaN(length(group_1),4);
[~,~,adj_ThData(:,1)] = regress(ThData(:,1),group1_TIV);
[~,~,adj_ThData(1:length(group_2),2)] = regress(ThData(1:length(group_2),2),group2_TIV);
[~,~,adj_ThData(:,3)] = regress(ThData(:,3),group1_TIV);
[~,~,adj_ThData(1:length(group_2),4)] = regress(ThData(1:length(group_2),4),group2_TIV);

adj_handData = NaN(length(group_1),4);
[~,~,adj_handData(:,1)] = regress(motorhandsData(:,1),group1_TIV);
[~,~,adj_handData(1:length(group_2),2)] = regress(motorhandsData(1:length(group_2),2),group2_TIV);
[~,~,adj_handData(:,3)] = regress(motorhandsData(:,3),group1_TIV);
[~,~,adj_handData(1:length(group_2),4)] = regress(motorhandsData(1:length(group_2),4),group2_TIV);

adj_feetData = NaN(length(group_1),4);
[~,~,adj_feetData(:,1)] = regress(motorfeetData(:,1),group1_TIV);
[~,~,adj_feetData(1:length(group_2),2)] = regress(motorfeetData(1:length(group_2),2),group2_TIV);
[~,~,adj_feetData(:,3)] = regress(motorfeetData(:,3),group1_TIV);
[~,~,adj_feetData(1:length(group_2),4)] = regress(motorfeetData(1:length(group_2),4),group2_TIV);

adj_cerebellData = NaN(length(group_1),4);
[~,~,adj_cerebellData(:,1)] = regress(VI_cerebellum_Data(:,1),group1_TIV);
[~,~,adj_cerebellData(1:length(group_2),2)] = regress(VI_cerebellum_Data(1:length(group_2),2),group2_TIV);
[~,~,adj_cerebellData(:,3)] = regress(VI_cerebellum_Data(:,3),group1_TIV);
[~,~,adj_cerebellData(1:length(group_2),4)] = regress(VI_cerebellum_Data(1:length(group_2),4),group2_TIV);

adj_insulaData = NaN(length(group_1),4);
[~,~,adj_insulaData(:,1)] = regress(Insula_Data(:,1),group1_TIV);
[~,~,adj_insulaData(1:length(group_2),2)] = regress(Insula_Data(1:length(group_2),2),group2_TIV);
[~,~,adj_insulaData(:,3)] = regress(Insula_Data(:,3),group1_TIV);
[~,~,adj_insulaData(1:length(group_2),4)] = regress(Insula_Data(1:length(group_2),4),group2_TIV);

adj_parietalData = NaN(length(group_1),4);
[~,~,adj_parietalData(:,1)] = regress(Parietal_Data(:,1),group1_TIV);
[~,~,adj_parietalData(1:length(group_2),2)] = regress(Parietal_Data(1:length(group_2),2),group2_TIV);
[~,~,adj_parietalData(:,3)] = regress(Parietal_Data(:,3),group1_TIV);
[~,~,adj_parietalData(1:length(group_2),4)] = regress(Parietal_Data(1:length(group_2),4),group2_TIV);

%test for statistical difference between groups at specific motor areas

adj_group1 = horzcat(adj_ThData(:,[1 3]),adj_handData(:,[1 3]),adj_feetData(:,[1 3]),adj_cerebellData(:,[1 3]),adj_insulaData(:,[1 3]),adj_parietalData(:,[1 3]));
adj_group2 = horzcat(adj_ThData(:,[2 4]),adj_handData(:,[2 4]),adj_feetData(:,[2 4]),adj_cerebellData(:,[2 4]),adj_insulaData(:,[2 4]),adj_parietalData(:,[2 4]));

[diff,CI,p,alphav,h]=rst_gp_multicompare(adj_group1,adj_group2,'alphav',0.05,'estimator','median','newfig','yes');

%% finding if tumour vol makes a difference

reordered_TIV = vertcat(group1_TIV,group2_TIV);

group1_tumour_vol = NaN(51,1);
index = 1;
for patient = group_1
    cd(local(patient+2).name)

    IMP = importdata([pwd filesep 'nbG1_tissue2' filesep 'volumes.csv']); nbG1_tiss2_vols = IMP.data;
    nbG1_tiss2_tumvol = nbG1_tiss2_vols(3);
    IMP = importdata([pwd filesep 'nbG1_tissue3' filesep 'volumes.csv']); nbG1_tiss3_vols = IMP.data;
    nbG1_tiss3_tumvol = nbG1_tiss3_vols(3);
    IMP = importdata([pwd filesep 'nbG2_tissue2' filesep 'volumes.csv']); nbG2_tiss2_vols = IMP.data;
    nbG2_tiss2_tumvol = nbG2_tiss2_vols(3);
    IMP = importdata([pwd filesep 'nbG2_tissue3' filesep 'volumes.csv']); nbG2_tiss3_vols = IMP.data;
    nbG2_tiss3_tumvol = nbG2_tiss3_vols(3);
    
    group1_tumour_vol(index,1) = (nbG1_tiss2_tumvol + nbG1_tiss3_tumvol + nbG2_tiss2_tumvol + nbG2_tiss3_tumvol)/4;
    index = index+1;
    cd ..
end
group1_tumour_vol(isnan(group1_tumour_vol))=0;

group2_tumour_vol = NaN(51,1);
index = 36;
for patient = group_2
    cd(local(patient+2).name)

    IMP = importdata([pwd filesep 'nbG1_tissue2' filesep 'volumes.csv']); nbG1_tiss2_vols = IMP.data;
    nbG1_tiss2_tumvol = nbG1_tiss2_vols(3);
    IMP = importdata([pwd filesep 'nbG1_tissue3' filesep 'volumes.csv']); nbG1_tiss3_vols = IMP.data;
    nbG1_tiss3_tumvol = nbG1_tiss3_vols(3);
    IMP = importdata([pwd filesep 'nbG2_tissue2' filesep 'volumes.csv']); nbG2_tiss2_vols = IMP.data;
    nbG2_tiss2_tumvol = nbG2_tiss2_vols(3);
    IMP = importdata([pwd filesep 'nbG2_tissue3' filesep 'volumes.csv']); nbG2_tiss3_vols = IMP.data;
    nbG2_tiss3_tumvol = nbG2_tiss3_vols(3);
    
    group2_tumour_vol(index,1) = (nbG1_tiss2_tumvol + nbG1_tiss3_tumvol + nbG2_tiss2_tumvol + nbG2_tiss3_tumvol)/4;
    index = index+1;
    cd ..
end

group2_tumour_vol(isnan(group2_tumour_vol))=0;
