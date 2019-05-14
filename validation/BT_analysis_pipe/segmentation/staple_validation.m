%% from the segmentation results, we now only generate 4 sets of images
% GM+WM/GM+WM+CSF * 2 or 3 Gaussians
% it is likely that one of those output is the best but without looking at
% the images, we don't know --> here we validate STAPLE to check if voting
% leads to the best image

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

cd(BRAT_dir);

%% Step 1: run the pipeline for segmentation

index = 1;
spm_jobman('initcfg')
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
        tmp = dir([patient_dir filesep 'VSD.Brain.XX.O.MR_T1.*']);
        T1name = [patient_dir filesep tmp.name filesep 'VSD.nii'];
        tmp = dir([patient_dir filesep 'VSD.Brain_*more*']);
        ground_truth = [patient_dir filesep tmp.name filesep 'VSD.nii'];
        tmp = dir([patient_dir filesep '*Flair*']);
        FLAIRname = [patient_dir filesep tmp.name filesep 'VSD.nii'];
        
        matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = {'/media/cpernet/SATA/BRAT/BRATS2015_Training/HGG/brats_2013_pat0001_1'};
        matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'lesion_masks';
        matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'mask';
        matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{ground_truth}};
        matlabbatch{3}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'T1w';
        matlabbatch{3}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{T1name}};
        matlabbatch{4}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'FLAIR';
        matlabbatch{4}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{FLAIRname}};
        matlabbatch{5}.spm.tools.USwLtools.uswl.imgMsk(1) = cfg_dep('Named File Selector: mask(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{5}.spm.tools.USwLtools.uswl.imgRef(1) = cfg_dep('Named File Selector: T1w(1) - Files', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{5}.spm.tools.USwLtools.uswl.imgStruc(1) = cfg_dep('Named File Selector: T1w(1) - Files', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{5}.spm.tools.USwLtools.uswl.imgStruc(2) = cfg_dep('Named File Selector: FLAIR(1) - Files', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{5}.spm.tools.USwLtools.uswl.imgOth = '';
        matlabbatch{5}.spm.tools.USwLtools.uswl.options.imgTpm = {'/home/cpernet/Documents/MATLAB/MyToolboxes/SPM/spm12/tpm/TPM.nii'};
        matlabbatch{5}.spm.tools.USwLtools.uswl.options.NbGaussian = [2 2 2 2 3 4 2];
        matlabbatch{5}.spm.tools.USwLtools.uswl.options.tpm4lesion = 2;
        matlabbatch{5}.spm.tools.USwLtools.uswl.options.bias.bias_no = 0;
        matlabbatch{5}.spm.tools.USwLtools.uswl.options.ICVmsk = 1;
        matlabbatch{5}.spm.tools.USwLtools.uswl.options.mrf = 2;
        matlabbatch{6}.spm.util.imcalc.input(1) = cfg_dep('US with lesion: c1 image', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c1'));
        matlabbatch{6}.spm.util.imcalc.input(2) = cfg_dep('US with lesion: c2 image', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c2'));
        matlabbatch{6}.spm.util.imcalc.input(3) = cfg_dep('US with lesion: c3 image', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c3'));
        matlabbatch{6}.spm.util.imcalc.input(4) = cfg_dep('US with lesion: c4 image', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c4'));
        matlabbatch{6}.spm.util.imcalc.output = 'maskc3_2GaussiansGMWM';
        matlabbatch{6}.spm.util.imcalc.outdir(1) = cfg_dep('Make Directory: Make Directory ''lesion_masks''', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));
        matlabbatch{6}.spm.util.imcalc.expression = '(i3>i1) & (i3>i2)&(i3>i4)';
        matlabbatch{6}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{6}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{6}.spm.util.imcalc.options.mask = 0;
        matlabbatch{6}.spm.util.imcalc.options.interp = 0;
        matlabbatch{6}.spm.util.imcalc.options.dtype = 4;
        matlabbatch{7}.spm.tools.USwLtools.USwLutils.FxLesMsk.fnMsk(1) = cfg_dep('Image Calculator: ImCalc Computed Image: maskc3_2GaussiansGMWM', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{7}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.minVol = Inf;
        matlabbatch{7}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.fnOth = '';
        matlabbatch{8}.spm.tools.USwLtools.uswl.imgMsk(1) = cfg_dep('Named File Selector: mask(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{8}.spm.tools.USwLtools.uswl.imgRef(1) = cfg_dep('Named File Selector: T1w(1) - Files', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{8}.spm.tools.USwLtools.uswl.imgStruc(1) = cfg_dep('Named File Selector: T1w(1) - Files', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{8}.spm.tools.USwLtools.uswl.imgStruc(2) = cfg_dep('Named File Selector: FLAIR(1) - Files', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{8}.spm.tools.USwLtools.uswl.imgOth = '';
        matlabbatch{8}.spm.tools.USwLtools.uswl.options.imgTpm = {'/home/cpernet/Documents/MATLAB/MyToolboxes/SPM/spm12/tpm/TPM.nii'};
        matlabbatch{8}.spm.tools.USwLtools.uswl.options.NbGaussian = [2 2 2 2 3 4 2];
        matlabbatch{8}.spm.tools.USwLtools.uswl.options.tpm4lesion = 3;
        matlabbatch{8}.spm.tools.USwLtools.uswl.options.bias.bias_no = 0;
        matlabbatch{8}.spm.tools.USwLtools.uswl.options.ICVmsk = 1;
        matlabbatch{8}.spm.tools.USwLtools.uswl.options.mrf = 2;
        matlabbatch{9}.spm.util.imcalc.input(1) = cfg_dep('US with lesion: c1 image', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c1'));
        matlabbatch{9}.spm.util.imcalc.input(2) = cfg_dep('US with lesion: c2 image', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c2'));
        matlabbatch{9}.spm.util.imcalc.input(3) = cfg_dep('US with lesion: c3 image', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c3'));
        matlabbatch{9}.spm.util.imcalc.input(4) = cfg_dep('US with lesion: c4 image', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c4'));
        matlabbatch{9}.spm.util.imcalc.output = 'maskc3_2GaussiansGMWMCSF';
        matlabbatch{9}.spm.util.imcalc.outdir(1) = cfg_dep('Make Directory: Make Directory ''lesion_masks''', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));
        matlabbatch{9}.spm.util.imcalc.expression = '(i3>i1) & (i3>i2)&(i3>i4)';
        matlabbatch{9}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{9}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{9}.spm.util.imcalc.options.mask = 0;
        matlabbatch{9}.spm.util.imcalc.options.interp = 0;
        matlabbatch{9}.spm.util.imcalc.options.dtype = 4;
        matlabbatch{10}.spm.tools.USwLtools.USwLutils.FxLesMsk.fnMsk(1) = cfg_dep('Image Calculator: ImCalc Computed Image: maskc3_2GaussiansGMWMCSF', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{10}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.minVol = Inf;
        matlabbatch{10}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.fnOth = '';
        matlabbatch{11}.spm.tools.USwLtools.uswl.imgMsk(1) = cfg_dep('Named File Selector: mask(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{11}.spm.tools.USwLtools.uswl.imgRef(1) = cfg_dep('Named File Selector: T1w(1) - Files', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{11}.spm.tools.USwLtools.uswl.imgStruc(1) = cfg_dep('Named File Selector: T1w(1) - Files', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{11}.spm.tools.USwLtools.uswl.imgStruc(2) = cfg_dep('Named File Selector: FLAIR(1) - Files', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{11}.spm.tools.USwLtools.uswl.imgOth = '';
        matlabbatch{11}.spm.tools.USwLtools.uswl.options.imgTpm = {'/home/cpernet/Documents/MATLAB/MyToolboxes/SPM/spm12/tpm/TPM.nii'};
        matlabbatch{11}.spm.tools.USwLtools.uswl.options.NbGaussian = [2 2 3 2 3 4 2];
        matlabbatch{11}.spm.tools.USwLtools.uswl.options.tpm4lesion = 2;
        matlabbatch{11}.spm.tools.USwLtools.uswl.options.bias.bias_no = 0;
        matlabbatch{11}.spm.tools.USwLtools.uswl.options.ICVmsk = 1;
        matlabbatch{11}.spm.tools.USwLtools.uswl.options.mrf = 2;
        matlabbatch{12}.spm.util.imcalc.input(1) = cfg_dep('US with lesion: c1 image', substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c1'));
        matlabbatch{12}.spm.util.imcalc.input(2) = cfg_dep('US with lesion: c2 image', substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c2'));
        matlabbatch{12}.spm.util.imcalc.input(3) = cfg_dep('US with lesion: c3 image', substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c3'));
        matlabbatch{12}.spm.util.imcalc.input(4) = cfg_dep('US with lesion: c4 image', substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c4'));
        matlabbatch{12}.spm.util.imcalc.output = 'maskc3_3GaussiansGMWM';
        matlabbatch{12}.spm.util.imcalc.outdir(1) = cfg_dep('Make Directory: Make Directory ''lesion_masks''', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));
        matlabbatch{12}.spm.util.imcalc.expression = '(i3>i1) & (i3>i2)&(i3>i4)';
        matlabbatch{12}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{12}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{12}.spm.util.imcalc.options.mask = 0;
        matlabbatch{12}.spm.util.imcalc.options.interp = 0;
        matlabbatch{12}.spm.util.imcalc.options.dtype = 4;
        matlabbatch{13}.spm.tools.USwLtools.USwLutils.FxLesMsk.fnMsk(1) = cfg_dep('Image Calculator: ImCalc Computed Image: maskc3_3GaussiansGMWM', substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{13}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.minVol = Inf;
        matlabbatch{13}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.fnOth = '';
        matlabbatch{14}.spm.tools.USwLtools.uswl.imgMsk(1) = cfg_dep('Named File Selector: mask(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{14}.spm.tools.USwLtools.uswl.imgRef(1) = cfg_dep('Named File Selector: T1w(1) - Files', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{14}.spm.tools.USwLtools.uswl.imgStruc(1) = cfg_dep('Named File Selector: T1w(1) - Files', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{14}.spm.tools.USwLtools.uswl.imgStruc(2) = cfg_dep('Named File Selector: FLAIR(1) - Files', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{14}.spm.tools.USwLtools.uswl.imgOth = '';
        matlabbatch{14}.spm.tools.USwLtools.uswl.options.imgTpm = {'/home/cpernet/Documents/MATLAB/MyToolboxes/SPM/spm12/tpm/TPM.nii'};
        matlabbatch{14}.spm.tools.USwLtools.uswl.options.NbGaussian = [2 2 3 2 3 4 2];
        matlabbatch{14}.spm.tools.USwLtools.uswl.options.tpm4lesion = 3;
        matlabbatch{14}.spm.tools.USwLtools.uswl.options.bias.bias_no = 0;
        matlabbatch{14}.spm.tools.USwLtools.uswl.options.ICVmsk = 1;
        matlabbatch{14}.spm.tools.USwLtools.uswl.options.mrf = 2;
        matlabbatch{15}.spm.util.imcalc.input(1) = cfg_dep('US with lesion: c1 image', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c1'));
        matlabbatch{15}.spm.util.imcalc.input(2) = cfg_dep('US with lesion: c2 image', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c2'));
        matlabbatch{15}.spm.util.imcalc.input(3) = cfg_dep('US with lesion: c3 image', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c3'));
        matlabbatch{15}.spm.util.imcalc.input(4) = cfg_dep('US with lesion: c4 image', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c4'));
        matlabbatch{15}.spm.util.imcalc.output = 'maskc3_3GaussiansGMWMCSF';
        matlabbatch{15}.spm.util.imcalc.outdir(1) = cfg_dep('Make Directory: Make Directory ''lesion_masks''', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));
        matlabbatch{15}.spm.util.imcalc.expression = '(i3>i1) & (i3>i2)&(i3>i4)';
        matlabbatch{15}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{15}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{15}.spm.util.imcalc.options.mask = 0;
        matlabbatch{15}.spm.util.imcalc.options.interp = 0;
        matlabbatch{15}.spm.util.imcalc.options.dtype = 4;
        matlabbatch{16}.spm.tools.USwLtools.USwLutils.FxLesMsk.fnMsk(1) = cfg_dep('Image Calculator: ImCalc Computed Image: maskc3_3GaussiansGMWMCSF', substruct('.','val', '{}',{15}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{16}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.minVol = Inf;
        matlabbatch{16}.spm.tools.USwLtools.USwLutils.FxLesMsk.options.fnOth = '';
        matlabbatch{17}.spm.util.imcalc.input(1) = cfg_dep('Fixing lesion mask image: Fixed mask image(s)', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{17}.spm.util.imcalc.input(2) = cfg_dep('Fixing lesion mask image: Fixed mask image(s)', substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{17}.spm.util.imcalc.input(3) = cfg_dep('Fixing lesion mask image: Fixed mask image(s)', substruct('.','val', '{}',{13}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{17}.spm.util.imcalc.input(4) = cfg_dep('Fixing lesion mask image: Fixed mask image(s)', substruct('.','val', '{}',{16}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
        matlabbatch{17}.spm.util.imcalc.output = 'sum_of_tmasks';
        matlabbatch{17}.spm.util.imcalc.outdir(1) = cfg_dep('Make Directory: Make Directory ''lesion_masks''', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));
        matlabbatch{17}.spm.util.imcalc.expression = '(i1+i2+i3+i4)';
        matlabbatch{17}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{17}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{17}.spm.util.imcalc.options.mask = 0;
        matlabbatch{17}.spm.util.imcalc.options.interp = 0;
        matlabbatch{17}.spm.util.imcalc.options.dtype = 4;
        seg_outputs{index} = spm_jobman('run', matlabbatch);
        index = index+1;
    end
end

%% Step 2: run Staple


%% Step 3: compare Staple to the best segmentation among the 4


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

% check volumes of the new lesion masks

mask_volumes = NaN(30,12,5);
ground_truth_volumes = NaN(30,1);

subj_index = 1;
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
        
        % find volumes of ground truth for each patient
        tmp                                = dir([patient_dir filesep 'VSD.Brain_*more*']);
        ground_truth                       = [patient_dir filesep tmp.name filesep 'VSD.nii'];
        ground_truth_vol                   = spm_read_vols(spm_vol(ground_truth));
        ground_truth_volumes(subj_index,:) = sum(ground_truth_vol(:)>0);
        writetable(array2table(ground_truth_volumes),[save_in filesep 'ground_truth_volumes.csv']);
        
        %find volume of lesion mask before and after segmentation
        tmp = dir([patient_dir filesep 'VSD.Brain.XX.O.MR_T1.*']);
        param_index = 1;
        for voi = 1:2
            for nbGaussian = 1:3
                for affectedtissue = 1:2
                    mask_folder = ['segmentation_voi' num2str(voi) '_nbG' num2str(Gaussian_param(nbGaussian)) '_tissue' num2str(affectedtissue+1)];
                    % c3_before  = [patient_dir filesep tmp.name filesep mask_folder filesep 'c3kVSD.nii'];
                    c3_thresh1 = [patient_dir filesep tmp.name filesep mask_folder filesep 'c3kVSD_bigger_than_sum_of_others.nii'];
                    c3_thresh2 = [patient_dir filesep tmp.name filesep mask_folder filesep 'c3kVSD_bigger_than_at_least_one.nii'];
                    c3_thresh3 = [patient_dir filesep tmp.name filesep mask_folder filesep 'c3kVSD_bigger_than_each.nii'];
                    c3_thresh4 = [patient_dir filesep tmp.name filesep mask_folder filesep 'c3kVSD_bigger_than_99.nii'];
                    c3_thresh5 = [patient_dir filesep tmp.name filesep mask_folder filesep 'c3kVSD_bigger_than_50.nii'];
                    c3_thresh1_vol = spm_read_vols(spm_vol(c3_thresh1));
                    c3_thresh1_vol = c3_thresh1_vol > 0;
                    c3_thresh2_vol = spm_read_vols(spm_vol(c3_thresh2));
                    c3_thresh2_vol = c3_thresh2_vol > 0;
                    c3_thresh3_vol = spm_read_vols(spm_vol(c3_thresh3));
                    c3_thresh3_vol = c3_thresh3_vol > 0;
                    c3_thresh4_vol = spm_read_vols(spm_vol(c3_thresh4));
                    c3_thresh4_vol = c3_thresh4_vol > 0;
                    c3_thresh5_vol = spm_read_vols(spm_vol(c3_thresh5));
                    c3_thresh5_vol = c3_thresh5_vol > 0;
                    mask_volumes(subj_index,param_index,1) = sum(c3_thresh1_vol(:));
                    mask_volumes(subj_index,param_index,2) = sum(c3_thresh2_vol(:));
                    mask_volumes(subj_index,param_index,3) = sum(c3_thresh3_vol(:));
                    mask_volumes(subj_index,param_index,4) = sum(c3_thresh4_vol(:));
                    mask_volumes(subj_index,param_index,5) = sum(c3_thresh5_vol(:));
                    
                    param_index = param_index+1;
                end
            end
        end
        subj_index = subj_index+1;
    end
    
    param_index = 1;
    mask_names = cell(1,60);
    for th = 1:5
        for voi = 1:2
            for nbGaussian = 1:3
                for affectedtissue = 1:2
                    mask_names{param_index} = ['voi' num2str(voi) '_nbG' num2str(nbGaussian) '_tissue' num2str(affectedtissue+1) '_threshold' num2str(th)];
                    param_index = param_index+1;
                end
            end
        end
    end
    
    mask_vols_reshape = reshape(mask_volumes,[30,60]);
    mask_vols_results = array2table(mask_vols_reshape,'VariableNames',mask_names);
    writetable(mask_vols_results,[save_in filesep 'mask_volumes.csv']);
end

gt_vs_mask = NaN(30,60);
param_index = 1;
for param = 1:60;
    gt_mask = mask_vols_reshape(:,param) - ground_truth_volumes;
    gt_vs_mask(:,param_index) = gt_mask;
    
    param_index = param_index+1;
end

mean_gt_vs_mask = mean(gt_vs_mask,1);
thresh1_gt_mask_mean = mean(mean_gt_vs_mask(:,1:12));
thresh2_gt_mask_mean = mean(mean_gt_vs_mask(:,13:24));
thresh3_gt_mask_mean = mean(mean_gt_vs_mask(:,25:36));
thresh4_gt_mask_mean = mean(mean_gt_vs_mask(:,37:48));
thresh5_gt_mask_mean = mean(mean_gt_vs_mask(:,49:50));

%% step 4: statistically test which masks are the best and in which conditions

[Perf_data, Perf_ranked_data, Perf_adjdata, Perf_ranked_adjdata,cluster_labels] = stats_analysis(save_in);


% calculating correlation between metrics

dice_data  = importdata([save_in filesep 'Dice_results.csv']);
mHd_data   = importdata([save_in filesep 'mHd_results.csv']);
mJ_data    = importdata([save_in filesep 'mJ_results.csv']);
kappa_data = importdata([save_in filesep 'overlap_kappa_results.csv']);
mcc_data   = importdata([save_in filesep 'overlap_mcc_results.csv']);

D = trimmean(dice_data.data ,0.2,'round',2);
H = trimmean(mHd_data.data  ,0.2,'round',2);
J = trimmean(mJ_data.data   ,0.2,'round',2);
M = trimmean(mcc_data.data  ,0.2,'round',2);
K = trimmean(kappa_data.data,0.2,'round',2);

[rg,tg,hg,outidg,hbootg,CIg] = skipped_correlation([D D J],[J H H],1); %compute correlations among global metrics
[rl,tl,hl,outidl,hbootl,CIl] = skipped_correlation(M,K,1); % compute correlation between local metrics

corr_D_J = table(rg(1),tg(1),hg(1),outidg(1),hbootg(1),CIg(1),'VariableNames',{'r','t','h','outid','hboot','CI'});
corr_D_H = table(rg(2),tg(2),hg(2),outidg(2),hbootg(2),CIg(2),'VariableNames',{'r','t','h','outid','hboot','CI'});
corr_J_H = table(rg(3),tg(3),hg(3),outidg(3),hbootg(3),CIg(3),'VariableNames',{'r','t','h','outid','hboot','CI'});
corr_local = table(rl.Spearman,tl.Spearman,hl.Spearman,outidl,hbootl,CIl.Spearman,'VariableNames',{'r','t','h','outid','hboot','CI'});
correlation_results = vertcat(corr_D_J,corr_D_H,corr_J_H,corr_local);
writetable(correlation_results,[save_in filesep 'correlation_results.csv']);
