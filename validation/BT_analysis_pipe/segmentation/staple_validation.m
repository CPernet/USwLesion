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
        
        matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = {patient_dir};
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
        matlabbatch{5}.spm.tools.USwLtools.uswl.options.imgTpm = {[fileparts(which('spm')) filesep 'tpm' filesep 'TPM.nii']};
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
        matlabbatch{8}.spm.tools.USwLtools.uswl.options.imgTpm = {[fileparts(which('spm')) filesep 'tpm' filesep 'TPM.nii']};
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
        matlabbatch{11}.spm.tools.USwLtools.uswl.options.imgTpm = {[fileparts(which('spm')) filesep 'tpm' filesep 'TPM.nii']};
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
        matlabbatch{14}.spm.tools.USwLtools.uswl.options.imgTpm = {[fileparts(which('spm')) filesep 'tpm' filesep 'TPM.nii']};
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
        matlabbatch{15}.spm.util.imcalc.expression = '(i3>i1) & (i3>i2) & (i3>i4)';
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

p = NaN(30,4);
q = NaN(30,4);
jobs = [7 10 13 16];
for patient = 1:30
    for n=1:4
        names{n}=seg_outputs{patient}{jobs(n)}.files{1};
    end
    [~, p(patient,:), q(patient,:)] = crc_STAPLE(names);
end

%% Step 3: compare Staple to the best segmentation among the 4

subject = 1;
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
        tmp = dir([patient_dir filesep 'lesion_masks' filesep 'staple_mask.nii']);
        staple_mask = [patient_dir filesep 'lesion_masks' filesep tmp.name];
        tmp = dir([patient_dir filesep 'VSD.Brain_*more*']);
        ground_truth = [patient_dir filesep tmp.name filesep 'VSD.nii'];
        
        % check staple vs ground truth
        img = spm_read_vols(spm_vol(staple_mask));
        volumes(subject) = sum(img(:)>0);
        [mJ(subject),mHd(subject),overlap] = image_overlap(ground_truth,staple_mask);
        MCC(subject) = overlap.voxel.mcc;
        subject = subject +1;
    end
end

% load previous data from segmentation and put it all together
cd([fileparts(which('crc_STAPLE')) filesep 'validation' filesep 'BT_analysis_pipe' filesep 'segmentation'])
IMP   = importdata('mJ_results.csv');  data  = IMP.data(:,31:34); results(:,1) = mJ; results(:,2) = max(data')';
IMP   = importdata('mHd_results.csv'); data  = IMP.data(:,31:34); results(:,3) = mHd; results(:,4) = min(data')';
IMP   = importdata('overlap_mcc_results.csv');  data  = IMP.data(:,31:34); results(:,5) = MCC; results(:,6) = max(data')';
IMP   = importdata('ground_truth_volumes.csv'); ref_vol = IMP.data;
IMP   = importdata('mask_volumes.csv'); data  = IMP.data(:,31:34); results(:,7) = volumes'-ref_vol; results(:,8) = min([data-ref_vol]')';

% save
staple_results = table(results(:,1),results(:,2),results(:,3),results(:,4),results(:,5),results(:,6), results(:,7),results(:,8),...
    'VariableNames',{'JaccardStaple', 'JaccardBest', 'HausdorffStaple', 'HausdorffBest', 'MCCStaple', 'MCCBest','VolumeStaple', 'VolumeBest'});
writetable(staple_results,'staple_results.csv')

% simple stats
IMP   = importdata('staple_results.csv'); results = IMP.data;
index = 1;figure;
for m=1:4
    subplot(1,4,m);
    [mediansim(m,:),CI(:,:)] = rst_data_plot(results(:,[index index+1]),'estimator','median','newfig','no');
    if m == 1; title('Jaccard');
    elseif m == 2; title('Hausdorff');
    elseif m == 3; title('MCC');
    else title('Volumes')
    end
    index = index+2;
end

[diff,CI,p,alphav,h]= rst_multicompare(results,[1 2;3 4; 5 6;7 8],'alphav',0.05,'estimator','median','newfig','yes')
