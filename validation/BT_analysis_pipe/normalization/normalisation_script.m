%% Normalisation

%% set the directory and defaults

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

Healthy_dir = uigetdir(pwd,'Select F1000 directory');

%% run the normalization on healthy brains

% SPM segment into standard space

cd(Healthy_dir); local = dir;
for subject = 1:30
     cd(local(subject+2).name)
    normalise = normalisation_batch_job();   
    cd ..
end

% intensity normalisation

if min(img(:)) < 0
    img = img + min(img(:));
else
    img = img - min(img(:));
end 

n = max(img(:)); 
if n < 500
    scaling = 500/n;
elseif n > 500
    scaling = n/500;
end
img = img.*scaling;

