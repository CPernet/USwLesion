% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'C:\Users\s1835343\mri_stuff\spm12\toolbox\USwLesion\normalization\normalisation_batch_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
