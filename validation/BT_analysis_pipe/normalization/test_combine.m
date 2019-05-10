% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'C:\Users\s1835343\mri_stuff\spm12\toolbox\USwLesion\validation\BT_analysis_pipe\normalization\test_combine_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
