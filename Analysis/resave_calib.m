%% load subject file and take all but stim pattern to have a lighter file

%% load subj data
[filename,filepath] = uigetfile('*.mat','Select calibration data');
calibration = importdata(fullfile(filepath,filename));


calibration.stim = rmfield(calibration.stim,'pattern');

fpath = '../Data/subj_calib_nopat';
fname = ['NOPAT_' , filename];
save(fullfile(fpath,fname),'calibration');