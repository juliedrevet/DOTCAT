%% load subject file and take all but stim pattern to have a lighter file

%% load subj data
[filename,filepath] = uigetfile('*.mat','Select subject data');
expe = importdata(fullfile(filepath,filename));


expe.stim = rmfield(expe.stim,'pattern');
fpath = '../Data/subj_data_nopat';
fname = ['NOPAT_' , filename];
save(fullfile(fpath,fname),'expe');