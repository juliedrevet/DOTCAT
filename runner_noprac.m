%%
close all
clear all java
close all hidden
clc

% ask for subject number, calibration (optional)
prompt={'Subject number (two-digit)','load calibration? (yes/no)'};

argindlg = inputdlg(prompt,'DOTCAT',1,{'','no'});

if isempty(argindlg)
    error('experiment cancelled!');
end
if isempty(argindlg{1})
    error('experiment cancelled!');
end

subj = str2num(argindlg{1}); % subject number
syncflip = true; % synchronize screen refreshs?

% run beginning calibration (or load if already existing)
if strcmp(argindlg{2},'yes')
    [FileName,PathName] = uigetfile('*.mat','Select a calibration file');
    calibration = importdata(fullfile(PathName,FileName));
    if isempty(calibration)
        error('no calibration available');
    end
else
    [calibration,aborted,errmsg] = run_calibration(subj,.2,.8,120);
    fpath = './Data';
    fname = sprintf('DOTCAT_calibration_S%02d_%s',subj,datestr(now,'yyyymmdd-HHMM'));
    fname = fullfile(fpath,fname);
    save([fname,'.mat'],'calibration');
end

expe  = gen_expe(subj);
shape = shuffle_shape(expe); % create subjet shape combinations

expe_raw       = [];
nprac = expe(1).cfg.nprac;
for i = (nprac+1):length(expe)
        expe_raw.expe(i-nprac).cfg.taskid = expe(i).cfg.taskid((nprac+1):end);
        expe_raw.expe(i-nprac).cfg.condtn = expe(i).cfg.condtn((nprac+1):end);
        expe_raw.expe(i-nprac).cfg.epimap = expe(i).cfg.epimap((nprac+1):end);
        expe_raw.expe(i-nprac).cfg.nprac  = 0;
        expe_raw.expe(i-nprac).blck = expe(i).blck;       
end

expe_raw.shape = shape(5:end,:);


foldname = sprintf('./Data/S%02d',subj);
% save raw expe structure in case of premature termination
fname = sprintf('DOTCAT_S%02d_%s_raw',subj,datestr(now,'yyyymmdd-HHMM'));
fname = fullfile(foldname,fname);
if ~exist(foldname,'dir')
    mkdir(foldname);
end
save([fname,'.mat'],'expe_raw');
start_blck = 1;
[expe,aborted] = run_expe(subj,calibration,syncflip,expe_raw,start_blck);
fpath = sprintf('./Data');
fname = sprintf('DOTCAT_S%02d_%s',subj,datestr(now,'yyyymmdd-HHMM'));
fname = fullfile(fpath,fname);
save([fname,'.mat'],'expe');


% run end calibration
[calibration2,aborted,errmsg] = run_calibration(subj,.35,.65,120);
fpath = './Data';
fname = sprintf('DOTCAT_calibration_S%02d_%s',subj,datestr(now,'yyyymmdd-HHMM'));
fname = fullfile(fpath,fname);
save([fname,'.mat'],'calibration2');

fprintf('\nBeginning Calibration:\n')
fprintf('\t minDots / maxDots: %d / %d\n\n',calibration.setup.minDots,calibration.setup.maxDots)
fprintf('End Calibration:\n')
fprintf('\t minDots / maxDots: %d / %d\n\n',calibration2.setup.minDots,calibration2.setup.maxDots)



