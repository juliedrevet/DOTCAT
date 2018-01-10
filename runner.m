%%
close all
clear all java
close all hidden
clc

% ask for subject number, calibration (optional) and raw experiment
% structure (optional)
prompt={'Subject number (two-digit)','load calibration? (yes/no)','load raw experiment structure? (enter starting block)'};

argindlg = inputdlg(prompt,'DOTCAT',1,{'','no','no'});

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


% run experiment (or start at specified block)
if ~strcmp(argindlg{3},'no')
    [FileName,PathName] = uigetfile('*.mat','Select a raw experiment structure');
    expe_raw = importdata(fullfile(PathName,FileName));
    if isempty(expe_raw)
        error('no raw experiment structure available');
    end
    start_blck = str2num(argindlg{3});
    [expe,aborted] = run_expe(subj,calibration,syncflip,expe_raw,start_blck);
    fpath = sprintf('./Data');
    fname = sprintf('DOTCAT_S%02d_%s',subj,datestr(now,'yyyymmdd-HHMM'));
    fname = fullfile(fpath,fname);
    save([fname,'.mat'],'expe');
else
    [expe,aborted] = run_expe(subj,calibration,syncflip);
    fpath = sprintf('./Data');
    fname = sprintf('DOTCAT_S%02d_%s',subj,datestr(now,'yyyymmdd-HHMM'));
    fname = fullfile(fpath,fname);
    save([fname,'.mat'],'expe');
end


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



