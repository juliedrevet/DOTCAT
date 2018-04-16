%%
close all
clear all java
close all hidden
clc

addpath ./Analysis
% ask for subject number and raw experiment structure

prompt={'Subject number (two-digit)'};
argindlg = inputdlg(prompt,'DOTCAT',1,{''});


if isempty(argindlg)
    error('experiment cancelled!');
end

subj = str2num(argindlg{1}); % subject number

% customizable options:
runbgncalib = true;
runmidcalib = true;
runendcalib = true;
start_blck  = 1;


%% create data folder for subject
foldname = sprintf('./Data/S%02d',subj);
if ~exist(foldname,'dir')
    mkdir(foldname);
end

diary(sprintf('./Data/S%02d/S%02d_Log_%s.txt',subj,subj,datestr(now,'yyyymmdd-HHMM')));

%%
fprintf('*** BEGINNING CALIBRATION ***\n\n')
fprintf('load calibration file...')
[calibfile,calibpath] = uigetfile('*.mat','Select a calibration file');
calibration = importdata(fullfile(calibpath,calibfile));
fprintf('done!\n\n')

if isempty(calibration)
    error('no calibration available');
elseif runbgncalib
    fprintf('Press any key to continue.\n')
    pause
    [calibration1,aborted,errmsg] = run_calibration(subj,calibration,true); % with training!
elseif ~runbgncalib && ~isfield(calibration, 'rslt') 
        fprintf('No calibration results available, run calibration!\n')
        fprintf('Press any key to continue.\n')
        pause
        [calibration1,aborted,errmsg] = run_calibration(subj,calibration,true);
elseif ~runbgncalib && isfield(calibration, 'rslt')
        calibration1 = calibration;   
end
fprintf('\t min / max proportion of Color 1: \t%d / %d\n\n',round(calibration1.rslt.range(1)*100),round(100*calibration1.rslt.range(2)))



fprintf('*** EXPERIMENT 1st Half ***\n\n')
% run experiment (1st half or start at specified block)
fprintf('load raw experiment structure...')
[FileName,PathName] = uigetfile('*.mat','Select a raw experiment structure');
expe_raw = importdata(fullfile(PathName,FileName));
fprintf('done!\n\n')
if isempty(expe_raw)
    error('no raw experiment structure available');
end

fprintf('load raw pattern structure...')
[FileName,PathName] = uigetfile('*.mat','Select a raw pattern structure');
pattern_raw = importdata(fullfile(PathName,FileName));
fprintf('done!\n\n')
if isempty(pattern_raw)
    error('no raw experiment structure available');
end

if start_blck <=6
    fprintf('Press any key to run experiment (1st half).\n')
    pause
    [expe,aborted] = run_expe(subj,calibration1,pattern_raw,expe_raw,start_blck,6);
    start_blck = 7;
end

fprintf('*** MID-CALIBRATION ***\n\n')
if isempty(calibration)
    error('no calibration available');
elseif runmidcalib
    fprintf('Press any key to run mid calibration.\n')
    pause
    [calibration2,aborted,errmsg] = run_calibration(subj,calibration,false);
elseif ~runmidcalib && ~isfield(calibration, 'rslt')
        fprintf('No calibration results available, run calibration!\n')
        fprintf('Press any key to continue.\n')
        pause
        [calibration2,aborted,errmsg] = run_calibration(subj,calibration,false);
elseif ~runmidcalib && isfield(calibration, 'rslt')
        calibration2 = calibration;   
end
fprintf('\t min / max proportion of Color 1: \t%d / %d\n\n',round(calibration2.rslt.range(1)*100),round(100*calibration2.rslt.range(2)))

fprintf('*** EXPERIMENT 2nd Half ***\n\n')
% run experiment (2nd half or start at specified block)
fprintf('Press any key to run experiment (2nd half).\n')
pause
[expe,aborted] = run_expe(subj,calibration2,pattern_raw,expe_raw,start_blck,10);

fprintf('*** END CALIBRATION ***\n\n')
% run end calibration
if isempty(calibration)
    error('no calibration available');
elseif runendcalib
    fprintf('Press any key to run end calibration.\n')
    pause
    [calibration3,aborted,errmsg] = run_calibration(subj,calibration,false); % no training
elseif ~runendcalib && ~isfield(calibration, 'rslt')
        fprintf('No calibration results available, run calibration!\n')
        fprintf('Press any key to continue.\n')
        pause
        [calibration3,aborted,errmsg] = run_calibration(subj,calibration,false);
elseif ~runendcalib && isfield(calibration, 'rslt')
        calibration3 = calibration;   
end

%
plot_calib(calibration);
savefig(sprintf('./Data/S%02d/calibration1',subj))
plot_calib(calibration2);
savefig(sprintf('./Data/S%02d/calibration2',subj))
plot_calib(calibration3);
savefig(sprintf('./Data/S%02d/calibration3',subj))
%
plot_calib(calibration);
plot_calib(calibration2,true);
plot_calib(calibration3,true);
savefig(sprintf('./Data/S%02d/calibration_both',subj))

fprintf('\nBeginning Calibration:\n')
fprintf('\t min / max proportion of Color 1: \t%d / %d\n\n',round(calibration.rslt.range(1)*100),round(100*calibration.rslt.range(2)))
fprintf('Mid Calibration:\n')
fprintf('\t  min / max proportion of Color 1: \t%d / %d\n\n',round(100*calibration2.rslt.range(1)),round(100*calibration2.rslt.range(2)))
fprintf('End Calibration:\n')
fprintf('\t  min / max proportion of Color 1: \t%d / %d\n\n',round(100*calibration3.rslt.range(1)),round(100*calibration3.rslt.range(2)))



diary('off')