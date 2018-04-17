function preconfig(subj)

% create subj folder
foldname = sprintf('./Raw/S%02d',subj);
if ~exist(foldname,'dir')
    mkdir(foldname);
end

tic
% generate raw calibration structure
fname = sprintf('DOTCAT_S%02d_%s_raw_calib',subj,datestr(now,'yyyymmdd-HHMM'));
fname = fullfile(foldname,fname);
calibration = gen_calibration(120,2,fname);
toc

cfg = calibration.cfg;

% generate raw experiment structure and save in subject folder
tic
expe  = gen_expe(subj);
shape = shuffle_shape(expe); % create subject shape combinations
expe_raw       = [];
expe_raw.expe  = expe;
expe_raw.shape = shape;
toc
fname = sprintf('DOTCAT_S%02d_%s_raw_expe',subj,datestr(now,'yyyymmdd-HHMM'));
fname = fullfile(foldname,fname);
save([fname,'.mat'],'expe_raw');

nblck = length(expe.blck);

% generate raw patterns (no color proportion applied yet) and save
pattern(nblck,1) = Pattern();
tic
for iblck = 1:nblck
    condtn = expe.blck(iblck).condtn;
    ntrl   = expe.blck(iblck).ntrl;
    if condtn == 2
        for itrl = 1:ntrl
            pattern(iblck,itrl) = gen_pattern(cfg);
        end
    else
        pattern(iblck,:) = Pattern();
    end
end
toc
fname = sprintf('DOTCAT_S%02d_%s_raw_pat',subj,datestr(now,'yyyymmdd-HHMM'));
fname = fullfile(foldname,fname);
save([fname,'.mat'],'pattern');


end