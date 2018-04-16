function [calibration] = training_calib
%% generate training calibration patterns

calibration = struct();
calibration.cfg  = struct(); % all configuration parameters needed (default)
calibration.stim = struct(); % stimulation patterns and color proportion presented

%% adapt default configuration structure for training
%cfg = gen_cfg(10,1);
ntrl = 20;
cfg  = gen_cfg(ntrl,1);
% assign NaN to unused fields
cfg.rangeStim = nan; 
cfg.rangePerc = nan;
cfg.mean = nan;
cfg.std  = nan;

calibration.cfg = cfg;

%color_prop = linspace(.3,.7,10);
color_prop = linspace(.3,.7,ntrl);
color_prop = Shuffle(color_prop);
calibration.stim.color_prop = color_prop;

pattern(cfg.nblck,cfg.ntrl) = Pattern();
for iblck = 1:cfg.nblck
    k = 1;
    for p = color_prop(iblck,:)
%         pattern(iblck,k) = gen_pattern(p,cfg);
%         k = k+1;
        pat = gen_pattern(cfg);
        pattern(iblck,k) = color_polygon(pat,p,1e4,1e-4);
        k = k+1;
    end
end

calibration.stim.pattern = pattern;

%% save training calibration structure
fpath = sprintf('./Calibration');
fname = fullfile(fpath,sprintf('training_calibration_%s',datestr(now,'yyyymmdd-HHMM')));
save([fname,'.mat'],'calibration');

end