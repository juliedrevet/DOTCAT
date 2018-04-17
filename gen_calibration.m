function [calibration] = gen_calibration(ntrl,nblck,fname)
%  GEN_CALIBRATION  Generate raw calibration structure
%
%  Usage: [calibration] = GEN_CALIBRATION(ntrl,nblck,fname)
%
%  with argument:
%   * ntrl  - number of trials for one block 
%   * nblck - number of blocks (should be 2 or at least even)
%   * fname - filename to save calibration structure (optional)
%
%  Output will be saved and should be used for run_calibration.m

if nargin<1
    ntrl = 120;
    nblck = 2;
end

calibration = struct();
calibration.cfg  = struct(); % all configuration parameters needed
calibration.stim = struct(); % stimulation patterns and color proportion presented

%% Generate configuration structure for calibration
cfg = gen_cfg(ntrl,nblck); % take default configuration
calibration.cfg = cfg;

%% Generate array of color proportion to apply to each stim pattern
% constraints: 
%       * proportions within rangeStim
%       * consecutive proportions >50% (<50%) w.r.t. pseudornd
%       NOT* remove proportions too close to 50%
%       * proportions <50% and >50% balanced accross the 2 blocks

color_prop = zeros(nblck,ntrl);
fprintf('Generating color proportions...')
for iblck = 1:nblck
    % generate indexing to control balancing constraints
    idx = ones(1,ntrl);
    while abs(diff(hist(idx,1:2)))~=0
        idx = ones(1,ntrl);
        idx(1:cfg.pseudornd) = ceil(2*rand(1,cfg.pseudornd));
        for i = (cfg.pseudornd+1):ntrl
            iidx = ceil(2*rand(1));
            if sum(idx((i-cfg.pseudornd):(i-1))==iidx)>=cfg.pseudornd
                idx(i) = 3-iidx;
            else
                idx(i) = iidx;
            end
        end
    end
        
    % generate proportions from gaussian distribution
    prop = 0;
    while length(prop)<4*ntrl
        % generate much more than necessary and remove step by step
        prop = normrnd(cfg.mean,cfg.std,10*ntrl,1);
        % remove proportions too close to cfg.mean
        %prop = prop((prop <= cfg.mean-.025)|(prop >= cfg.mean+.025));
        % remove proportions < rangeStim(1)
        prop = prop(prop > cfg.rangeStim(1));
        % remove proportions > rangeStim(2)
        prop = prop(prop < cfg.rangeStim(2));
    end
    prop = Shuffle(prop);
    prop1 = prop(prop<.5);
    prop2 = prop(prop>.5);
    color_prop(iblck,idx == 1) = prop1(1:sum(idx == 1));
    color_prop(iblck,idx == 2) = prop2(1:sum(idx == 2));  
end
fprintf('done!\n')
calibration.stim.color_prop = color_prop;

% % plot color proportion distribution
% h = histogram(color_prop(:),cfg.rangeStim(1):0.025:cfg.rangeStim(2));
% xlim([0 1]);
% xticks([0 cfg.rangeStim(1) mean([cfg.rangeStim(1) cfg.mean]) cfg.mean  mean([cfg.rangeStim(2) cfg.mean]) cfg.rangeStim(2) 1])
% title('Proportion of color 1 in DOT pattern')é
% xlabel('proportion')
% ylabel('#occurence')

%% Generate stim patterns according to color proportion
fprintf('Generating patterns...')
pattern(nblck,ntrl) = Pattern();
for iblck = 1:nblck
    k = 1;
    for p = color_prop(iblck,:)
        pat = gen_pattern(cfg);
        pattern(iblck,k) = color_polygon(pat,p,1e4,1e-4);
        k = k+1;
    end
end
fprintf('done!\n')

calibration.stim.pattern = pattern;

%% save raw calibration structure (to be used for run_calibration)
if nargin<3
    fpath = sprintf('./Calibration/Raw');
    fname = fullfile(fpath,sprintf('DOTCAT_calibration_raw_%s_%s',cfg.type,datestr(now,'yyyymmdd-HHMM')));
end
save([fname,'.mat'],'calibration');
end