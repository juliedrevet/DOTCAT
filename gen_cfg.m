function cfg = gen_cfg(ntrl,nblck)
% generate default configuration structure needed for gen_calibration

if nargin<1
    ntrl = 120;
    nblck = 2;
end

cfg = struct();

% calibration experiment parameters
cfg.ntrl  = ntrl;
cfg.nblck = nblck; % counterbalance left|right buttons for classification

% pattern generation parameters
cfg.dmtr = 160;      % DOT diameter
cfg.nang = 22;       % nb of angles to approximate DOT with polygon
cfg.xy   = [0,0];    % center of pattern coordinates
cfg.type = 'size';   % type of pattern ('size'>> based on polygon size, else based on area)
cfg.dpth = 30;       % depth of recursion for binary space partitioning algorithm

% color proportions parameters 
cfg.pseudornd = 8; % max consecutive proportions >50%, resp. <50% (pseudo-random generation)
cfg.rangeStim = [.35 .65];%[.4 .6];% % only show color proportions within this range
cfg.rangePerc = [.2 .8]; % range of perception uncertainty (to be matched with false feedbacks rate in single dot condition!)
cfg.mean = .5; % mean of Gaussian distribution
cfg.std  = .05; % standard deviation of Gaussian distribution

end