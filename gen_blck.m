function [blck] = gen_blck(isprac)
%  GEN_BLCK  Generate block of DOTCAT experiment
%
%  Usage: [blck] = GEN_BLCK(cfg)
%
%  Output: BLCK structure containing following fields
%   NEW FIELDS:
%   * volatility - block volatility (double)
%   * nepi       - number of episodes
%   * ntrl       - number of trial for this block
%   * p_reward   - reward probability => p_reward in [0,1]
%   * reward_seq - probability of reward at each trial => array(1,ntrl)
%   * switch_seq - reward prob switch => 1: switch or 0: no switch
%   * color_seq  - actual rewarding sequence => 1: epimap color or 0: other color
%   * false_seq  - false positives sequence  => 1: false positive or 0: no

if nargin < 1
    isprac = false;
end

% set block parameters
if  ~isprac
    davg     = 10;
    nepi     = 12;
    p_reward = 0.80; % 1/5 false positive
else % practice
    davg     = 15;
    nepi     = 2;    
    p_reward = 0.90; % 1/10 false positive
end

ntrl   = davg * nepi;  % number of trials per volatility level
dlim   = [5 ntrl];     % min|max number of trials before reversal

% generate episodes (Valentin's function)
b = gen_epi(davg,dlim,nepi);

% create block structure
blck             = b;
blck.ntrl        = ntrl;
blck.p_reward    = p_reward;

blck = orderfields(blck,['ntrl'; fieldnames(b);'p_reward']);

% switch sequence
idx_switch = cumsum(b.xs)+1;
switch_seq = zeros(1,ntrl);
switch_seq(idx_switch(1:end-1)) = 1;

% reward probablity sequence
reward_seq = zeros(1,ntrl);
reward_seq(b.ys==1) = p_reward;
reward_seq(b.ys~=1) = 1-p_reward;

%%
blck.reward_seq   = reward_seq;
blck.switch_seq   = switch_seq;

%false positive sequence
if ~isprac
    false_seq = gen_ffb(blck);
else % fix false positives to go faster
    false_seq = zeros(1,ntrl);
    false_seq([8 20 25]) = 1;
end

%actual rewarding sequence
color_seq = b.ys;
color_seq(color_seq == 2) = 0;
colorffb_seq = color_seq;
colorffb_seq(false_seq==1) = 1-colorffb_seq(false_seq==1);


blck.color_seq    = color_seq;
blck.colorffb_seq = colorffb_seq;
blck.false_seq    = false_seq;

end