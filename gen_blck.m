function [blck] = gen_blck(cfg, isprac)
%  GEN_BLCK  Generate block of DOTCAT experiment
%
%  Usage: [blck] = GEN_BLCK(cfg)
%
%  where the configuration structure cfg should contain the following fields:
%   * taskid - task identifier => 1:Observer or 2:Agent
%   * condtn - task condition  => 1:single Dot or 2: handful of Dots
%   * epimap - episode mapping => 1:green=target|starting or 2:blue=target|starting
%
%  Output: BLCK structure containing following fields
%   AS IN INPUT CFG:
%   * taskid - task identifier => 1:Observer or 2:Agent
%   * condtn - task condition  => 1:single Dot or 2: handful of Dots
%   * epimap - episode mapping => 1:green=target|starting or 2:blue=target|starting
%   NEW FIELDS:
%   * volatility - block volatility (double)
%   * nepie      - number of episodes
%   * ntrl       - number of trial for this block
%   * p_reward   - reward probability => p_reward in [0,1]
%   * reward_seq - probability of reward at each trial => array(1,ntrl)
%   * switch_seq - reward prob switch => 1: switch or 0: no switch
%   * color_seq  - actual rewarding sequence => 1: epimap color or 0: other color
%   * false_seq  - false positives sequence => 1: false positive or 0: no

if nargin < 1
    error('Missing configuration structure!');
elseif nargin == 1
    isprac = false;
end

% get configuration parameters
taskid = cfg.taskid; % task identifier => 1:observer or 2:agent
condtn = cfg.condtn; % task condition  => 1:single Dot or 2: handful of Dots

gentype = 1; % switch sequence generation 1: based on volatility | 2: base on nepi

% set block parameters
if  ~isprac
    v        = 0.19; % volatility
    nswitch  = 11;    % expected number of switches (should be coherent with volatility and ntrl)
    nepi     = nswitch+1;
    ntrl     = 120;  % number of trials
    vol_csrt = 5;    % wait 5 trials before volatility induces a new switch
    p_reward = 0.80; % 1/5 false positive
    rfalse   = 1-p_reward; % rfalse*ntrl must be integer
else % practice
    v        = 0.12;
    nswitch  = 1;
    nepi     = nswitch+1;
    ntrl     = 30;
    vol_csrt = 5;
    p_reward = 0.90; % 1/10 false positive
    rfalse   = 1-p_reward; % rfalse*ntrl must be integer
end

% create block structure
blck             = [];
blck.taskid      = taskid;
blck.condtn      = condtn;
blck.epimap      = cfg.epimap;
blck.volatility  = v;
blck.nepi        = nepi;
blck.ntrl        = ntrl;
blck.p_reward    = p_reward;

color_seq = zeros(1,ntrl);
skip = false;

while(~(((sum(color_seq)/ntrl)>=0.47)&&((sum(color_seq)/ntrl)<=0.53))&&~skip)
    if isprac
        skip = true;
    else
        skip = false;
    end
    
    %% generate block sequences
    %reward switch sequence {0;1}         
    if gentype == 1
        if isprac
            switch_seq = gen_switch_seq(1000,v);
        else
            switch_seq = gen_switch_seq(1000,v,nswitch);
        end
    else
        if nepi == 2
            b = gen_draws_julie(ntrl/nepi,[vol_csrt ntrl],nepi);
        else
           b = gen_draws_julie(ntrl/nepi,[vol_csrt ntrl/2],nepi);
        end
        sidx = cumsum(b.xs)+1; % idx of switch
        switch_seq = zeros(1,ntrl);
        switch_seq(sidx(1:end-1)) = 1;
    end
    
    %reward probability sequence
    reward_seq = zeros(1,ntrl);
    reward_seq(1:vol_csrt) = p_reward;
    
    for i = vol_csrt+1:ntrl
        if switch_seq(i) == 1
            reward_seq(i) = 1-reward_seq(i-1);
        else
            reward_seq(i) = reward_seq(i-1);
        end
    end
    
    %false positive sequence
    false_seq = sample_false;
    idx_false = find(false_seq);
    % NOT: more than 2 consecutives false positives, false positive at the beginning
    % or at the end of the block
    while any(diff(idx_false(diff(idx_false)==1))==1)||ismember(1,idx_false)||ismember(ntrl,idx_false)... % maybe remove following conditions
            ||any(switch_seq(idx_false))||any(switch_seq(idx_false+1))
        false_seq = sample_false;
        idx_false = find(false_seq);
    end
    
    %actual rewarding sequence
    color_seq = (reward_seq==p_reward)*1.;
    colorffb_seq = color_seq;
    colorffb_seq(false_seq==1) = 1-colorffb_seq(false_seq==1);
    
end

blck.reward_seq   = reward_seq;
blck.switch_seq   = switch_seq;
blck.color_seq    = color_seq;
blck.colorffb_seq = colorffb_seq;
blck.false_seq    = false_seq;


%% nested function to generate switch sequence
    function y = gen_switch_seq(n,v,nswitch)
        % output: switch sequence for volatility v
        if nargin < 3
            nswitch = [];
        end
        
        seq  = zeros(n,ntrl);
        ml   = zeros(n,1); % mean length before switch
        stdl = zeros(n,1); % std length before switch
        
        for m = 1:n
            % generate reward switch sequence
            for l = vol_csrt+1:ntrl
                seq(m,l) =  (rand(1)<=v);
                if (seq(m,l) == 1) && sum(seq(m,(l-vol_csrt):(l-1))~=0)
                    seq(m,l) = 0;
                end
            end
            
            idx = find(seq(m,:));
            
            if isempty(idx)       % nswitch = 0: length is ntrl
                ml(m) = ntrl;
            elseif length(idx)==1 % nswitch = 1: length is the longest part
                ml(m) = max(idx, ntrl-idx);
                stdl(m) = std([idx,ntrl-idx]);
            else                  % nswitch >=2 : mean length taking beginning and end of block into account
                idx2 = [1 idx ntrl];
                ml(m) = mean(diff(idx2));
                stdl(m) = std(diff(idx2));
            end
        end
        
        zs    = zscore(ml);
        if isprac
            idx   = find(stdl<5); %only look for acceptable std
        else
            idx   = find(stdl<10); %only look for acceptable std
        end
        zstd  = zs(idx).';
        [~,J] = min(abs(zstd)); % take the min zscore
        
        y = seq(idx(J),:);
        
        if ~isempty(nswitch) % control switch number
            if sum(y)~=nswitch
                error('Number of episodes doesn''t match, please retry or check parameters! (%d episodes instead of %d)',sum(y)+1,nswitch+1)
            end
        end
        
    end

    function s = sample_false
        if ~mod(int8(rfalse*ntrl),2)
            s = [Shuffle(kron([1 zeros(1,int8(1/rfalse-1))],ones(1,int8(rfalse*ntrl/2))))...
                Shuffle(kron([1 zeros(1,int8(1/rfalse-1))],ones(1,int8(rfalse*ntrl/2))))];
        else % ntrl*rfalse odd
            s = Shuffle([ones(1,int8(ntrl*rfalse)) zeros(1,ntrl*(1-rfalse))]);
            while ((sum(s(1:(ntrl/2)))+1)~=sum(s((ntrl/2)+1:end)))&&(sum(s(1:(ntrl/2)))~=(sum(s((ntrl/2)+1:end))+1))
                s = Shuffle([ones(1,int8(ntrl*rfalse)) zeros(1,ntrl*(1-rfalse))]);
            end
        end
    end

end