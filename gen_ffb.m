function ffb_seq = gen_ffb(blck)
% generate acceptable ffb sequence


%% fix RL "arbitrarly" RL parameters
alphabeta = [0.52 7];

%% get block configuration
ntrl      = blck.ntrl;
nfalse    = round((1-blck.p_reward)*ntrl);
stim_cor  = blck.ys; % no ffb
idxrev    = find(blck.switch_seq); % reversal index
nrev      = length(idxrev);
win       = 4; % size of one-sided window around reversal 

%% simulation parameters
nseq   = 5e3;
nsim   = 1e3;
r0     = 1 + (rand(nsim,1,nseq) < 0.5); % random 1st response
sftrnd = rand(nsim,ntrl+1,nseq);        % random softmax

%% generate FFB sequences
stim_ffb = repmat(blck.ys,nseq,1);
for iseq = 1:nseq
    false_seq = sample_falsebis;
    idx_false = find(false_seq);
    % no false positive at beginning of the block (?)
    while idx_false(1) == 1
        false_seq = sample_falsebis;
        idx_false = find(false_seq);
    end
    % inject ffb in stim sequence
    stim_ffb(iseq,false_seq == 1) = 3-stim_ffb(iseq,false_seq == 1);
end
sgnstim_ffb = sign(stim_ffb-1.5);

% evidence x accumulated throughout all ffb sequences and all simulations
x = zeros(nsim,ntrl,nseq);
for iseq = 1:nseq
    x(:,:,iseq) = repmat(sgnstim_ffb(iseq,:),nsim,1);
end

%% Simulate Models
% simulate CatAcc Model
[~,~,~,~,~,~,seqperfsim] = modelsimulation(3,alphabeta,0);

%% Find ffb_seq giving a rev_curve next to the mean rev_curve
% compute zscore for all sequences and take the min zscore
zs = zscore(squeeze(mean(seqperfsim,1)),0,2);
[~,J] = min(mean(zs.^2));


ffb_seq = stim_ffb(J,:) ~= stim_cor;

    function s = sample_falsebis
        if ~mod(nfalse,2)
            s1 = zeros(1,ntrl/2);
            s1(1:nfalse/2) = 1;
            s2 = zeros(1,ntrl/2);
            s2(1:nfalse/2) = 1;
            s = [Shuffle(s1) Shuffle(s2)];
        else
            s1 = zeros(1,ntrl/2);
            s2 = zeros(1,ntrl/2);
            m = randi(2);
            if m == 1
                s1(1:(nfalse+1)/2) = 1;
                s2(1:(nfalse-1)/2) = 1;
            else
                s1(1:(nfalse-1)/2) = 1;
                s2(1:(nfalse+1)/2) = 1;
            end
            s = [Shuffle(s1) Shuffle(s2)];
        end
    end


    function [eperf,mperf,actn,pcor,z,prep,seqperf] = modelsimulation(model,alphabeta,perfrev)
        alpha = alphabeta(1);
        beta  = alphabeta(2);
        
        z = zeros(nsim,ntrl+1,nseq);
        r = zeros(nsim,ntrl+1,nseq);
        r(:,1,:) = r0;
        for t = 1:ntrl
            if model == 1
                z(:,t+1,:) = z(:,t,:) + alpha * (x(:,t,:) - z(:,t,:));
            elseif model == 2
                z(:,t+1,:) = alpha * z(:,t,:) + x(:,t,:);
            elseif model == 3
                z(:,t+1,:) = alpha * z(:,t,:) + sign(x(:,t,:));
            end
            r(:,t+1,:) = 1 + (z(:,t+1,:) > 0);
        end
        actn = (sftrnd>(1./(1+exp(beta*z))))+1;
        pcor = mean(mean(mean(bsxfun(@eq,actn(:,2:end,:),stim_cor),2),1),3);
        
        iprep = zeros(1,nseq);
        for k = 1:nseq
            iprep(k) = sum(sum(diff(actn(:,:,k),1,2)==0))/numel(diff(actn(:,:,k),1,2));
        end
        prep = mean(iprep);
        
        iperf = zeros(nrev,2*win);
        seqperf = zeros(nrev,2*win,nseq);
        
        k = 1;
        for irev = idxrev
            iwin = (irev-win+1):(irev+win);
            iresp   = actn(:,iwin,:);
            icorrect = stim_cor(irev);
            iperf(k,:) = mean(mean(iresp==icorrect),3);
            seqperf(k,:,:) = mean(iresp==icorrect);
            k = k+1;
        end
        
         mperf = mean(iperf,1);
         eperf = mean((mperf-perfrev).^2);
    end

end
