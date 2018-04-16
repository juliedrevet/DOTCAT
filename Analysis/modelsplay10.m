function modelsplay10
% look if the reversal curve generated by a CatAcc Model with "randomly"
% fixed alpha and beta and internal FFB is exactly the same as the one that
% can be generated with a CatAcc Model with explicit FFB (without any
% constraints or with coarse constraints) using the same alpha and beta
% LONG RUNNING TIME!!
% CONCLUSION: we take FFB sequences balanced but we abandon the constraint
% "max 2 consecutives"

addpath ../

davg = 10;
nepi = 12;
ntrl = davg * nepi;
dlim   = [5 ntrl];     % min|max number of trials before reversal
blck = gen_epi(davg,dlim,nepi);

% subj
% 0.0795    0.0682    0.0511    0.0341    0.5795    0.8068    0.8807    0.9091


%% customizable option
alphabeta1 = [0.54 7];
% alphabeta2 = [0.46 7];
% alphabeta3 = [0.46 2];
% alphabeta4 = [0.46 100];

%% get block configuration
ntrl     = length(blck.ys);
rfalse   = .2;
stim_cor = blck.ys; % no ffb
idxrev   = cumsum(blck.xs)+1; % reversal index
idxrev   = idxrev(1:end-1);
nrev     = length(idxrev);
win      = 4;
nfalse   = rfalse*ntrl;

%% simulation parameters
mu     = abs(fzero(@(x)(normcdf(x)-rfalse),-1)); % x drawn from normal distribution mean +/-mu
nseq   = 1e4;
nsim   = 1e3;
r0     = 1 + (rand(nsim,1,nseq) < 0.5); % random 1st response
sftrnd = rand(nsim,ntrl+1,nseq);        % random softmax
x      = zeros(nsim,ntrl,nseq);

%% generate CatAcc Model simulation with evidence randn (internal FFB)
stim_ffb    = repmat(stim_cor,nseq,1); %take seq without ffb and add nseq various ffb sequences
sgnstim_ffb = sign(stim_ffb-1.5);
for iseq = 1:nseq
    x(:,:,iseq) = repmat(sgnstim_ffb(iseq,:),nsim,1)*mu + randn(nsim,ntrl);
end
% compute reversal curve
[~,perfrev1sim,~,~,~,~,~] = modelsimulation(3,alphabeta1,0);

%% generate CatAcc Model simulation with evidence dizcretized (external FFB NO constraints)
for iseq = 1:nseq
    false_seq = zeros(1,ntrl);
    false_seq(1:rfalse*ntrl) = 1;
    false_seq = Shuffle(false_seq);
    % inject ffb in stim sequence
    stim_ffb(iseq,false_seq == 1) = 3-stim_ffb(iseq,false_seq == 1);
end
sgnstim_ffb = sign(stim_ffb-1.5);
for iseq = 1:nseq
    x(:,:,iseq) = repmat(sgnstim_ffb(iseq,:),nsim,1);
end
% compute reversal curve
[~,perfrev2sim,~,~,~,~,~] = modelsimulation(3,alphabeta1,0);

%% generate CatAcc Model simulation with evidence dizcretized (external FFB balanced)
stim_ffb    = repmat(stim_cor,nseq,1); %take seq without ffb and add nseq various ffb sequences
for iseq = 1:nseq
    false_seq = sample_falsebis;
    % inject ffb in stim sequence
    stim_ffb(iseq,false_seq == 1) = 3-stim_ffb(iseq,false_seq == 1);
end
sgnstim_ffb = sign(stim_ffb-1.5);
for iseq = 1:nseq
    x(:,:,iseq) = repmat(sgnstim_ffb(iseq,:),nsim,1);
end
% compute reversal curve
[~,perfrev3sim,~,~,~,~,~] = modelsimulation(3,alphabeta1,0);

%% generate CatAcc Model simulation with evidence dizcretized (external FFB balanced + 2 consecutive max)
stim_ffb    = repmat(stim_cor,nseq,1); %take seq without ffb and add nseq various ffb sequences
for iseq = 1:nseq
    false_seq = sample_false;
    idx_false = find(false_seq);
    % NOT: more than 2 consecutives false positives, false positive at the beginning
    % or at the end of the block
    while any(diff(idx_false(diff(idx_false)==1))==1)||ismember(1,idx_false)||ismember(ntrl,idx_false)
        false_seq = sample_false;
        idx_false = find(false_seq);
    end
    % inject ffb in stim sequence
    stim_ffb(iseq,false_seq == 1) = 3-stim_ffb(iseq,false_seq == 1);
end
sgnstim_ffb = sign(stim_ffb-1.5);
for iseq = 1:nseq
    x(:,:,iseq) = repmat(sgnstim_ffb(iseq,:),nsim,1);
end
% compute reversal curve
[~,perfrev4sim,~,~,~,~,~] = modelsimulation(3,alphabeta1,0);



%% plot
c = lines;
figure()
hold on
plot(perfrev2sim,'LineWidth',2,'Color',c(1,:))
plot(perfrev3sim,'LineWidth',2,'Color',c(2,:))
plot(perfrev4sim,'LineWidth',2,'Color',c(3,:))
plot(perfrev1sim,'k.','MarkerSize',20)
ylabel('p(corr)')
xticklabels({'-4','','-2','','','2','','4'})
xlabel('trial position from reversal')
set(gca,'FontSize',20);
line([4.5 4.5],[0 1],'LineStyle','--','Color','k')
text(4.5,0.9,'reversal','HorizontalAlignment','center','FontAngle','italic','FontSize',16);
legend({'explicit FFB no constraints' 'explicit FFB balanced' 'explicit FFB balanced + 2 consec max' 'implicit FFB'},'Location','southeast')
title(sprintf('Performance around reversal for CatAcc model \nalpha = %4.2f; beta = %d; %d ffb sequences; %d simulations per sequence',alphabeta1(1),alphabeta1(2),nseq,nsim))

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
end
