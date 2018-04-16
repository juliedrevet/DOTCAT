function gen_ffb3(blck)
% generate acceptable ffb sequence

addpath ../

if nargin<1
    davg = 10;
    nepi = 12;
    ntrl = davg * nepi;
    dlim   = [5 ntrl];     % min|max number of trials before reversal
    blck = gen_epi(davg,dlim,nepi);
end

% subj
% 0.0795    0.0682    0.0511    0.0341    0.5795    0.8068    0.8807    0.9091
% explicit ffb
% 0.1717    0.1512    0.1401    0.1280    0.5270    0.7196    0.7877    0.8305


%% customizable option
alphabeta1 = [0.54 7];
alphabeta2 = [0.46 7];
alphabeta3 = [0.46 2];
alphabeta4 = [0.46 100];

%% get block configuration
ntrl     = length(blck.ys);
rfalse   = .2;%1-blck.p_reward;
stim_cor = blck.ys; % no ffb
idxrev   = cumsum(blck.xs)+1; % reversal index
idxrev   = idxrev(1:end-1);
nrev     = length(idxrev);
win      = 4;
nfalse   = rfalse*ntrl;

%% simulation parameters
mu     = abs(fzero(@(x)(normcdf(x)-rfalse),-1)); % x drawn from normal distribution mean +/-mu
nseq   = 1e4;
nsim   = 5e2;
r0     = 1 + (rand(nsim,1,nseq) < 0.5); % random 1st response
sftrnd = rand(nsim,ntrl+1,nseq);        % random softmax
x      = zeros(nsim,ntrl,nseq);

%% generate model simulation with evidence randn (internal FFB)
stim_ffb    = repmat(stim_cor,nseq,1); %take seq without ffb and add nseq various ffb sequences
sgnstim_ffb = sign(stim_ffb-1.5);
% randn evidence x accumulated throughout all ffb sequences and all simulations
for iseq = 1:nseq
    x(:,:,iseq) = repmat(sgnstim_ffb(iseq,:),nsim,1)*mu + randn(nsim,ntrl);
end

tic
%% compute reversal curve from simulated CatAcc Model
[~,perfrev0sim,~,~,~,~,~] = modelsimulation(3,alphabeta1,0);
toc

%% generate CatAcc Model simulation with evidence dizcretized (external FFB balanced + 2 consecutive max)
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
tic
% compute reversal curve
[~,perfrev1sim,~,~,~,~,seqperf1sim] = modelsimulation(3,alphabeta1,0);
toc

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
tic
% compute reversal curve
[~,perfrev2sim,~,~,~,~,seqperf2sim] = modelsimulation(3,alphabeta1,0);
toc

%% Find ffb_seq giving a rev_curve next to the mean rev_curve
% compute zscore for all sequences and take the min zscore
zs1 = zscore(squeeze(mean(seqperf1sim,1)),0,2);
zs2 = zscore(squeeze(mean(seqperf2sim,1)),0,2);

[~,J1] = min(mean(zs1.^2));
[~,J2] = min(mean(zs2.^2));

c = lines;
figure()
hold on
plot(perfrev1sim,'LineWidth',2,'Color',c(1,:))
plot(mean(seqperf1sim(:,:,J1)),'--','Color',c(1,:))
plot(perfrev2sim,'LineWidth',2,'Color',c(2,:))
plot(mean(seqperf2sim(:,:,J2)),'--','Color',c(2,:))
plot(perfrev0sim,'k.','MarkerSize',20)
title('Generate false feedback sequence based on performance curve around reversal')
ylabel('p(corr)')
xticklabels({'-4','','-2','','','2','','4'})
xlabel('trial position from reversal')
set(gca,'FontSize',20);
line([4.5 4.5],[0 1],'LineStyle','--','Color','k')
text(4.5,0.9,'reversal','HorizontalAlignment','center','FontAngle','italic','FontSize',16);

hleg = legend({'Balanced + max 2 consecutives' 'best ffb sequence' 'Balanced only' 'best ffb sequence' 'implicit'},'Location','southeast');
title(hleg, 'False feedback sequences');
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

function [err] = modelsimulation2minimize(model,alphabeta,p,err2minimize)
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
        
        % compute error depending on what is to minimize
        if strcmp(err2minimize,'perfrev')
            iperf = zeros(nrev,2*win);
            
            k = 1;
            for irev = idxrev
                iwin = (irev-win+1):(irev+win);
                iresp   = actn(:,iwin,:);
                icorrect = stim_cor(irev);
                iperf(k,:) = mean(mean(iresp==icorrect),3);
                k = k+1;
            end
            mperf = mean(iperf,1);
            err = mean((mperf-p).^2);
            
        elseif strcmp(err2minimize,'pcor')
            err = abs(mean(mean(mean(bsxfun(@eq,actn(:,2:end,:),stim_cor),2),1),3)-p);
        elseif strcmp(err2minimize,'prep')
            iprep = zeros(1,nseq);
            for k = 1:nseq
                iprep(k) = sum(sum(diff(actn(:,:,k),1,2)==0))/numel(diff(actn(:,:,k),1,2));
            end
            prep = mean(iprep);
            err = abs(prep-p);
        end
    end
end
