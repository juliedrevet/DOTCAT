function gen_ffb2(blck)
% generate acceptable ffb sequence

addpath ../

if nargin<1
    davg = 10;
    nepi = 12;
    ntrl = davg * nepi;
    dlim   = [5 ntrl];     % min|max number of trials before reversal
    blck = gen_epi(davg,dlim,nepi);
end


%% reversal curve from our subjects 1 2 3 6
perf = [0.0795    0.0682    0.0511    0.0341    0.5795    0.8068    0.8807    0.9091];

%% customizable option
alphabeta1 = [0.54 7];
alphabeta2 = [0.46 7];
alphabeta3 = [0.46 2];
alphabeta4 = [0.46 100];

%% get block configuration
ntrl       = length(blck.ys);
rfalse     = .2;%1-blck.p_reward;
stim_cor   = blck.ys; % no ffb
idxrev     = cumsum(blck.xs)+1; % reversal index
idxrev     = idxrev(1:end-1);
nrev       = length(idxrev);
win        = 4;

%% simulation parameters
mu     = abs(fzero(@(x)(normcdf(x)-rfalse),-1)); % x drawn from normal distribution mean +/-mu
nseq   = 1e3;
nsim   = 1e3;
r0     = 1 + (rand(nsim,1,nseq) < 0.5); % random 1st response
sftrnd = rand(nsim,ntrl+1,nseq);        % random softmax

%% generate FFB sequences
stim_ffb = repmat(blck.ys,nseq,1);
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

% evidence x accumulated throughout all ffb sequences and all simulations
x = zeros(nsim,ntrl,nseq);
for iseq = 1:nseq
    x(:,:,iseq) = repmat(sgnstim_ffb(iseq,:),nsim,1);
end

%% Simulate Models
% simulate CatAcc Model
[~,perfrev1sim,actn1sim,pcor1sim,z1sim,prep1sim,seqperf1sim] = modelsimulation(3,alphabeta1,perf);
[~,perfrev2sim,actn2sim,pcor2sim,z2sim,prep2sim,seqperf2sim] = modelsimulation(3,alphabeta2,perf);
[~,perfrev3sim,actn3sim,pcor3sim,z3sim,prep3sim,seqperf3sim] = modelsimulation(3,alphabeta3,perf);
[~,perfrev4sim,actn4sim,pcor4sim,z4sim,prep4sim,seqperf4sim] = modelsimulation(3,alphabeta4,perf);

%% Find ffb_seq giving a rev_curve next to the mean rev_curve
% compute zscore for all sequences and take the min zscore
zs1 = zscore(squeeze(mean(seqperf1sim,1)),0,2);
zs2 = zscore(squeeze(mean(seqperf2sim,1)),0,2);
zs3 = zscore(squeeze(mean(seqperf3sim,1)),0,2);
zs4 = zscore(squeeze(mean(seqperf4sim,1)),0,2);

[~,J1] = min(mean(zs1.^2));
[~,J2] = min(mean(zs2.^2));
[~,J3] = min(mean(zs3.^2));
[~,J4] = min(mean(zs4.^2));


c = lines;
figure()
plot(perfrev1sim,'LineWidth',2,'Color',c(1,:))
hold on
plot(mean(seqperf1sim(:,:,J1)),'--','Color',c(1,:))
plot(perfrev2sim,'LineWidth',2,'Color',c(2,:))
plot(mean(seqperf2sim(:,:,J2)),'--','Color',c(2,:))
plot(perfrev3sim,'LineWidth',2,'Color',c(3,:))
plot(mean(seqperf3sim(:,:,J3)),'--','Color',c(3,:))
plot(perfrev4sim,'LineWidth',2,'Color',c(4,:))
plot(mean(seqperf4sim(:,:,J4)),'--','Color',c(4,:))
title('Generate false feedback sequence based on performance curve at reversal')
ylabel('p(corr)')
xticklabels({'-4','','-2','','','2','','4'})
xlabel('trial position from reversal')
set(gca,'FontSize',20);
line([4.5 4.5],[0 1],'LineStyle','--','Color','k')
text(4.5,0.9,'reversal','HorizontalAlignment','center','FontAngle','italic','FontSize',16);
legend({sprintf('mean curve alpha = %4.2f beta = %d',alphabeta1(1),alphabeta1(2)) ...
sprintf('best FFB seq curve alpha = %4.2f beta = %d',alphabeta1(1),alphabeta1(2)) ...
sprintf('mean curve alpha = %4.2f beta = %d',alphabeta2(1),alphabeta2(2)) ...
sprintf('best FFB seq curve alpha = %4.2f beta = %d',alphabeta2(1),alphabeta2(2))...
sprintf('mean curve alpha = %4.2f beta = %d',alphabeta3(1),alphabeta3(2)) ...
sprintf('best FFB seq curve alpha = %4.2f beta = %d',alphabeta3(1),alphabeta3(2))...
sprintf('mean curve alpha = %4.2f beta = %d',alphabeta4(1),alphabeta4(2)) ...
sprintf('best FFB seq curve alpha = %4.2f beta = %d',alphabeta4(1),alphabeta4(2))},'FontSize',16,'Location','southeast')

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
