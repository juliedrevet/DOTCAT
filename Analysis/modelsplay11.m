function s = modelsplay11
% verify that if we fix "arbitrarily" alpha and beta to select the
% "best" FFB sequence, the reversal curve found for other alpha beta and
% this selected ffb sequence is not completely absurd
% NOTE: the FFB sequence constraint is BALANCED only (1st half / 2nd half)
% LONG RUNNING TIME (~40 min)!

tic

addpath ../

davg = 10;
nepi = 12;
ntrl = davg * nepi;
dlim = [5 ntrl];     % min|max number of trials before reversal
blck = gen_epi(davg,dlim,nepi);

%% fix alpha and beta to fit and to test
alphabeta0 = [0.54 7];
alphabetatest1 = [0.46 7];
alphabetatest2 = [0.54 2];
alphabetatest3 = [0.46 2];

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
[~,perfrevimpl0,~,~,~,~,~] = modelsimulation(3,alphabeta0,0);
[~,perfrevimpl1,~,~,~,~,~] = modelsimulation(3,alphabetatest1,0);
[~,perfrevimpl2,~,~,~,~,~] = modelsimulation(3,alphabetatest2,0);
[~,perfrevimpl3,~,~,~,~,~] = modelsimulation(3,alphabetatest3,0);
% 
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
[~,perfrev0,~,~,~,~,seqperf0sim] = modelsimulation(3,alphabeta0,0);
[~,perfrev1,~,~,~,~,seqperf1sim] = modelsimulation(3,alphabetatest1,0);
[~,perfrev2,~,~,~,~,seqperf2sim] = modelsimulation(3,alphabetatest2,0);
[~,perfrev3,~,~,~,~,seqperf3sim] = modelsimulation(3,alphabetatest3,0);

toc

zs0 = zscore(squeeze(mean(seqperf0sim,1)),0,2);
[~,J0] = min(mean(zs0.^2));
zs1 = zscore(squeeze(mean(seqperf1sim,1)),0,2);
[~,J1] = min(mean(zs1.^2));
zs2 = zscore(squeeze(mean(seqperf2sim,1)),0,2);
[~,J2] = min(mean(zs2.^2));
zs3 = zscore(squeeze(mean(seqperf3sim,1)),0,2);
[~,J3] = min(mean(zs3.^2));

%% create struct with results to return
s = struct();
s(1).alphabeta    = alphabeta0;
s(1).impl_perfrev = perfrevimpl0;
s(1).expl_perfrev = perfrev0;
s(1).all_perfrev  = seqperf0sim;
s(2).alphabeta    = alphabetatest1;
s(2).impl_perfrev = perfrevimpl1;
s(2).expl_perfrev = perfrev1;
s(2).all_perfrev  = seqperf1sim;
s(3).alphabeta    = alphabetatest2;
s(3).impl_perfrev = perfrevimpl2;
s(3).expl_perfrev = perfrev2;
s(3).all_perfrev  = seqperf2sim;
s(4).alphabeta    = alphabetatest3;
s(4).impl_perfrev = perfrevimpl3;
s(4).expl_perfrev = perfrev3;
s(4).all_perfrev  = seqperf3sim;

%% plot
c = lines;
figure()
subplot(2,2,1)
hold on
plot(perfrev0,'LineWidth',3,'Color',c(1,:))
plot(mean(seqperf0sim(:,:,J0)),'--','Color',c(1,:),'LineWidth',2)
plot(perfrevimpl0,'o','MarkerSize',10,'Color',c(1,:))
set(gca,'FontSize',16);
ylabel('p(corr)')
xlabel('trial position from reversal')
xticks(1:2*win+1)
xticklabels({'-4','','-2','','','2','','4'})
line([4.5 4.5],[0 1],'LineStyle','--','Color','k')
text(4.5,0.9,'reversal','HorizontalAlignment','center','FontAngle','italic','FontSize',20);
title('TRAINING : alpha = 0.54 beta = 7','FontSize',20);
hleg = legend({'mean','closest to mean (selected)','implicit'},'Location','southeast');
title(hleg,'False feedback sequences');


subplot(2,2,2)
hold on
plot(perfrev1,'LineWidth',3,'Color',c(2,:))
plot(mean(seqperf1sim(:,:,J0)),'--','Color',c(2,:),'LineWidth',2)
plot(mean(seqperf1sim(:,:,J1)),'LineWidth',.5,'Color',c(2,:))
plot(perfrevimpl1,'o','MarkerSize',10,'Color',c(2,:))
set(gca,'FontSize',16);
ylabel('p(corr)')
xlabel('trial position from reversal')
xticks(1:2*win+1)
xticklabels({'-4','','-2','','','2','','4'})
line([4.5 4.5],[0 1],'LineStyle','--','Color','k')
text(4.5,0.9,'reversal','HorizontalAlignment','center','FontAngle','italic','FontSize',20);
title('TESTING : alpha = 0.46 beta = 7','FontSize',20);
hleg = legend({'mean','selected in training','actually closest to mean','implicit'},'Location','southeast');
title(hleg,'False feedback sequences');

subplot(2,2,3)
hold on
plot(perfrev2,'LineWidth',3,'Color',c(3,:))
plot(mean(seqperf2sim(:,:,J0)),'--','Color',c(3,:),'LineWidth',2)
plot(mean(seqperf2sim(:,:,J2)),'LineWidth',.5,'Color',c(3,:))
plot(perfrevimpl2,'o','MarkerSize',10,'Color',c(3,:))
set(gca,'FontSize',16);
ylabel('p(corr)')
xlabel('trial position from reversal')
xticks(1:2*win+1)
xticklabels({'-4','','-2','','','2','','4'})
line([4.5 4.5],[0 1],'LineStyle','--','Color','k')
text(4.5,0.9,'reversal','HorizontalAlignment','center','FontAngle','italic','FontSize',20);
title('TESTING : alpha = 0.46 beta = 2','FontSize',20);
hleg = legend({'mean','selected in training','actually closest to mean','implicit'},'Location','southeast');
title(hleg,'False feedback sequences');


subplot(2,2,4)
hold on
plot(perfrev3,'LineWidth',3,'Color',c(4,:))
plot(mean(seqperf3sim(:,:,J0)),'--','Color',c(4,:),'LineWidth',2)
plot(mean(seqperf3sim(:,:,J3)),'LineWidth',.5,'Color',c(4,:))
plot(perfrevimpl3,'o','MarkerSize',10,'Color',c(4,:))
set(gca,'FontSize',16);
ylabel('p(corr)')
xlabel('trial position from reversal')
xticks(1:2*win+1)
xticklabels({'-4','','-2','','','2','','4'})
line([4.5 4.5],[0 1],'LineStyle','--','Color','k')
text(4.5,0.9,'reversal','HorizontalAlignment','center','FontAngle','italic','FontSize',20);
title('TESTING : alpha = 0.54 beta = 2','FontSize',20);
hleg = legend({'mean','selected in training','actually closest to mean','implicit'},'Location','southeast');
title(hleg,'False feedback sequences');




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
