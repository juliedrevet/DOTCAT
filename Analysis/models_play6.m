function models_play6
% with fixed alpha & beta, examinate the reversal curve with ffb sequences
% thrown without constraints

%% reversal curve from our subjects 1 2 3 6
perf = [0.0795    0.0682    0.0511    0.0341    0.5795    0.8068    0.8807    0.9091];

%% customizable options
alphabeta = [0.55 200];
add_ffb   = false; % ffb or no ffb at all
ffb_ctrl  = false; % ffb sequence generation controlled (constraints)
x_randn   = true; % x thrown from 2 normal distributions around +/- mu
plotornot = true;

%% simulation parameters
mu     = abs(fzero(@(x)(normcdf(x)-.2),-1)); % x drawn from normal distribution mean +/-mu
nseq   = 1e2;
nsim   = 1e3;
rfalse = .2;

%% generate 1 block configuration (as used during experiment)
addpath ../
cfg = struct();
cfg.taskid = 1; % doesn't matter
cfg.condtn = 1; % unicolor dots
cfg.epimap = 1; % doesn't matter

blck = gen_blck(cfg);
while (sum(blck.ys==2)~=.5*blck.ntrl)
    blck = gen_blck(cfg);
end
switch_seq = blck.switch_seq;
ntrl = blck.ntrl;
stim_cor = blck.ys; % no ffb
idxrev = find(blck.switch_seq); % reversal index
nrev = length(idxrev);
win = 4;
r0     = 1 + (rand(nsim,1,nseq) < 0.5); % random 1st response
sftrnd = rand(nsim,ntrl+1,nseq);        % random softmax

stim_ffb = repmat(blck.ys,nseq,1);
% our ffb sequences
for iseq = 1:nseq
    if ffb_ctrl && add_ffb
        false_seq = sample_false;
        idx_false = find(false_seq);
        % NOT: more than 2 consecutives false positives, false positive at the beginning
        % or at the end of the block
        while any(diff(idx_false(diff(idx_false)==1))==1)||ismember(1,idx_false)||ismember(ntrl,idx_false)... % maybe remove following conditions
                ||any(switch_seq(idx_false))||any(switch_seq(idx_false+1))
            false_seq = sample_false;
            idx_false = find(false_seq);
        end
    else
        false_seq = zeros(1,ntrl);
        if add_ffb
            false_seq(1:24) = 1; % comment if you want no ffb at all
            false_seq = Shuffle(false_seq); % comment if you want no ffb at all
        end
    end

    % inject ffb in stim sequence
    stim_ffb(iseq,false_seq == 1) = 3-stim_ffb(iseq,false_seq == 1);
end
sgnstim_ffb = sign(stim_ffb-1.5);

% evidence x accumulated throughout all ffb sequences and all simulations
x = zeros(nsim,ntrl,nseq);
for iseq = 1:nseq
    if x_randn
        x(:,:,iseq) = repmat(sgnstim_ffb(iseq,:),nsim,1)*mu + randn(nsim,ntrl);
    else
        x(:,:,iseq) = repmat(sgnstim_ffb(iseq,:),nsim,1);
    end
end


if plotornot
    color1 = [255,170,128]/255;
    color2 = [128,191,255]/255;
    xidx = zeros(1,ntrl,nseq);
    xidx(1,:,:) = (sgnstim_ffb ==  1)';
    figure()
    hold on
    histogram(x(:),'BinWidth',.1,'FaceColor',[.8 .8 .8]);
    h1 = histogram(x(logical(repmat(xidx(1,:,:),nsim,1,1))),'BinWidth',.1,'FaceColor',color1);
    h2 = histogram(x(logical(repmat(1-xidx(1,:,:),nsim,1,1))),'BinWidth',.1,'FaceColor',color2);
    set(gca,'TickDir','out','XTick',[-mu,mu],'XTickLabel',{'-mu','+mu'},'FontSize',20);
    set(gca,'YTick',[],'YTickLabel',{''},'YColor',[1 1 1]);
    set(gca,'Layer','top','Box','off')
    line([mu mu],[0 max(h1.Values)],'LineStyle','--','LineWidth',2,'Color','k')
    line([-mu -mu],[0 max(h2.Values)],'LineStyle','--','LineWidth',2,'Color','k')
    
    figure()
    hold on
    h1 = histogram(x(logical(repmat(xidx(1,:,:),nsim,1,1))),'BinWidth',.1,'FaceColor',color1);
    h2 = histogram(x(logical(repmat(1-xidx(1,:,:),nsim,1,1))),'BinWidth',.1,'FaceColor',color2);
    set(gca,'TickDir','out','XTick',[-mu,mu],'XTickLabel',{'-mu','+mu'},'FontSize',20);
    set(gca,'YTick',[],'YTickLabel',{''},'YColor',[1 1 1]);
    set(gca,'Layer','top','Box','off')
    line([mu mu],[0 max(h1.Values)],'LineStyle','--','LineWidth',2,'Color','k')
    line([-mu -mu],[0 max(h2.Values)],'LineStyle','--','LineWidth',2,'Color','k')
end

%% simulate CatAcc Model
[~,perfrev1,actn1,pcor1,z1] = modelsimulation(3,alphabeta);

%% simulate ContAcc Model
[~,perfrev2,actn2,pcor2,z2] = modelsimulation(2,alphabeta);
%% simulate RL Continous Model
[~,perfrev3,actn3,pcor3,z3] = modelsimulation(1,alphabeta);

figure()
plot(perf,'LineWidth',2)
hold on
plot(perfrev1)
plot(perfrev2,'*-')
plot(perfrev3,'--')
title('Performance at reversal')
ylabel('p(corr)')
xticklabels({'-4','','-2','','','2','','4'})
xlabel('trial position from reversal')
line([4.5 4.5],[0 1],'LineStyle','--','Color','k')
text(4.5,0.9,'reversal','HorizontalAlignment','center','FontAngle','italic','FontSize',16);
legend({'Subjects 1,2,3,6','CatAcc Model','ContAcc Model','ContRL Model'},'FontSize',16,'Location','southeast')

figure()
plot(stim_cor,'LineWidth',2);
hold on
plot(mean(mean(actn1,1),3))
plot(mean(mean(actn2,1),3))
plot(mean(mean(actn3,1),3))
title('mean actions taken');
legend({'correct responses','CatAcc Model actions','ContAcc Model actions','ContRL Model actions'},'FontSize',14,'Location','southeast')

figure()
hold on
plot(mean(mean(z1,1),3))
plot(mean(mean(z2,1),3))
plot(mean(mean(z3,1),3))
title('mean accumulated evidence');

figure()
pbar = 1;
hold on;
bar(1,pcor1,'LineWidth',1);
bar(2,pcor2,'LineWidth',1);
bar(3,pcor3,'LineWidth',1);
%errorbar(1:2, [mean(ma(idx1yn,1)); mean(ma(idx2yn,1))],[std(ma(idx1yn,1))/sqrt(length(ma(idx1yn,1))); std(ma(idx2yn,1))/sqrt(length(ma(idx2yn,1)))],'k.'); hold on 
ylim([0 1]);
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]); 
set(gca,'TickDir','out');%,'TickLength',[1,1]*0.02/max(pbar,1)); 
set(gca,'FontName','Helvetica','FontSize',16); 
set(gca,'YTick',0:.1:1,'FontSize',12);
set(gca,'XTick',[1:3]);
set(gca,'XTickLabel',{'CatAcc','ContAcc','ContRL'},'FontSize',18,'FontWeight','bold');
ylabel('p(corr)','FontSize',18,'FontWeight','bold')

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


    function [eperf,mperf,actn,pcor,z,prep] = modelsimulation(model,alphabeta)
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
        
        k = 1;
        for irev = idxrev
            iwin = (irev-win+1):(irev+win);
            iresp   = actn(:,iwin,:);
            icorrect = stim_cor(irev);
            iperf(k,:) = mean(mean(iresp==icorrect),3);
            k = k+1;
        end
        
         mperf = mean(iperf,1);
         eperf = mean((mperf-perf).^2);
    end
end
