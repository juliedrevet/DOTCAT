function models_play5
% function trying to find alpha and beta parameter of our 3 Models when
% fitted to the reversal curve of bicolor condition of subjects 1 2 3 6
% 1. it generates nseq different false feedback sequences and nsim
% simulations of the model, 2. computes the reversal curves for our 3
% models and compared to subject reversal curve (least squares) to find 
% the best fitting alpha and beta.

%% TO FIT: perf reversal curve for subjects 1 2 3 6 or p_cor
% bicolor
perfrev = [0.0795    0.0682    0.0511    0.0341    0.5795    0.8068    0.8807    0.9091];
% unicolor
%perfrev = [0.3125    0.1989    0.1534    0.0568    0.1193    0.3920    0.6307    0.7386];
pcor =  0.9005;

%% customizable options
add_ffb   = true; % ffb or no ffb at all
ffb_ctrl  = true; % ffb sequence generation controlled (constraints)
x_randn   = false; % x thrown from 2 normal distributions around +/- mu
plotornot = true;

%% simulation parameters
nseq = 1e3;
nsim = 1e2;

%% generate block configuration as used during experiment
addpath ../
cfg = struct();
cfg.taskid = 1; % doesn't matter
cfg.condtn = 1; % unicolor dots
cfg.epimap = 1; % doesn't matter

blck = gen_blck(cfg);
switch_seq = blck.switch_seq;
ntrl   = blck.ntrl;
rfalse = .2;
stim_cor = blck.ys; % no ffb
idxrev = find(blck.switch_seq); % reversal index
nrev   = length(idxrev);
win    = 4;
r0     = 1 + (rand(nsim,1,nseq) < 0.5); % random 1st response
sftrnd = rand(nsim,ntrl+1,nseq);        % random softmax
mu = abs(fzero(@(x)(normcdf(x)-rfalse),-1)); % x drawn from normal distribution mean +/-mu



% our ffb sequences
stim_ffb = repmat(stim_cor,nseq,1); %take seq without ffb and add nseq various ffb sequences
for iseq = 1:nseq
    if ffb_ctrl && add_ffb
        false_seq = sample_false;
        idx_false = find(false_seq);
        % NOT: more than 2 consecutives false positives, false positive at the beginning
        % or at the end of the block
        while any(diff(idx_false(diff(idx_false)==1))==1)||ismember(1,idx_false)||ismember(ntrl,idx_false)%... % maybe remove following conditions
                %||any(switch_seq(idx_false))||any(switch_seq(idx_false+1))
            false_seq = sample_false;
            idx_false = find(false_seq);
        end
    else
        false_seq = zeros(1,ntrl);
        if add_ffb
            false_seq(1:24) = 1; 
            false_seq = Shuffle(false_seq);
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
    histogram(x(logical(repmat(xidx(1,:,:),nsim,1,1))),'BinWidth',.1,'FaceColor',color1);
    histogram(x(logical(repmat(1-xidx(1,:,:),nsim,1,1))),'BinWidth',.1,'FaceColor',color2);
    set(gca,'TickDir','out','XTick',[-mu,mu],'XTickLabel',{'-mu','+mu'},'FontSize',20);
    set(gca,'YTick',[],'YTickLabel',{''},'YColor',[1 1 1]);
    set(gca,'Layer','top','Box','off')
end

%% find alpha and beta for CatAcc Model that fits subj bicolor reversal curve
Eq1 = @(alphabeta1)(modelsimulation2minimize(3,alphabeta1,'perfrev'));
alphabeta1 = fminsearch(Eq1,[.5,5]);
[eperf1,perfrev1,actn1,pcor1,z1,prep1] = modelsimulation(3,alphabeta1);
fprintf('Model CatAcc:\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alphabeta1(1),alphabeta1(2),100*pcor1,100*prep1)

%% find alpha and beta for ContAcc Model that fits subj bicolor reversal curve
Eq2 = @(alphabeta2)(modelsimulation2minimize(2,alphabeta2,'perfrev'));
alphabeta2 = fminsearch(Eq2,[.5,5]);
[eperf2,perfrev2,actn2,pcor2,z2,prep2] = modelsimulation(2,alphabeta2);
fprintf('Model ContAcc:\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alphabeta2(1),alphabeta2(2),100*pcor2,100*prep2)

%% find alpha and beta for ContRL Model that fits subj bicolor reversal curve
Eq3 = @(alphabeta3)(modelsimulation2minimize(1,alphabeta3,'perfrev'));
alphabeta3 = fminsearch(Eq3,[.5,5]);
[eperf3,perfrev3,actn3,pcor3,z3,prep3] = modelsimulation(1,alphabeta3);
fprintf('Model ContRL:\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alphabeta3(1),alphabeta3(2),100*pcor3,100*prep3)


figure()
plot(perfrev,'LineWidth',2)
hold on
plot(perfrev1)
plot(perfrev2,'*-')
plot(perfrev3,'--')
title('Performance around reversal')
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
legend({'correct responses','CatAcc Model','ContAcc Model','ContRL Model'},'FontSize',14,'Location','southeast')

figure()
hold on
plot(mean(mean(z1,1),3))
plot(mean(mean(z2,1),3))
plot(mean(mean(z3,1),3))
title('mean accumulated evidence');
legend({'CatAcc Model','ContAcc Model','ContRL Model'},'FontSize',14,'Location','southeast')


pbar = 1;
figure()
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
         eperf = mean((mperf-perfrev).^2);
    end
    function [err] = modelsimulation2minimize(model,alphabeta,err2minimize)
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
            err = mean((mperf-perfrev).^2);
            
        elseif strcmp(err2minimize,'pcor')
            err = abs(mean(mean(mean(bsxfun(@eq,actn(:,2:end,:),stim_cor),2),1),3)-pcor);
        elseif strcmp(err2minimize,'prep')
            iprep = zeros(1,nseq);
            for k = 1:nseq
                iprep(k) = sum(sum(diff(actn(:,:,k),1,2)==0))/numel(diff(actn(:,:,k),1,2));
            end
            prep = mean(iprep);
            err = abs(prep-pcor);
        end
    end
end
