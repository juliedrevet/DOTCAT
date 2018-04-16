function models_play8
% function comparing fitted reversal curves when changing the false
% positive rate WITH EVIDENCE X THROWN RANDN

%% TO FIT: perf reversal curve for subjects 1 2 3 6 or p_cor
perfrev = [0.0795    0.0682    0.0511    0.0341    0.5795    0.8068    0.8807    0.9091];
pcor    =  0.9005;

%% generate block configuration as used during experiment
addpath ../
cfg = struct();
cfg.taskid = 1; % doesn't matter
cfg.condtn = 1; % unicolor dots
cfg.epimap = 1; % doesn't matter

blck = gen_blck(cfg);
ntrl     = blck.ntrl;
stim_cor = blck.ys; % no ffb
idxrev   = find(blck.switch_seq); % reversal index
nrev     = length(idxrev);
win      = 4;

c = lines;
l = 1;
legcell = cell(1);

for rfalse = [.2 .19 .18 .17 .16 .15]
    %% simulation parameters
    nseq = 1e3;
    nsim = 1e2;
    r0     = 1 + (rand(nsim,1,nseq) < 0.5); % random 1st response
    sftrnd = rand(nsim,ntrl+1,nseq);        % random softmax
    mu = abs(fzero(@(x)(normcdf(x)-rfalse),-1)); % x drawn from normal distribution mean +/-mu
    
    % our ffb sequences
    stim_ffb    = repmat(stim_cor,nseq,1); %take seq without ffb and add nseq various ffb sequences
    sgnstim_ffb = sign(stim_ffb-1.5);
    
    % evidence x accumulated throughout all ffb sequences and all simulations
    x = zeros(nsim,ntrl,nseq);
    for iseq = 1:nseq
        x(:,:,iseq) = repmat(sgnstim_ffb(iseq,:),nsim,1)*mu + randn(nsim,ntrl);
    end
    color1 = [255,170,128]/255;
    color2 = [128,191,255]/255;
    xidx = zeros(1,ntrl,nseq);
    xidx(1,:,:) = (sgnstim_ffb ==  1)';
    
    figure(1)
    hold on
    histogram(x(:),'BinWidth',.1,'FaceColor',[.8 .8 .8]);
    histogram(x(logical(repmat(xidx(1,:,:),nsim,1,1))),'BinWidth',.1,'FaceColor',color1);
    histogram(x(logical(repmat(1-xidx(1,:,:),nsim,1,1))),'BinWidth',.1,'FaceColor',color2);
    set(gca,'TickDir','out','XTick',[-mu,mu],'XTickLabel',{'-mu','+mu'},'FontSize',20);
    set(gca,'YTick',[],'YTickLabel',{''},'YColor',[1 1 1]);
    set(gca,'Layer','top','Box','off')
    
    
    %% find alpha and beta for ContAcc Model that fits subj bicolor reversal curve
    Eq = @(alphabeta)(modelsimulation2minimize(2,alphabeta,'perfrev'));
    alphabeta = fminsearch(Eq,[.5,5]);
    [~,perfrev2,~,pcor2,~,prep2] = modelsimulation(2,alphabeta);
    fprintf('Model ContAcc (%d%% false positives):\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',rfalse*100,alphabeta(1),alphabeta(2),100*pcor2,100*prep2)
    
    figure(2)
    if l == 1
        plot(perfrev,'LineWidth',2)
    end
    hold on
    plot(perfrev2,'--','Color',c(l+1,:))
    legcell{l} = sprintf('%d%%',rfalse*100);
    l = l+1;
end
    title('Fitting Performance around reversal with ContAcc Model','FontSize',30)
    ylabel('p(corr)')
    xticklabels({'-4','','-2','','','2','','4'})
    xlabel('trial position from reversal')
    set(gca,'FontSize',20);
    line([4.5 4.5],[0 1],'LineStyle','--','Color','k')
    text(4.5,0.9,'reversal','HorizontalAlignment','center','FontAngle','italic','FontSize',16);
    
    hleg = legend({'Subjects 1,2,3,6' legcell{:}},'FontSize',16,'Location','southeast');
    title(hleg, 'False positives rate');
    
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
