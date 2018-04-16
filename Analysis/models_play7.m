function models_play7
% simulate and find back parameters alpha and beta

%% reversal curve from our subjects 1 2 3 6
perf = [0.0795    0.0682    0.0511    0.0341    0.5795    0.8068    0.8807    0.9091];

%% customizable options
alphabeta = [0.3 1];
add_ffb   = true; % ffb or no ffb at all
ffb_ctrl  = false; % ffb sequence generation controlled (constraints)
x_randn   = true; % x thrown from 2 normal distributions around +/- mu
plotornot = false;

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

%% Simulate Models
% simulate CatAcc Model
[~,perfrev1sim,actn1sim,pcor1sim,z1sim,prep1sim] = modelsimulation(3,alphabeta,perf);
% simulate ContAcc Model
[~,perfrev2sim,actn2sim,pcor2sim,z2sim,prep2sim] = modelsimulation(2,alphabeta,perf);
% simulate RL Continous Model
[~,perfrev3sim,actn3sim,pcor3sim,z3sim,prep3sim] = modelsimulation(1,alphabeta,perf);

%% Find back parameters
% CatAcc Model
Eq1 = @(alphabeta1)(modelsimulation2minimize(3,alphabeta1,perfrev1sim,'perfrev'));
alphabeta1 = fminsearch(Eq1,[.1,5]);
[~,perfrev1fit,actn1fit,pcor1fit,z1fit,prep1fit] = modelsimulation(3,alphabeta1,perf);
fprintf('Model CatAcc (sim):\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alphabeta(1),alphabeta(2),100*pcor1sim,100*prep1sim)
fprintf('Model CatAcc (fit):\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alphabeta1(1),alphabeta1(2),100*pcor1fit,100*prep1fit)
% ContAcc Model
Eq2 = @(alphabeta2)(modelsimulation2minimize(2,alphabeta2,perfrev2sim,'perfrev'));
alphabeta2 = fminsearch(Eq2,[.1,5]);
[~,perfrev2fit,actn2fit,pcor2fit,z2fit,prep2fit] = modelsimulation(2,alphabeta2,perf);
fprintf('Model ContAcc (sim):\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alphabeta(1),alphabeta(2),100*pcor2sim,100*prep2sim)
fprintf('Model ContAcc (fit):\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alphabeta2(1),alphabeta2(2),100*pcor2fit,100*prep2fit)
% ContRL Model that fits subj bicolor reversal curve
Eq3 = @(alphabeta3)(modelsimulation2minimize(1,alphabeta3,perfrev3sim,'perfrev'));
alphabeta3 = fminsearch(Eq3,[.1,5]);
[~,perfrev3fit,actn3fit,pcor3fit,z3fit,prep3fit] = modelsimulation(1,alphabeta3,perf);
fprintf('Model ContRL (sim):\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alphabeta(1),alphabeta(2),100*pcor3sim,100*prep3sim)
fprintf('Model ContRL (fit):\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alphabeta3(1),alphabeta3(2),100*pcor3fit,100*prep3fit)

c = lines;
figure()
hold on
plot(perfrev1sim,'LineWidth',2)
plot(perfrev2sim,'LineWidth',2)
plot(perfrev3sim,'LineWidth',2)
plot(perfrev1fit,'o','MarkerSize',10,'Color',c(1,:))
plot(perfrev2fit,'o','MarkerSize',10,'Color',c(2,:))
plot(perfrev3fit,'o','MarkerSize',10,'Color',c(3,:))
plot(perf,'LineWidth',3,'Color',[.5 .5 .5])
title(sprintf('Performance at reversal with alpha = %4.2f and beta =  %4.2f',alphabeta(1),alphabeta(2)))
ylabel('p(corr)')
xticklabels({'-4','','-2','','','2','','4'})
xlabel('trial position from reversal')
line([4.5 4.5],[0 1],'LineStyle','--','Color','k')
text(4.5,0.9,'reversal','HorizontalAlignment','center','FontAngle','italic','FontSize',16);
legend({'CatAcc Model','ContAcc Model','ContRL Model','CatAcc Model (fit)','ContAcc Model (fit)','ContRL Model (fit)','subjects 1 2 3 6'},'FontSize',16,'Location','southeast')




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


    function [eperf,mperf,actn,pcor,z,prep] = modelsimulation(model,alphabeta,perfrev)
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
