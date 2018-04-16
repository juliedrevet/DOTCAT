function models_play3
%
% modèle 1: "RL"
% z_t = z_{t-1} + alpha * (x - z_{t-1})
%
% modèle 2: "continuous accumulator"
% z_t = alpha * z_{t-1} + x
%
% modèle 3: "categorizing accumulator"
% z_t = alpha * z_{t-1} + sign(x)
%
% z_0 = 0

addpath ../
cfg = struct();
cfg.taskid = 1; % doesn't matter
cfg.condtn = 1; % unicolor dots
cfg.epimap = 1; % doesn't matter

%% fix simulation parameters
nsim   = 1e3;
mu     = abs(fzero(@(x)(normcdf(x)-.2),-1));

blck = gen_blck(cfg);
switch_seq = blck.switch_seq;
ntrl = blck.ntrl;
rfalse = .2;
stim_cor = blck.ys; % no ffb
sgnstim_cor = sign(stim_cor-1.5);
alpha3 = .65;
beta3  = 2;
nffb = 100;
allprev3 = zeros(nffb,8);
ipcor3 = zeros(1,nffb);
iprep3 = zeros(1,nffb);
idxrev = find(blck.switch_seq); % reversal index

r0     = 1 + (rand(nsim,1) < 0.5); % random 1st response
sftrnd = rand(nsim, ntrl+1);

for i = 1:nffb
    % generate ffb sequence
    false_seq = sample_false;
    idx_false = find(false_seq);
    % NOT: more than 2 consecutives false positives, false positive at the beginning
    % or at the end of the block
    while any(diff(idx_false(diff(idx_false)==1))==1)||ismember(1,idx_false)||ismember(ntrl,idx_false)... % maybe remove following conditions
            ||any(switch_seq(idx_false))||any(switch_seq(idx_false+1))
        false_seq = sample_false;
        idx_false = find(false_seq);
    end
    
    stim_ffb = blck.ys; %with ffb
    stim_ffb(false_seq == 1) = 3-stim_ffb(false_seq == 1);

    sgnstim1 = sign(stim_ffb-1.5);
    %% simulate Model 3 with ffb
    stim = stim_ffb;
    x    = mu*sgnstim1+randn(nsim,ntrl);
    actions3 = sim_model(3,alpha3,beta3);
    
    % compute Model 3 correct and reversal rates
    stim = stim_cor; % reversal curves with respect to correct (no ffb!)
    ipcor3(i) = mean(mean(bsxfun(@eq,actions3(:,2:end),stim),2),1);
    iprep3(i) = sum(sum(diff(actions3')==0))/numel(diff(actions3'));
    iprev3    = zeros(nsim,8);
    for isim = 1:nsim
        iprev3(isim,:) = rev_curve(actions3(isim,:),4,true);
    end
    allprev3(i,:) = mean(iprev3);

end
pcor3 = mean(ipcor3);
prev3 = mean(allprev3);
prep3 = mean(iprep3);

%% find alpha beta back with RL without the ffb
stim = stim_cor;
%x    = x-mu*sgnstim1+mu*sgnstim_cor;
x    = mu*sgnstim_cor+randn(nsim,ntrl);

%% find alpha and beta for Model RL that fits Model 3 reversal curve
Eq11 = @(alphabeta1)(mean((revcurve(1,alphabeta1(1),alphabeta1(2))-prev3).^2));
alphabeta1 = fminsearch(Eq11,[alpha3,beta3]);
prev1    = revcurve(1,alphabeta1(1),alphabeta1(2));
actions1 = sim_model(1,alphabeta1(1),alphabeta1(2));
pcor1    = mean(mean(bsxfun(@eq,actions1(:,2:end),stim),2),1);
prep1    = sum(sum(diff(actions1')==0))/numel(diff(actions1'));

fprintf('Model 3 - categorizing accumulator:\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alpha3,beta3,100*pcor3,100*prep3)
fprintf('Model 1 - RL (fit):\t\t\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alphabeta1(1),alphabeta1(2),100*pcor1,100*prep1)



figure()
plot(prev3,'LineWidth',2)
hold on
plot(prev1)


    function perf = rev_curve(resp,win,flp)
        nrev = length(idxrev);
        itokeep = zeros(nrev,2*win);
        iperf = zeros(nrev,2*win);
        
        k = 1;
        for irev = idxrev
            iwin = (irev-win):(irev+win-1);
            itokeep(k,:) = logical(iwin>0 & iwin<ntrl+1);
            iresp   = resp(iwin(logical(itokeep(k,:)))+1);
            if flp
                icorrect = stim(irev);
            else
                icorrect = stim(iwin);
            end
            iperf(k,logical(itokeep(k,:))) = (iresp == icorrect);
            k = k+1;
        end
        
        perf = sum(iperf)./sum(itokeep);
    end

    function rev = revcurve(model,alpha,beta)
        if model == 1
            actions = sim_model(1,alpha,beta);
        elseif model == 2
            actions = sim_model(2,alpha,beta);
        elseif model == 3
            actions = sim_model(3,alpha,beta);
        end
        prev = zeros(nsim,2*4);
        for s = 1:nsim
            prev(s,:) = rev_curve(actions(s,:),4,true);
        end
        rev = mean(prev);
    end

    function actn = sim_model(model,alpha,beta)
        z = zeros(nsim,ntrl+1);
        r = zeros(nsim,ntrl+1);
        r(:,1) = r0;
        for t = 1:ntrl
            if model == 1
                z(:,t+1) = z(:,t) + alpha * (x(:,t) - z(:,t));
            elseif model == 2
                z(:,t+1) = alpha * z(:,t) + x(:,t);
            elseif model == 3
                z(:,t+1) = alpha * z(:,t) + sign(x(:,t));
            end
            r(:,t+1) = 1 + (z(:,t+1) > 0);
        end
        actn = (sftrnd>(1./(1+exp(beta*z))))+1;
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

% Model 3 - categorizing accumulator:	 alpha: 0.65	 beta: 2.00 	 p_cor: 66.51%
% Model 1 - RL (fit):                    alpha: 0.37	 beta: 1.89 	 p_cor: 66.68%
% models_play3
% Model 3 - categorizing accumulator:	 alpha: 0.50	 beta: 1.00 	 p_cor: 61.33%
% Model 1 - RL (fit):                    alpha: 0.60	 beta: 0.67 	 p_cor: 61.08%
% models_play3
% Model 3 - categorizing accumulator:	 alpha: 0.50	 beta: 1.00 	 p_cor: 61.34%
% Model 1 - RL (fit):                    alpha: 0.57	 beta: 0.73 	 p_cor: 61.43%
% models_play3
% Model 3 - categorizing accumulator:	 alpha: 0.50	 beta: 2.00 	 p_cor: 66.20%
% Model 1 - RL (fit):                    alpha: 0.64	 beta: 1.15 	 p_cor: 66.74%
% models_play3
% Model 3 - categorizing accumulator:	 alpha: 0.50	 beta: 2.00 	 p_cor: 65.95%
% Model 1 - RL (fit):                    alpha: 0.65	 beta: 1.10 	 p_cor: 66.19%
% models_play3
% Model 3 - categorizing accumulator:	 alpha: 0.65	 beta: 2.00 	 p_cor: 65.84%
% Model 1 - RL (fit):                    alpha: 0.39	 beta: 1.77 	 p_cor: 66.38%
% models_play3
% Model 3 - categorizing accumulator:	 alpha: 0.65	 beta: 2.00 	 p_cor: 66.19% 	 p_rep: 81.83%
% Model 1 - RL (fit):                    alpha: 0.40	 beta: 1.69 	 p_cor: 66.61% 	 p_rep: 69.56%

% after changing the x
% models_play3
% Model 3 - categorizing accumulator:	 alpha: 0.65	 beta: 10.00 	 p_cor: 69.28% 	 p_rep: 67.69%
% Model 1 - RL (fit):                    alpha: 0.39	 beta: 1.54 	 p_cor: 68.80% 	 p_rep: 60.67%
% models_play3
% Model 3 - categorizing accumulator:	 alpha: 0.65	 beta: 2.00 	 p_cor: 66.70% 	 p_rep: 66.38%
% Model 1 - RL (fit):                    alpha: 0.41	 beta: 1.26 	 p_cor: 66.64% 	 p_rep: 58.16%






