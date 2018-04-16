function models_play2
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
% with ffb
stim1 = blck.colorffb_seq;
stim1(stim1 == 0) = 2;
% without ffb
stim2 = blck.color_seq;
stim2(stim2 == 0) = 2;

% reversal index
idxrev = find(blck.switch_seq);

sgnstim1 = sign(stim1-1.5);
sgnstim2 = sign(stim2-1.5);

ntrl   = length(stim1);
r0     = 1 + (rand(nsim,1) < 0.5); % random 1st response
sftrnd = rand(nsim, ntrl+1);

alpha3 = .65;
beta3  = 2;

%% simulate Model 3 with ffb
stim = stim1;
x    = mu*sgnstim1+randn(nsim,1);
actions3 = sim_model(3,alpha3,beta3);

% compute Model 3 correct and reversal rates
stim = stim2; % reversal curves with respect to correct (no ffb!)
pcor3 = mean(mean(bsxfun(@eq,actions3(:,2:end),stim),2),1);
iprev = zeros(nsim,8);
for isim = 1:nsim
    iprev3(isim,:) = rev_curve(actions3(isim,:),4,true);
end
prev3 = mean(iprev3);

%% find alpha beta back with RL without the ffb
stim = stim2;
%x    = x-mu*sgnstim1+mu*sgnstim2;
x    = mu*sgnstim2+randn(nsim,1);

%% find alpha and beta for Model RL that fits Model 3 reversal curve
Eq11=@(alphabeta1)(mean((revcurve(1,alphabeta1(1),alphabeta1(2))-prev3).^2));
alphabeta1 = fminsearch(Eq11,[alpha3,beta3]);
prev11    = revcurve(1,alphabeta1(1),alphabeta1(2));
actions11 = sim_model(1,alphabeta1(1),alphabeta1(2));
pcor11    = mean(mean(bsxfun(@eq,actions11(:,2:end),stim),2),1);




%% find alpha and beta for Model RL that fits Model 3 reversal curve
stim = stim1;
x    = x-mu*sgnstim2+mu*sgnstim1;
Eq12=@(alphabeta2)(mean((revcurve(1,alphabeta2(1),alphabeta2(2))-prev3).^2));
alphabeta2 = fminsearch(Eq12,[alpha3,beta3]);
prev12    = revcurve(1,alphabeta2(1),alphabeta2(2));
actions12 = sim_model(1,alphabeta2(1),alphabeta2(2));
pcor12    = mean(mean(bsxfun(@eq,actions12(:,2:end),stim),2),1);


figure()
plot(prev3,'LineWidth',2)
hold on
plot(prev11)
plot(prev12)




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

end







