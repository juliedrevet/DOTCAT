function simulmodel3(expe,iblck)
%% this function simulates responses from categorizing accumulator model
%  (model 3) based on bicolor outcome (no ffb) with fixed alpha and beta
%  parameters and fit continuous accumulator and RL model only based on the
%  reversal curve resulting from model 3 simulations.
%
% modèle simulation billes marbrées
% variables fixées: sigma = 1, mu de sorte à avoir 20% de cdf quand < 0
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

blck = gen_blck(cfg);
stim_ffb = blck.colorffb_seq;
stim_ffb(stim_ffb == 0) = 2;
stim = blck.color_seq;
stim(stim == 0) = 2;

ntrl    = length(stim);
sgnstim = sign(stim-1.5);

%% fix simulation parameters
nsim   = 1e4;
mu     = abs(fzero(@(x)(normcdf(x)-.2),-1));
alpha3 = .5;
beta3  = 1;
x      = mu*sgnstim+randn(nsim,1);
r0     = 1 + (rand(nsim,1) < 0.5); % random 1st response
sftrnd = rand(nsim, ntrl+1);

%% simulate Model 3 (categorizing accumulator)
actions3 = sim_model(3,alpha3,beta3);
pcor3    = mean(mean(bsxfun(@eq,actions3(:,2:end),stim),2),1);

%% compute reversal curve for Model 3
for isim = 1:nsim
    iprev3(isim,:) = rev_curve(actions3(isim,:),4,true);
end
prev3 = mean(iprev3);

%% find alpha and beta for Model 1 and 2 that fits Model 3 reversal curve
Eq1=@(alphabeta1)(mean((revcurve(1,alphabeta1(1),alphabeta1(2))-prev3).^2));
Eq2=@(alphabeta2)(mean((revcurve(2,alphabeta2(1),alphabeta2(2))-prev3).^2));

alphabeta1 = fminsearch(Eq1,[alpha3,beta3]);
alphabeta2 = fminsearch(Eq2,[alpha3,beta3]);
prev1    = revcurve(1,alphabeta1(1),alphabeta1(2));
actions1 = sim_model(1,alphabeta1(1),alphabeta1(2));
pcor1    = mean(mean(bsxfun(@eq,actions1(:,2:end),stim),2),1);

prev2    = revcurve(2,alphabeta2(1),alphabeta2(2));
actions2 = sim_model(2,alphabeta2(1),alphabeta2(2));
pcor2    = mean(mean(bsxfun(@eq,actions2(:,2:end),stim),2),1);

prep1 = sum(sum(diff(actions1')==0))/numel(diff(actions1'));
prep2 = sum(sum(diff(actions2')==0))/numel(diff(actions2'));
prep3 = sum(sum(diff(actions3')==0))/numel(diff(actions3'));

fprintf('Model 3 - categorizing accumulator:\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alpha3,beta3,100*pcor3,100*prep3)
fprintf('Model 1 - RL (fit):\t\t\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alphabeta1(1),alphabeta1(2),100*pcor1,100*prep1)
fprintf('Model 2 - continuous accumulator (fit):\t alpha: %4.2f\t beta: %4.2f \t p_cor: %5.2f%% \t p_rep: %5.2f%%\n',alphabeta2(1),alphabeta2(2),100*pcor2,100*prep2)


figure()
plot(prev3,'LineWidth',2);
hold on
plot(prev1);
plot(prev2);
legend({'Categorizer','RL','accumulator'})



    function perf = rev_curve(resp,win,flp)
        nrev = length((find(diff(stim))+1));
        itokeep = zeros(nrev,2*win);
        iperf = zeros(nrev,2*win);
        
        k = 1;
        for irev = (find(diff(stim))+1)
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







