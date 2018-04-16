function models_play
%% play with models to have a better insight
% cannot be runned as is, onoly in mode debug
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
nblck = 10;

pcor3 = zeros(nblck,11);
prev3 = zeros(nblck,8,11);
pcor2 = zeros(nblck,11);
prev2 = zeros(nblck,8,11);

for iblck = 1:nblck
    %% generate 10 outcome sequences
    blck(iblck) = gen_blck(cfg);
    outc_ffb(iblck,:) = blck.colorffb_seq;
    outc_ffb(iblck,outc_ffb(iblck,:) == 0) = 2;
    outc(iblck,:) = blck(iblck).color_seq;
    outc(iblck,outc(iblck,:) == 0) = 2;
    sgnoutc_ffb(iblck,:) = sign(outc_ffb(iblck,:)-1.5);
    sgnoutc(iblck,:) = sign(outc(iblck,:)-1.5);
    
    stim = outc(iblck,:);
    sgnstim = sign(stim-1.5);
    ntrl   = length(stim);
    x      = mu*sgnstim+randn(nsim,1);
    r0     = 1 + (rand(nsim,1) < 0.5); % random 1st response
    sftrnd = rand(nsim, ntrl+1);
    
    kalpha = 1;
    for alpha3 = .5% .25:.05:.75
        kbeta = 1;
        for beta3 = .5:.1:1.5 %1
            %% simulate Model 3 (categorizing accumulator)
            actions3 = sim_model(3,alpha3,beta3);
            actions2 = sim_model(2,alpha3,beta3);
            %% correct responses
            pcor3(iblck,kbeta) = mean(mean(bsxfun(@eq,actions3(:,2:end),stim),2),1);
            pcor2(iblck,kbeta) = mean(mean(bsxfun(@eq,actions2(:,2:end),stim),2),1);
            %% compute reversal curve for Model 3
            iprev3 = zeros(nsim,8);
            iprev2 = zeros(nsim,8);
            for isim = 1:nsim
                iprev3(isim,:) = rev_curve(actions3(isim,:),4,true);
                iprev2(isim,:) = rev_curve(actions2(isim,:),4,true);
            end
            prev3(iblck,:,kbeta) = mean(iprev3);
            prev2(iblck,:,kbeta) = mean(iprev2);
            kbeta = kbeta+1;
        end
        kalpha = kalpha+1;
    end
    
end

figure()
subplot(2,1,1)
plot(mean(prev3(:,:,1)))
hold on
plot(mean(prev3(:,:,2)))
plot(mean(prev3(:,:,3)))
plot(mean(prev3(:,:,4)))
plot(mean(prev3(:,:,5)))
plot(mean(prev3(:,:,6)))
plot(mean(prev3(:,:,7)))
plot(mean(prev3(:,:,8)))
plot(mean(prev3(:,:,9)))
plot(mean(prev3(:,:,10)))
plot(mean(prev3(:,:,11)))
%%
subplot(2,1,2)
plot(mean(prev2(:,:,1)))
hold on
plot(mean(prev2(:,:,2)))
plot(mean(prev2(:,:,3)))
plot(mean(prev2(:,:,4)))
plot(mean(prev2(:,:,5)))
plot(mean(prev2(:,:,6)))
plot(mean(prev2(:,:,7)))
plot(mean(prev2(:,:,8)))
plot(mean(prev2(:,:,9)))
plot(mean(prev2(:,:,10)))
plot(mean(prev2(:,:,11)))
%%
legend({'0.5'    '0.6'    '0.7'    '0.8'    '0.9'    '1.0'    '1.1'    '1.2'    '1.3'    '1.4'    '1.5'})
title('alpha = 0.5, changing beta')
%%
legend({'0.25'    '0.30'    '0.35'    '0.40'    '0.45'    '0.50'    '0.55'    '0.60'    '0.65'    '0.70'    '0.75'});
title('changing alpha, beta = 1')


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







