function models_play4
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

%% TO FIT: perf reversal curve for subjects 1 2 3 6
perf = [0.0795    0.0682    0.0511    0.0341    0.5795    0.8068    0.8807    0.9091];

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
idxrev = find(blck.switch_seq); % reversal index

r0     = 1 + (rand(nsim,1) < 0.5); % random 1st response
sftrnd = rand(nsim, ntrl+1);

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

% false_seq = zeros(1,ntrl);
% false_seq(1:24) = 1;
% false_seq = Shuffle(false_seq);

stim_ffb = blck.ys; %with ffb
stim_ffb(false_seq == 1) = 3-stim_ffb(false_seq == 1);

sgnstim1 = sign(stim_ffb-1.5);

stim = stim_ffb;
x    = mu*sgnstim1+randn(nsim,ntrl);


%% find alpha and beta for CatAcc Model that fits subj bicolor reversal curve
Eq1 = @(alphabeta1)(fct2minimize(stim_cor,3,alphabeta1));
alphabeta1 = fminsearch(Eq1,[.5,5]);
prev1      = revcurve(stim_cor,3,alphabeta1(1),alphabeta1(2));
actions1   = sim_model(3,alphabeta1(1),alphabeta1(2));
pcor1      = mean(mean(bsxfun(@eq,actions1(:,2:end),stim_cor),2),1);

%% find alpha and beta for ContAcc Model that fits subj bicolor reversal curve
Eq2 = @(alphabeta2)(fct2minimize(stim_cor,2,alphabeta2));
alphabeta2 = fminsearch(Eq2,[.5,5]);
prev2      = revcurve(stim_cor,2,alphabeta2(1),alphabeta2(2));
actions2   = sim_model(2,alphabeta2(1),alphabeta2(2));
pcor2      = mean(mean(bsxfun(@eq,actions2(:,2:end),stim_cor),2),1);

figure()
plot(perf,'LineWidth',2)
hold on
plot(prev1)
plot(prev2)


    function perf = rev_curve(stim,resp,win,flp)
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

    function rev = revcurve(stim,model,alpha,beta)
        if model == 1
            actions = sim_model(1,alpha,beta);
        elseif model == 2
            actions = sim_model(2,alpha,beta);
        elseif model == 3
            actions = sim_model(3,alpha,beta);
        end
        prev = zeros(nsim,2*4);
        for s = 1:nsim
            prev(s,:) = rev_curve(stim,actions(s,:),4,true);
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

    function eperf = fct2minimize(stim,model,alphabeta)
        alpha = alphabeta(1);
        beta  = alphabeta(2);
%         actions3 = sim_model(model,alpha,beta);
%         % compute Model 3 reversal curve
%         iprev3    = zeros(nsim,8);
%         for isim = 1:nsim
%             iprev3(isim,:) = rev_curve(stim,actions3(isim,:),4,true);
%         end
%         allprev3 = mean(iprev3);
%         eperf = mean((allprev3-perf).^2);
        eperf = mean((revcurve(stim,model,alpha,beta)-perf).^2);
    end
end
