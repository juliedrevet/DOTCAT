function simulmodel(expe,iblck)
% modèle simulation billes marbrées
% variables fixées: sigma = 1, mu de sorte à avoir 20% de cdf quand < 0
%
% modèle 1: sorte de RL
% z_t = z_{t-1} + alpha * (x - z_{t-1})
%
% modèle 2: "leaky accumulator"
% z_t = alpha * z_{t-1} + x
%
% modèle 3: categorizer
% z_t = alpha * z_{t-1} + sign(x)
%
% z_0 = 0

nsim = 1e3;
% sig = 1 >> find +mu
mu = abs(fzero(@(x)(normcdf(x)-.2),-1));

stim = expe.stim(iblck).outcome(1,:);
ntrl    = length(stim);
sgnstim = sign(stim-1.5);

alpha_vec = 0:0.01:1;

for i = 1:length(alpha_vec)
    alpha = alpha_vec(i);
    z1 = zeros(nsim,1);
    z2 = zeros(nsim,1);
    z3 = zeros(nsim,1);
    r1 = zeros(nsim,ntrl);
    r2 = zeros(nsim,ntrl);
    r3 = zeros(nsim,ntrl);
    
    for itrl = 2:ntrl
        x = mu*sgnstim(itrl)+randn(nsim,1);
        z1 = z1 + alpha * (x - z1);
        z2 = alpha * z2 + x;
        z3 = alpha * z3 + sign(x);
        r1(:,itrl) = 1 + (z1 > 0);
        r2(:,itrl) = 1 + (z2 > 0);
        r3(:,itrl) = 1 + (z3 > 0);
    end
    % get accuracy
    pcor1(i) = mean(mean(bsxfun(@eq,r1,stim),2),1);
    pcor2(i) = mean(mean(bsxfun(@eq,r2,stim),2),1);
    pcor3(i) = mean(mean(bsxfun(@eq,r3,stim),2),1);
end
alpha1 = .5;%fminbnd(@(l)-interp1(alpha_vec,pcor1,l,'spline'),0,1);
alpha2 = .5;%fminbnd(@(l)-interp1(alpha_vec,pcor2,l,'spline'),0,1)
alpha3 = .5;%fminbnd(@(l)-interp1(alpha_vec,pcor3,l,'spline'),0,1)

% re-simulate models with best alphas
z1 = zeros(nsim,1);
z2 = zeros(nsim,1);
z3 = zeros(nsim,1);
r1 = zeros(nsim,ntrl+1);
r2 = zeros(nsim,ntrl+1);
r3 = zeros(nsim,ntrl+1);
r1(:,1) = 1 + (rand(nsim,1) < 0.5); % random 1st response
r2(:,1) = 1 + (rand(nsim,1) < 0.5); % random 1st response
r3(:,1) = 1 + (rand(nsim,1) < 0.5); % random 1st response
r11 = r1;
r22 = r2;
r33 = r3;
z12(nsim,1) = 0;
z22(nsim,1) = 0;
z32(nsim,1) = 0;

for itrl = 1:ntrl
    x = mu*sgnstim(itrl)+randn(nsim,1);
    z1 = z1 + alpha1 * (x - z1);
    z2 = alpha2 * z2 + x;
    z3 = alpha3 * z3 + sign(x);
    r1(:,itrl+1) = 1 + (z1 > 0);
    r2(:,itrl+1) = 1 + (z2 > 0);
    r3(:,itrl+1) = 1 + (z3 > 0);
    % or follow time course
    z12(:,itrl+1) = z12(:,itrl) + alpha1 * (x - z12(:,itrl));
    z22(:,itrl+1) = alpha2 * z22(:,itrl) + x;
    z32(:,itrl+1) = alpha3 * z32(:,itrl) + sign(x);
    r11(:,itrl+1) = 1 + (z12(:,itrl+1) > 0);
    r22(:,itrl+1) = 1 + (z22(:,itrl+1) > 0);
    r33(:,itrl+1) = 1 + (z32(:,itrl+1) > 0);
end

pcor1 = mean(mean(bsxfun(@eq,r1(:,2:end),stim),2),1);
pcor2 = mean(mean(bsxfun(@eq,r2(:,2:end),stim),2),1);
pcor3 = mean(mean(bsxfun(@eq,r3(:,2:end),stim),2),1);

%% add softmax
beta = 1;
% same softmax rnd for all models
sftrnd = rand(size(z32));
actions12 = (sftrnd>(1./(1+exp(beta*z12))))+1;
actions22 = (sftrnd>(1./(1+exp(beta*z22))))+1;
actions32 = (sftrnd>(1./(1+exp(beta*z32))))+1;

pcor12 = mean(mean(bsxfun(@eq,actions12(:,2:end),stim),2),1);
pcor22 = mean(mean(bsxfun(@eq,actions22(:,2:end),stim),2),1);
pcor32 = mean(mean(bsxfun(@eq,actions32(:,2:end),stim),2),1);

% prep1 = 
% prep2 = 
% prep3 = 

% prev1 = 
% prev2 =
for isim = 1:nsim
    prev1(isim,:) = rev_curve(stim,actions12(isim,:),4,true);
    prev2(isim,:) = rev_curve(stim,actions22(isim,:),4,true);
    prev3(isim,:) = rev_curve(stim,actions32(isim,:),4,true);
end


figure()
plot(stim)
hold on
plot(mean(r1(:,2:end)))
plot(mean(r2(:,2:end)))
plot(mean(r3(:,2:end)))
plot(mean(actions12(:,2:end)))
plot(mean(actions22(:,2:end)))
plot(mean(actions32(:,2:end)))
ylim([.9 2.1])
% title(sprintf('alpha_{RL} = %4.2f; alpha_{LA} = %4.2f\n pcor_{RL} = %02d, pcor_{LA} = %02d',alpha1,alpha2,round(pcor1*100),round(pcor2*100)))
% legend({'stim','RL','leaky accumulator'})

figure()
hold on
plot(mean(prev1))
plot(mean(prev2))
plot(mean(prev3))

end

function perf = rev_curve(stim,resp,win,flp)

nrev = length((find(diff(stim))+1));
itokeep = zeros(nrev,2*win);
iperf = zeros(nrev,2*win);
ntrl = length(stim);


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





