function [out] = gen_draws_julie(davg,dlim,nepi)
%  GEN_DRAWS  Generate random draws
%
%  Usage: [out] = GEN_DRAWS(davg,dlim,nepi)
%
%  where davg - average episode length
%        dlim - minimum/maximum episode length
%        nepi - number of episodes to be generated
%
%  Valentin Wyart <valentin.wyart@ens.fr>

%%

if nargin < 3
    error('Missing input arguments!');
end

nper = 5e3; %5e3; % number of random permutations

plotornot = false; % plot results or not?

dmin = dlim(1); % minimum distance between reversals
dmax = dlim(2); % maximum distance between reversals

if davg < get_davg(1-eps) || davg > get_davg(eps)
    % eps, with no arguments, is the distance from 1.0 to the next larger double
    % precision number, that is eps with no arguments returns 2^(-52).
    error('Requested distance between reversals out of bounds!');
end

% get p(reversal)
prev = fzero(@(p)get_davg(p)-davg,0.5);

% plot p(reversal) wrt distance between reversals
if plotornot
    pvec = 0.01:0.01:0.50;
    dvec = nan(size(pvec));
    for i = 1:length(pvec)
        dvec(i) = get_davg(pvec(i));
    end
    figure('Color','white');
    hold on
    xlim([0,0.5]);
    plot(pvec,dvec,'k-','LineWidth',2);
    ylim(ylim);
    plot(prev*[1,1],ylim,'r-');
    plot(prev,davg,'ko','MarkerSize',12,'MarkerFaceColor',[1,0.5,0.5]);
    hold off
    set(gca,'Layer','top','Box','off','TickDir','out');
    set(gca,'FontName','Helvetica','FontSize',16);
    set(gca,'XTick',0:0.1:0.5);
    set(gca,'YTick',0:2:20);
    xlabel('p(reversal)');
    ylabel('distance between reversals');
end

% get theoretical pdf/cdf %loi géométrique "bornée" à dmin:dmax
% probabilité p(k) correspond à la probabilité d'obtenir dans une succession 
% de k épreuves de Bernoulli, k-1 échecs suivis d'un succès. 
% Les épreuves étant indépendantes, cette probabilité est de p*(1-p)^(k-1)
ps = prev*(1-prev).^(1:dmax); 
ps(1:dmin-1) = 0; % on vire les trop petits episodes
ps = ps/sum(ps); % normalisation
cs = cumsum(ps); % cdf to sample episode lengths

% sample distances between reversals
xs = nan(nper,nepi);
pr = nan(nper,1);
for iper = 1:nper
    % sample distances between reversals %from cs
    for iepi = 1:nepi
        xs(iper,iepi) = find(rand < cs,1,'first');
    end
    % fit p(reversal) wrt sampled distances between reversals
    fm = @(p,x)-sum(log(p)+log(1-p)*x)+nepi*log(sum(p*(1-p).^(dmin:dmax))); %2nd terme = constante de normalisation correspond à sum(ps) d'au-dessus
    pr(iper) = fminbnd(@(p)fm(p,xs(iper,:)),0.001,0.999);
end
ns = sum(xs,2); % nb total de trials par iper dim (nper,1)

% sort wrt offset between desired and observed p(reversal)
[~,is] = sort(abs(log(pr)-log(prev)));
xs_all = xs(is,:);
ns_all = ns(is);
pr_all = pr(is);

% filter sequences with desired number of sequences
is = 1:find(ns_all ~= nepi*davg,1)-1; %comment on sait que c'est consécutif?
xs_all = xs_all(is,:);
ns_all = ns_all(is);
pr_all = pr_all(is);

nseq = size(xs_all,1);

found = false;
for iseq = 1:nseq
    xs = xs_all(iseq,:);
    ns = ns_all(iseq);
    pr = pr_all(iseq);
    for iper = 1:nper
        xs = xs(randperm(nepi)); % repermute
        ys = unpack(xs);
        if ...  %difference entre nombre de 1 et nombre de 2 dans la 1ere et la 2eme moitié doit être <dmin >>? pourquoi dmin?
                (abs(diff(hist(ys(1:ns/2),[1,2]))) < round(davg/2) && ... %dmin round(davg/2)
                abs(diff(hist(ys(ns/2+1:end),[1,2]))) < round(davg/2))... %dmin
                || (length(xs)==2)
            if mod(nepi,2) > 0 %nombre impair d'episodes >> ? pourquoi c'est ok? du coup on a une seq = 6 7 8 16 13 ca m'a l'air bof (ptet de tte facon jamais de seq impaires)
                found = true;
                break
            elseif abs(sum(xs(1:nepi/2))-ns/2) < round(davg/3) %dmin>> ne pas avoir longueur 1ere moitié >>> 2eme moitié?
                found = true;
                break
            end
        end
    end
    if found
        break
    end
end
if ~found
    error('could not generate draws!');
end

% get biased p(reversal)
fm = @(p,x)-sum(log(p)+log(1-p)*x); % on vire la constante de normalisation
pb = fminbnd(@(p)fm(p,xs),0.001,0.999);

out      = [];
out.davg = davg;
out.dlim = dlim;
out.nepi = nepi;
out.prev = prev;
out.xs   = xs;
out.pr   = pr;
out.pb   = pb;
out.ys   = ys;

if plotornot
    % show sequence using imagesc
    % with pink and blue colors
end

    function [d] = get_davg(p) % average length of episode given p(reversal)
    ps2 = p*(1-p).^(dmin:dmax);
    ps2 = ps2/sum(ps2);
    d = sum(ps2.*(dmin:dmax));
    end

    function [y] = unpack(x) % déploie la séquence de 1 et 2 étant données les longueurs d'épisode
    n = sum(x);
    y = ones(1,n);
    for k = 1:length(x)-1
        j = sum(x(1:k))+1;
        y(j:end) = 3-y(j);
    end
    end

end