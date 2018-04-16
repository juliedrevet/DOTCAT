function [perf] = all_rev_perf(win,flp,clem)
% reponse itrl+1 par rapport à correct itrl (autre config donne meme résultat!)

if nargin<1
    win = 4;
    flp = true;
    clem = false;
elseif nargin==1
    flp = true;
    clem = false;
elseif nargin == 2
    clem = false;
end

if clem
    figname = './fig/DOTCAT_pilot_clem_';
else
    figname = './fig/DOTCAT_pilot_bis_';
end

%% load data
expes = load_data(clem);
keepsubj = [1 2 3 6];
%keepsubj = [6];
expes = expes(keepsubj);
nsubj = length(expes);


perf = struct();
p = 1;
for expe = expes
    %% get configuration parameters
    nprac = expe.cfg.nprac;
    nblck = length(expe.cfg.taskid);
    ntrl  = expe.blck(nprac+1).ntrl;
    nrev = sum(expe.blck(nprac+1).switch_seq);
    %% run analysis
    perfob1 = [];
    tokeepob1 = [];
    perfob2 = [];
    tokeepob2 = [];
    perfag1 = [];
    tokeepag1 = [];
    perfag2 = [];
    tokeepag2 = [];
    itokeep = zeros(nrev,2*win,nblck-nprac);
    iperf = zeros(nrev,2*win,nblck-nprac);
    for iblck = nprac+1:nblck
        rslt = expe.rslt(iblck);
        k = 1;
        for irev = find(expe.blck(iblck).switch_seq)
            iwin = (irev-win):(irev+win-1);
            itokeep(k,:,iblck-nprac) = logical(iwin>0 & iwin<ntrl+1);
            iresp   = rslt.resp(iwin(logical(itokeep(k,:,iblck-nprac)))+1);
            if flp
                icorrect = rslt.correct(irev);
            else
                icorrect = rslt.correct(iwin);
            end
            iperf(k,logical(itokeep(k,:,iblck-nprac)),iblck-nprac) = (iresp == icorrect);
            k = k+1;
        end
        if (expe.cfg.taskid(iblck) == 1) && (expe.cfg.condtn(iblck) == 1)
            perfob1 = [perfob1; iperf(:,:,iblck-nprac)];
            tokeepob1 = [tokeepob1; itokeep(:,:,iblck-nprac)];
        elseif (expe.cfg.taskid(iblck) == 1) && (expe.cfg.condtn(iblck) == 2)
            perfob2 = [perfob2; iperf(:,:,iblck-nprac)];
            tokeepob2 = [tokeepob2; itokeep(:,:,iblck-nprac)];
        elseif (expe.cfg.taskid(iblck) == 2) && (expe.cfg.condtn(iblck) == 1)
            perfag1 = [perfag1; iperf(:,:,iblck-nprac)];
            tokeepag1 = [tokeepag1; itokeep(:,:,iblck-nprac)];
        elseif (expe.cfg.taskid(iblck) == 2) && (expe.cfg.condtn(iblck) == 2)
            perfag2 = [perfag2; iperf(:,:,iblck-nprac)];
            tokeepag2 = [tokeepag2; itokeep(:,:,iblck-nprac)];
        end
    end
    
    perf(p).ob1 = sum(perfob1)./sum(tokeepob1);
    perf(p).ob2 = sum(perfob2)./sum(tokeepob2);
    perf(p).ag1 = sum(perfag1)./sum(tokeepag1);
    perf(p).ag2 = sum(perfag2)./sum(tokeepag2);
    
    p = p+1;
end

% reshape perf
perfob1 = reshape([perf.ob1],2*win,nsubj)';
perfob2 = reshape([perf.ob2],2*win,nsubj)';
perfag1 = reshape([perf.ag1],2*win,nsubj)';
perfag2 = reshape([perf.ag2],2*win,nsubj)';

%% plot
c = lines;
figure() % OBS/AG
hold on
plot((-win+.5:win),mean([perfob1;perfob2]),'Linewidth',2,'Color',c(1,:))
plot((-win+.5:win),mean([perfag1;perfag2]),'Linewidth',2,'Color',c(2,:))
xlabel('sequence position from reversal')
ylabel('p(correct)')
line([-win win],[.5 .5],'LineStyle',':','Color','k')
line([0 0],[0 .95],'LineStyle','--','Color','k')
xlim([-win win])
box off
legend('observer','agent','Location','southeast')
legend('boxoff')
text(0,0.97,'reversal','HorizontalAlignment','center','VerticalAlignment','baseline','FontAngle','italic','FontSize',18);
savefig([figname,'rev_obag'])

figure() % Uni/Bicolor
hold on
plot((-win+.5:win),mean([perfob1;perfag1]),'Linewidth',2,'Color',c(3,:))
plot((-win+.5:win),mean([perfob2;perfag2]),'Linewidth',2,'Color',c(4,:))
xlabel('sequence position from reversal')
ylabel('p(correct)')
line([-win win],[.5 .5],'LineStyle',':','Color','k')
line([0 0],[0 .95],'LineStyle','--','Color','k')
xlim([-win win])
box off
legend('unicolor dots','bicolor dots','Location','southeast')
legend('boxoff')
text(0,0.97,'reversal','HorizontalAlignment','center','VerticalAlignment','baseline','FontAngle','italic','FontSize',18);
savefig([figname,'rev_12'])

figure() % all
hold on
plot((-win+.5:win),mean(perfob1),'-','Linewidth',2,'Color',c(1,:))
plot((-win+.5:win),mean(perfag1),'-','Linewidth',2,'Color',c(2,:))
plot((-win+.5:win),mean(perfob2),'--','Linewidth',2,'Color',c(1,:))
plot((-win+.5:win),mean(perfag2),'--','Linewidth',2,'Color',c(2,:))
xlabel('sequence position from reversal')
ylabel('p(correct)')
line([-win win],[.5 .5],'LineStyle',':','Color','k')
line([0 0],[0 .95],'LineStyle','--','Color','k')
xlim([-win win])
box off
legend('OBS unicolor','AGENT unicolor','OBS bicolor','AGENT bicolor','Location','southeast')
legend('boxoff')
text(0,0.97,'reversal','HorizontalAlignment','center','VerticalAlignment','baseline','FontAngle','italic','FontSize',18);
savefig([figname,'rev_all'])

end