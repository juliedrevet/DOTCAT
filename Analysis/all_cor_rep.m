function [p_cor,p_corafter,p_rep] = all_cor_rep(clem)

if nargin<1
    clem = false;
end

%% load data
expes = load_data(clem);
if clem
    figname = './fig/DOTCAT_pilot_clem_';
else
    figname = './fig/DOTCAT_pilot_';
end

keepsubj = [1 2 3 6];
expes = expes(keepsubj);

p_cor = struct();
p_corafter = struct();
p_rep = struct();

p = 1;
for expe = expes
    
    %% run analysis
    p_cor(p).ob1 = 0;
    p_cor(p).ob2 = 0;
    p_cor(p).ag1 = 0;
    p_cor(p).ag2 = 0;
    p_corafter(p).ob1 = 0;
    p_corafter(p).ob2 = 0;
    p_corafter(p).ag1 = 0;
    p_corafter(p).ag2 = 0;
    p_rep(p).ob1 = 0;
    p_rep(p).ob2 = 0;
    p_rep(p).ag1 = 0;
    p_rep(p).ag2 = 0;
    
    j1 = 0;
    j2 = 0;
    j3 = 0;
    j4 = 0;
    
    for i = (expe(1).cfg.nprac+1):length(expe.cfg.taskid)
        if expe.blck(i).taskid == 1
            if expe.blck(i).condtn == 1
                p_cor(p).ob1 = p_cor(p).ob1+expe.rslt(i).perf;
                if isfield(expe.rslt(i),'perfafter')
                    p_corafter(p).ob1 = p_corafter(p).ob1+expe.rslt(i).perfafter;
                else
                    p_corafter(p).ob1 = p_corafter(p).ob1+...
                        sum(expe.rslt(i).resp(2:end) == expe.rslt(i).correct(1:(end-1)))/length(expe.rslt(i).resp(2:end));
                end
                p_rep(p).ob1 = p_rep(p).ob1+sum(diff(expe.rslt(i).resp)==0)/length(diff(expe.rslt(i).resp));
                j1 = j1+1;
            else
                p_cor(p).ob2 = p_cor(p).ob2+expe.rslt(i).perf;
                if isfield(expe.rslt(i),'perfafter')
                    p_corafter(p).ob2 = p_corafter(p).ob2+expe.rslt(i).perfafter;
                else
                    p_corafter(p).ob2 = p_corafter(p).ob2+...
                        sum(expe.rslt(i).resp(2:end) == expe.rslt(i).correct(1:(end-1)))/length(expe.rslt(i).resp(2:end));
                end
                p_rep(p).ob2 = p_rep(p).ob2+sum(diff(expe.rslt(i).resp)==0)/length(diff(expe.rslt(i).resp));
                j2 = j2+1;
            end
        else
            if expe.blck(i).condtn == 1
                p_cor(p).ag1 = p_cor(p).ag1+expe.rslt(i).perf;
                if isfield(expe.rslt(i),'perfafter')
                    p_corafter(p).ag1 = p_corafter(p).ag1+expe.rslt(i).perfafter;
                else
                    p_corafter(p).ag1 = p_corafter(p).ag1+...
                        sum(expe.rslt(i).resp(2:end) == expe.rslt(i).correct(1:(end-1)))/length(expe.rslt(i).resp(2:end));
                end
                p_rep(p).ag1 = p_rep(p).ag1+sum(diff(expe.rslt(i).resp)==0)/length(diff(expe.rslt(i).resp));
                j3 = j3+1;
            else
                p_cor(p).ag2 = p_cor(p).ag2+expe.rslt(i).perf;
                if isfield(expe.rslt(i),'perfafter')
                    p_corafter(p).ag2 = p_corafter(p).ag2+expe.rslt(i).perfafter;
                else
                    p_corafter(p).ag2 = p_corafter(p).ag2+...
                        sum(expe.rslt(i).resp(2:end) == expe.rslt(i).correct(1:(end-1)))/length(expe.rslt(i).resp(2:end));
                end
                p_rep(p).ag2 = p_rep(p).ag2+sum(diff(expe.rslt(i).resp)==0)/length(diff(expe.rslt(i).resp));
                j4 = j4+1;
            end
        end
        
    end
    
    p_cor(p).ob1 = p_cor(p).ob1/j1;
    p_cor(p).ob2 = p_cor(p).ob2/j2;
    p_cor(p).ag1 = p_cor(p).ag1/j3;
    p_cor(p).ag2 = p_cor(p).ag2/j4;
    p_corafter(p).ob1 = p_corafter(p).ob1/j1;
    p_corafter(p).ob2 = p_corafter(p).ob2/j2;
    p_corafter(p).ag1 = p_corafter(p).ag1/j3;
    p_corafter(p).ag2 = p_corafter(p).ag2/j4;
    p_rep(p).ob1 = p_rep(p).ob1/j1;
    p_rep(p).ob2 = p_rep(p).ob2/j2;
    p_rep(p).ag1 = p_rep(p).ag1/j3;
    p_rep(p).ag2 = p_rep(p).ag2/j4;
    
    p = p+1;
end

c = lines;

figure()
bar(1,mean([p_cor.ob1]),'FaceColor',c(1,:),'LineWidth',1); hold on
bar(2,mean([p_cor.ag1]),'FaceColor',c(7,:),'LineWidth',1); hold on
bar(3,mean([p_cor.ob2]),'FaceColor',c(1,:),'LineWidth',1); hold on
bar(4,mean([p_cor.ag2]),'FaceColor',c(7,:),'LineWidth',1); hold on
xticks([])
ylim([0 1])
ylabel('p(corr)')
savefig([figname,'p_cor'])

figure()
bar(1,mean([p_corafter.ob1]),'FaceColor',c(1,:),'LineWidth',1); hold on
bar(2,mean([p_corafter.ag1]),'FaceColor',c(7,:),'LineWidth',1); hold on
bar(3,mean([p_corafter.ob2]),'FaceColor',c(1,:),'LineWidth',1); hold on
bar(4,mean([p_corafter.ag2]),'FaceColor',c(7,:),'LineWidth',1); hold on
xticks([])
ylim([0 1])
ylabel('p(corr_{after})')
savefig([figname,'p_corafter'])

figure()
bar(1,mean([p_rep.ob1]),'FaceColor',c(1,:),'LineWidth',1); hold on
bar(2,mean([p_rep.ag1]),'FaceColor',c(7,:),'LineWidth',1); hold on
bar(3,mean([p_rep.ob2]),'FaceColor',c(1,:),'LineWidth',1); hold on
bar(4,mean([p_rep.ag2]),'FaceColor',c(7,:),'LineWidth',1); hold on
xticks([])
ylim([0 1])
ylabel('p_{rep}')
savefig([figname,'p_rep'])


end
