function [p_cor,p_corafter,p_rep] = cor_rep(rmb)

if nargin<1
    rmb = 0;
end
[FileName,PathName] = uigetfile('*.mat','Select subject experiment data');
expe = importdata(fullfile(PathName,FileName));


p_cor.obsdot  = 0;
p_cor.obshand = 0;
p_cor.agdot   = 0;
p_cor.aghand  = 0;
p_corafter.obsdot  = 0;
p_corafter.obshand = 0;
p_corafter.agdot   = 0;
p_corafter.aghand  = 0;
p_rep.obsdot  = 0;
p_rep.obshand = 0;
p_rep.agdot   = 0;
p_rep.aghand  = 0;

j1 = 0;
j2 = 0;
j3 = 0;
j4 = 0;

for i = (expe(1).cfg.nprac+1):length(expe)
    if i~=rmb
        if expe(i).blck.taskid == 1
            if expe(i).blck.condtn == 1
                p_cor.obsdot = p_cor.obsdot+expe(i).rslt.perf;
                if isfield(expe(i).rslt,'perfafter')
                    p_corafter.obsdot = p_corafter.obsdot+expe(i).rslt.perfafter;
                else
                    p_corafter.obsdot = p_corafter.obsdot+...
                        sum(expe(i).rslt.resp(2:end) == expe(i).rslt.correct(1:(end-1)))/length(expe(i).rslt.resp(2:end));
                end
                p_rep.obsdot = p_rep.obsdot+sum(diff(expe(i).rslt.resp)==0)/length(diff(expe(i).rslt.resp));
                j1 = j1+1;
            else
                p_cor.obshand = p_cor.obshand+expe(i).rslt.perf;
                if isfield(expe(i).rslt,'perfafter')
                    p_corafter.obshand = p_corafter.obshand+expe(i).rslt.perfafter;
                else
                    p_corafter.obshand = p_corafter.obshand+...
                        sum(expe(i).rslt.resp(2:end) == expe(i).rslt.correct(1:(end-1)))/length(expe(i).rslt.resp(2:end));
                end
                p_rep.obshand = p_rep.obshand+sum(diff(expe(i).rslt.resp)==0)/length(diff(expe(i).rslt.resp));
                j2 = j2+1;
            end
        else
            if expe(i).blck.condtn == 1
                p_cor.agdot = p_cor.agdot+expe(i).rslt.perf;
                if isfield(expe(i).rslt,'perfafter')
                    p_corafter.agdot = p_corafter.agdot+expe(i).rslt.perfafter;
                else
                    p_corafter.agdot = p_corafter.agdot+...
                        sum(expe(i).rslt.resp(2:end) == expe(i).rslt.correct(1:(end-1)))/length(expe(i).rslt.resp(2:end));
                end
                p_rep.agdot = p_rep.agdot+sum(diff(expe(i).rslt.resp)==0)/length(diff(expe(i).rslt.resp));
                j3 = j3+1;
            else
                p_cor.aghand = p_cor.aghand+expe(i).rslt.perf;
                if isfield(expe(i).rslt,'perfafter')
                    p_corafter.aghand = p_corafter.aghand+expe(i).rslt.perfafter;
                else
                    p_corafter.aghand = p_corafter.aghand+...
                        sum(expe(i).rslt.resp(2:end) == expe(i).rslt.correct(1:(end-1)))/length(expe(i).rslt.resp(2:end));
                end
                p_rep.aghand = p_rep.aghand+sum(diff(expe(i).rslt.resp)==0)/length(diff(expe(i).rslt.resp));
                j4 = j4+1;
            end
        end
    end
end

p_cor.obsdot  = p_cor.obsdot/j1;
p_cor.obshand = p_cor.obshand/j2;
p_cor.agdot   = p_cor.agdot/j3;
p_cor.aghand  = p_cor.aghand/j4;
p_corafter.obsdot  = p_corafter.obsdot/j1;
p_corafter.obshand = p_corafter.obshand/j2;
p_corafter.agdot   = p_corafter.agdot/j3;
p_corafter.aghand  = p_corafter.aghand/j4;
p_rep.obsdot  = p_rep.obsdot/j1;
p_rep.obshand = p_rep.obshand/j2;
p_rep.agdot   = p_rep.agdot/j3;
p_rep.aghand  = p_rep.aghand/j4;

end
