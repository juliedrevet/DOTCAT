function [alphabeta,taskid,condtn] = findab(subj, fun)
% estimate alpha and beta per block for subject subj

if nargin<2
    fun = @fminsearch;
end

% load data
subj_data = load_data(~ismember(1:23,subj));
expe = subj_data.expe;

taskid = expe(1).cfg.taskid;
condtn = expe(1).cfg.condtn;
alphabeta = zeros(length(expe),2);
alphabeta2 = zeros(length(expe),2);


for b = 1:length(expe)
    
    choice = expe(b).rslt.resp;
    
    if expe(b).cfg.taskid(b) == 1 % not dependent on choice
        outcome = expe(b).blck.outcome(1,:);
    elseif expe(b).cfg.taskid(b) == 2 % dependent on choice
            outcome = expe(b).blck.epimap*expe(b).blck.color_seq+(3-expe(b).blck.epimap)*(~expe(b).blck.color_seq);
    end
    
    Eq=@(alpha,beta)(neglog(deltaQlearn(alpha, outcome), choice,beta));
    
    LLH=@(X)Eq(X(1),X(2));
    
    if isequal(fun,@fminsearch)
        alphabeta(b,:)  = fminsearch(LLH,[0.01,5]);
    elseif isequal(fun,@fmincon)
        alphabeta(b,:) = fmincon(LLH,[0.01,5],[],[],[],[],[0.001 0],[1 50],[], ...
            optimset('Display','notify','FunValCheck','on','Algorithm','interior-point','TolX',1e-20,'MaxFunEvals',1e6));
    else
        error('not a valid function')
    end
end
end