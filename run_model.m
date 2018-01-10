function perf = run_model(blck, nsim, ntrl, perc)

if blck.condtn == 1 
    % single dot outcome: include false feedbacks
    color_outcome = blck.epimap*blck.colorffb_seq+(3-blck.epimap)*(~blck.colorffb_seq);
    % run only one wsls model
    wsls_ans = zeros(1,blck.ntrl);
    wsls_ans(1) = randi(2);
    for i = 2:blck.ntrl
        if wsls_ans(i-1) == color_outcome(i-1)
            wsls_ans(i) = wsls_ans(i-1);
        else
            wsls_ans(i) = 3-wsls_ans(i-1);
        end
    end
elseif blck.condtn == 2
    % handful dots outcome: no false feedbacks
    color_outcome = blck.epimap*blck.color_seq+(3-blck.epimap)*(~blck.color_seq);
    % run nsim wsls models and take mean performance
    perceived = rand(nsim,blck.ntrl)<perc; % percentage of right outcome categorization
    perceived_outcome = color_outcome.*perceived+(3-color_outcome).*(~perceived);
    wsls_ans = zeros(nsim,blck.ntrl);
    wsls_ans(:,1) = randi(2,nsim,1);
    
    for i = 2:blck.ntrl
        wsls_ans((wsls_ans(:,i-1) == perceived_outcome(:,i-1)),i)=wsls_ans((wsls_ans(:,i-1) == perceived_outcome(:,i-1)),i-1);
        wsls_ans((wsls_ans(:,i-1) ~= perceived_outcome(:,i-1)),i)=3-wsls_ans((wsls_ans(:,i-1) ~= perceived_outcome(:,i-1)),i-1);
    end
end

perf = mean(sum((wsls_ans(:,1:ntrl)==blck.correct(1:ntrl)),2)/ntrl);

end