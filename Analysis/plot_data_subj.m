%% plot subj data
function plot_data_subj

%% load subj data
[filename,filepath] = uigetfile('*.mat','Select subject data');
expe = importdata(fullfile(filepath,filename));


%% plot subj resp versus correct
decale = true;

condtntxt = {'unicolor' 'bicolor'};
taskidtxt = {'OBSERVER' 'AGENT'};

figure()
for i = (expe.cfg.nprac+1):length(expe.cfg.taskid)
    subplot(2,4,i-expe.cfg.nprac);
    taskid = expe.cfg.taskid(i);
    condtn = expe.cfg.condtn(i);
    correct = expe.rslt(i).correct;
    if decale
        resp = expe.rslt(i).resp(2:end);
    else
        resp = expe.rslt(i).resp(1:end-1);
    end
    
    plot(correct); 
    hold on; 
    plot(resp,'.'); 
    title(sprintf('Block %d - %s - %s dots',i-expe.cfg.nprac,taskidtxt{taskid},condtntxt{condtn}));
    xlim([1 120]);
    ylim([0.9 2.1]);
    box off
end



end