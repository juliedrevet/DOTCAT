%% plot calib without hist for all subj

calib = load_calib;

for subj = 1:7
    
    figure()
    plot_calib(calib((subj-1)*3+1),true,true)
    plot_calib(calib((subj-1)*3+2),true,true)
    plot_calib(calib((subj-1)*3+3),true,true)
    
    fpath = '../Data/subj_calib_nopat';
    fname = sprintf('calib_%02d',subj);
    savefig(fullfile(fpath,fname));
end