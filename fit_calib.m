function calibration = fit_calib(calibration,rtmax,plotornot,hold_on)
% fit calibration data to have psychometric curve
% use rtmax as max allowed response time, rt > rtmax ignored for fitting

if nargin<2
    rtmax = 4;%seconds
    plotornot = false;
    hold_on   = false;
elseif nargin == 2
    plotornot = false;
    hold_on   = false;
elseif nargin == 3
    hold_on   = false;
end

nblck = calibration.cfg.nblck;
ntrl  = calibration.cfg.ntrl;
color_prop = calibration.stim.color_prop;
rt = calibration.rslt.rt;
resp = calibration.rslt.resp;
rangeStim = calibration.cfg.rangeStim;

% remove data if response time too long
calibration.rslt.rtmax = rtmax;

try
    % reshape rslt for fit purposes
    rt_fit         = reshape(rt',nblck*ntrl,1);
    color_prop_fit = reshape(color_prop',nblck*ntrl,1);
    resp_fit       = reshape(resp',nblck*ntrl,1);
    idx_keep = (rt_fit < rtmax);
    fprintf('%d responses not used for fitting because response time > %3.2f sec\n',length(idx_keep)-sum(idx_keep),rtmax)
    color_prop_fit = color_prop_fit(idx_keep);
    resp_fit = resp_fit(idx_keep);
    [b,~,stat]     = glmfit(color_prop_fit,resp_fit == 1,'binomial','link','probit');
    calibration.rslt.b    = b;
    calibration.rslt.stat = stat;
    
    tmin = calibration.cfg.rangePerc(1);
    tmax = calibration.cfg.rangePerc(2);
    try
        rhat = [fzero(@(xx)normcdf(b(1)+xx*b(2))-tmin,[0,1]) ...
            fzero(@(xx)normcdf(b(1)+xx*b(2))-tmax,[0,1])];
        calibration.rslt.range = rhat;
    catch
        fprintf('No values found for %02d%%|%02d%% correct perception!\n',round(tmin*100),round(tmax*100))
        calibration.rslt.range = nan(1,2);
    end
    
    % plot calibration results
    if plotornot
        if ~hold_on
            figure()
        end
        subplot(2,1,1)
        h = histogram(color_prop_fit(:),rangeStim(1):0.025:rangeStim(2));
        xlim([rangeStim(1) rangeStim(2)]);        
        ylabel('#occurence')
        
        % plot subject response and fit
        subplot(2,1,2)
        hold on
        xlim([rangeStim(1) rangeStim(2)])
        plot(0:0.01:1,normcdf(b(1)+(0:0.01:1)*b(2)));
        bins = discretize(color_prop_fit,h.BinEdges);
        categ = zeros(1,h.NumBins);
        for ibin = 1:h.NumBins
            idx = find(bins == ibin);
            categ(ibin) = sum(resp_fit(idx) == 1)/length(idx);
        end
        xabs = h.BinEdges+h.BinWidth/2;
        plot(xabs(1:end-1),categ,'o');
        
        line([0 rhat(1)],[tmin tmin],'LineStyle',':')
        line([0 rhat(2)],[tmax tmax],'LineStyle',':')
        line([rhat(1) rhat(1)],[0 tmin],'LineStyle',':')
        line([rhat(2) rhat(2)],[0 tmax],'LineStyle',':')
        plot(rhat(1),tmin,'b+')
        plot(rhat(2),tmax,'b+')
        ylabel('p(correct)')
        
        xlabel('Proportion of color 1 in DOT pattern')
    end
    
catch
    fprintf('No fit for psychometric curve found!\n')
end

end