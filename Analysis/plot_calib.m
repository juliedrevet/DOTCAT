function plot_calib(calibration,hold_on,fit_only)
% plot subject calibration
% calibration must contain following fields
%   * hdr
%   * cfg
%   * stim
%   * rslt


if nargin<2
    hold_on = false;
    fit_only = false;
elseif nargin == 2
    fit_only = false;
end

%% get data from calibration structure
ntrl  = calibration.cfg.ntrl;
nblck = calibration.cfg.nblck;
tmin  = calibration.cfg.rangePerc(1);
tmax  = calibration.cfg.rangePerc(2);
rhat  = calibration.rslt.range;
% in pilot data rtmax was not saved with calibration structure
try
    rtmax = calibration.rslt.rtmax;
catch
    rtmax = 3;
end
color_prop = calibration.stim.color_prop;
rt   = calibration.rslt.rt;
resp = calibration.rslt.resp;
rangeStim = calibration.cfg.rangeStim;
b = calibration.rslt.b;

%% reshape/adapt rslt
rt_fit         = reshape(rt',nblck*ntrl,1);
color_prop_fit = reshape(color_prop',nblck*ntrl,1);
resp_fit       = reshape(resp',nblck*ntrl,1);
idx_keep = (rt_fit < rtmax);
fprintf('%d responses not used for fitting because response time > %3.2f sec\n',length(idx_keep)-sum(idx_keep),rtmax)
color_prop_fit = color_prop_fit(idx_keep);
resp_fit       = resp_fit(idx_keep);

%% plot color proportion distribution
if ~hold_on
    f = figure();
else
    f = gcf;
end
if ~fit_only
    subplot(2,1,1)
    h = histogram(color_prop_fit(:),rangeStim(1):0.025:rangeStim(2));
    xlim([rangeStim(1) rangeStim(2)]);
    
    %xlabel('proportion')
    ylabel('#occurence')
    
    % plot subject response and fit
    subplot(2,1,2)
else
    tempf = figure();
    h = histogram(color_prop_fit(:),rangeStim(1):0.025:rangeStim(2));
    set(tempf,'Visible','Off')
end
figure(f)
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

suptitle(sprintf('Subject %02d',calibration.hdr.subj)); % keep?

if fit_only
    close(tempf);
end
end
