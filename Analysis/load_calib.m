%% load all expe data
function calib = load_calib(clem)

if nargin<1
    clem = false;
end

if clem
    d = dir('../Data/subj_calib_nopat_with_clem');
else
    d = dir('../Data/subj_calib_nopat');
end
if isempty(d)
    error('no data can be imported!');
end

calib = struct();
calib.hdr  = [];
calib.cfg  = [];
calib.stim = [];
calib.rslt = [];

k = 1;
for i = 1:length(d)
    if length(d(i).name)>20
        ex = importdata(fullfile(d(i).folder,d(i).name));
        calib(k) = ex;
        k = k+1;
    end
end

end