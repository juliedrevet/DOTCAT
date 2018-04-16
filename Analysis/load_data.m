%% load all expe data
function expe = load_data(clem)

if nargin<1
    clem = false;
end

if clem
    d = dir('../Data/subj_data_nopat_with_clem');
else
    d = dir('../Data/subj_data_nopat');
end
if isempty(d)
    error('no data can be imported!');
end

expe = struct();
expe.hdr  = [];
expe.cfg  = [];
expe.blck = [];
expe.stim = [];
expe.rslt = [];

k = 1;
for i = 1:length(d)
    if length(d(i).name)>10
        ex = importdata(fullfile(d(i).folder,d(i).name));
        expe(k) = ex;
        k = k+1;
    end
end

end