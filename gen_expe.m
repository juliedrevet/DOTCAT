function [expe] = gen_expe(subj, taskeq)
%  GEN_EXPE  Generate DOTCAT experiment
%
%  Usage: [expe] = GEN_EXPE(subj,taskeq)
%
%  where subj is the subject number. Experiment parameters are counter-balanced
%  after groups of 8 subjects. The experiment structure expe contains the block
%  substructure blck required to run the experiment and the configuration
%  cfg for this subject
%
%  The optional parameter taskeq controls the degree of equalization between
%  the two tasks (true for full theoretical equalization between the two tasks 
%  (mirror)). The default is true.


addpath ./Toolboxes/Rand
addpath ./Toolboxes/IO
addpath ./Toolboxes/Stimuli/

expe  = struct(); %12 => 4 practice blocks, 8 experimental blocks
cfg   = struct();
nprac = 2;

% check input arguments
if nargin < 2
    % apply full equalization between tasks
    taskeq = true;
elseif nargin < 1
    error('Missing subject number!');
end

%% configuration structure
% cfg.condtn: task condition  => 1:simple Dot or 2: bicolor dot
% cfg.epimap: episode mapping => 1:color1=target|starting or 2:color2=target|starting

% task condition => 1:simple dot or 2:bicolor dot
% 1st bit of subject number
if bitget(subj,1)
    cfg.condtn = [1 2  1 2 1 2  1 2 1 2];
else
    cfg.condtn = [2 1  2 1 2 1  2 1 2 1];
end
% episode mapping => 1:color1=target|starting or 2:color2=target|starting
% 2rd bit of subject number
if bitget(subj,2)
    cfg.epimap = [2 1  2 1 1 2  1 2 2 1];
else
    cfg.epimap = [1 2  1 2 2 1  2 1 1 2];
end

cfg.nprac = nprac;
expe.cfg  = cfg;

nblck = length(cfg.condtn);

%%
b = 1;
fprintf('Generating practice blocks...');
if taskeq
    % generate only 4 and distribute it to equalize between condtn UNI and BICOLOR
    isub = find(cfg.condtn(nprac+1:end) == cfg.condtn(1))+nprac; % block indices of 1st condtn
    jsub = find(cfg.condtn(nprac+1:end) ~= cfg.condtn(1))+nprac; % block indices of 2nd condtn
    
    % idx 1st condt: 1:4
    % idx 2nd condt: mod(((1:4)+1),4)+1 (mirror!)
    
    % generate practice blocks
    for iprac = 1:nprac
        blck = gen_blck(true);
        f = fieldnames(blck);
        blck.condtn = cfg.condtn(iprac);
        blck.epimap = cfg.condtn(iprac);
        blck = orderfields(blck,['condtn'; 'epimap';f]);
        expe.blck(iprac) = blck;
    end
    fprintf('done.\nGenerating testing blocks...\n');
    
    % generate expe blocks (4 needed)
    for iblck = 1:(nblck-nprac)/2
        fprintf('  * generating episodes for block %d/%d...\n',b,(nblck-nprac));
        blck = gen_blck;
        f = fieldnames(blck);
        blck.condtn = cfg.condtn(isub(iblck));
        blck.epimap = cfg.epimap(isub(iblck));
        blck = orderfields(blck,['condtn'; 'epimap';f]);
        expe.blck(isub(iblck)) = blck;
        b = b+1;
        fprintf('  * generating episodes for block %d/%d...\n',b,(nblck-nprac));
        blck.condtn = cfg.condtn(jsub(mod(iblck+1,4)+1));
        blck.epimap = cfg.epimap(jsub(mod(iblck+1,4)+1));
        blck = orderfields(blck,['condtn'; 'epimap';f]);
        expe.blck(jsub(mod(iblck+1,4)+1)) = blck;
        b = b+1;
    end
    fprintf('...done.\n\n');
else
    % no equalization between condtn
    for iblck = 1:nblck
        % generate block
        if iblck <= nprac
            blck = gen_blck(true);
        else
            fprintf('  * generating episodes for block %d/%d...\n',b-nprac,nblck-nprac);
            blck = gen_blck;
        end
        f = fieldnames(blck);
        blck.condtn = cfg.condtn(iblck);
        blck.epimap = cfg.condtn(iblck);
        blck = orderfields(blck,['condtn'; 'epimap';f]);
        expe.blck(iblck) = blck;
        if b == nprac
            fprintf('done.\nGenerating testing blocks...\n');
        elseif b == nblck
            fprintf('...done.\n\n')
        end
        b = b+1;
    end
end

end