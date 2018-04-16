function [expe,aborted,errmsg] = run_expe(subj,calibration,pattern_raw,expe_raw,start_blck,end_blck)

addpath ./Toolboxes/Rand
addpath ./Toolboxes/IO
addpath ./Toolboxes/Stimuli/Visual

syncflip = true;
% start_blck include practice blocks!!

% check input arguments
if nargin < 1
    error('Missing subject number!');
elseif nargin == 1
    warning('Will take default calibration!');
    calibration = importdata('./Data/DOTCAT_default_calibration.mat');
    warning('Will take default patterns!');
    pattern_raw = importdata('./Data/DOTCAT_default_pattern_raw.mat');
    warning('Will take default raw experiment!');
    expe_raw    = importdata('./Data/DOTCAT_default_expe_raw.mat');
    start_blck  = 1;
    end_blck    = 10;
elseif nargin == 2
    warning('Will take default patterns!');
    pattern_raw = importdata('./Data/DOTCAT_default_pattern_raw.mat');
    warning('Will take default raw experiment!');
    expe_raw   = importdata('./Data/DOTCAT_default_expe_raw.mat');
    start_blck = 1;
    end_blck   = 10;
elseif nargin == 3
    warning('Will take default raw experiment!');
    expe_raw   = importdata('./Data/DOTCAT_default_expe_raw.mat');
    start_blck = 1;
    end_blck   = 10;
elseif nargin == 4
    start_blck = 1;
    end_blck   = 10;
end

if ~isscalar(subj) || mod(subj,1) ~= 0
    error('Invalid subject number!');
end

%% create data folder for interim subject files
foldname = sprintf('./Data/S%02d',subj);
if ~exist(foldname,'dir')
    mkdir(foldname);
end

hdr = struct();
hdr.subj = subj;
hdr.date = datestr(now,'yyyymmdd-HHMM');

%% options (customizable)
make_scrshots = false; % make screenshots?
if make_scrshots
    warning('Will make screenshots in ./scrshots/S%02d !',subj);
    mkdir(sprintf('./scrshots/Data_S%02d',subj))
end
perf_screen = false;  % show performance screen at the end of the block?
keepshape   = false; % show shapes and choice framed during outcome?
respprob    = true;  % show response probe that disappear just before stim?

%% generate experiment for subject
if isempty(expe_raw)
    expe  = gen_expe(subj);
    shape = shuffle_shape(expe); % create subject shape combinations
    expe_raw       = [];
    expe_raw.expe  = expe;
    expe_raw.shape = shape;
    % save raw expe structure in case of premature termination
    fname = sprintf('DOTCAT_S%02d_%s_raw',subj,datestr(now,'yyyymmdd-HHMM'));
    fname = fullfile(foldname,fname);
    save([fname,'.mat'],'expe_raw');
else % in case of premature termination use raw expe struct already generated
    expe  = expe_raw.expe; % contient deja cfg
    shape = expe_raw.shape;
    nprac = expe.cfg.nprac;
    for b = 1:(start_blck-1) % load already saved blocks in general expe struct
        if b <= nprac
            partname = sprintf('DOTCAT_S%02d_t%02d',subj,b);
        else
            partname = sprintf('DOTCAT_S%02d_b%02d',subj,b-nprac);
        end
        d = dir(sprintf('./Data/S%02d',subj));
        if isempty(d)
            error('no data can be imported!');
        end
        for i = 1:length(d)
            if strncmp(d(i).name,partname,length(partname))
                blck2merge   = importdata(fullfile(d(i).folder,d(i).name));
                expe.blck(b) = blck2merge.blck;
                expe.stim(b) = blck2merge.stim;
                expe.rslt(b) = blck2merge.rslt;
                break;
            end
        end
    end
end

%%
% define output arguments
aborted = false; % aborted prematurely?
errmsg  = []; % error message

% get cfg from calibration structure
cfg = calibration.cfg;

% set screen parameters
iscr = 0;%2  % screen index
res  = []; % screen resolution
fps  = []; % screen refresh rate
ppd  = cfg.ppd; % number of screen pixels per degree of visual angle

% set setup parameters from calibration data
lumibg   = cfg.lumibg;   % background luminance (grey)
fix_siz  = cfg.fix_siz;  % fixation point size
dot_siz  = cfg.dmtr;  % DOT size
probwdth = cfg.probwdth; % response probe width

% set stimulation parameters
shape_siz = 6.0*ppd;        % shape size
shape_off = 7.0*ppd;          % shape offset
shape_add = round(0.1*ppd); % choice rectangle offset

% set instructions parameters
vert_off  = round(1.2*ppd); % vertical offset for instructions screen
info_fac  = 2.5;            % informational text magnification factor
instr_fac = 1.15;           % instructions text magnification factor
instr_siz = 3.5*ppd;        % instructions shapes size
bag_siz = round(6.5*ppd);   % instructions bag size

% set instructions labels
label_esc   = 'appuyez sur [espace] pour continuer';
label_obs   = 'Dans quel sac pioche l''ordinateur?';
label_agent = 'Piochez un maximum de fois dans le sac';
label_wait  = 'veuillez patienter...';

% set list of colors (R/G/B values)
color_frame = [96,96,96]/255;
color_shape = [175,175,175]/255;
colors      = cfg.colors;

% create video structure
video = [];
%%
try
    % hide cursor and stop spilling key presses into MATLAB windows
    HideCursor;
    FlushEvents;
    ListenChar(2);
    
    % check keyboard responsiveness before doing anything
    fprintf('\n');
    fprintf('Press any key to check keyboard responsiveness... ');
    if WaitKeyPress([],30) == 0
        fprintf('\n\n');
        error('No key press detected after 30 seconds.');
    else
        fprintf('Good.\n\n');
    end
    
    % set keys
    KbName('UnifyKeyNames');
    keywait = KbName('space');
    keyquit = KbName('ESCAPE');
    keyresp = KbName({'E','P'});
    
    % open main window
    % set screen resolution and refresh rate
    if ~isempty(res) && ~isempty(fps)
        r = Screen('Resolutions',iscr);
        i = find([r.width] == res(1) & [r.height] == res(2));
        if isempty(i) || ~any([r(i).hz] == fps)
            error('Cannot set screen to %d x %d at %d Hz.',res(1),res(2),fps);
        end
        Screen('Resolution',iscr,res(1),res(2),fps);
    end
    % set screen synchronization properties
    % see 'help SyncTrouble',
    %     'help BeampositionQueries' or
    %     'help ConserveVRAMSettings' for more information
    if syncflip
        if ispc
            % soften synchronization test requirements
            Screen('Preference','SyncTestSettings',[],[],0.2,10);
            % enforce beamposition workaround for missing VBL interval
            Screen('Preference','ConserveVRAM',bitor(4096,Screen('Preference','ConserveVRAM')));
        end
        Screen('Preference','VisualDebuglevel',3);
    else
        % skip synchronization tests altogether
        Screen('Preference','SkipSyncTests',1);
        Screen('Preference','VisualDebuglevel',0);
        Screen('Preference','SuppressAllWarnings',1);
    end
    % set font properties
    if ismac
        txtfnt = 'Arial';
        txtsiz = round(1.0*ppd);
    elseif ispc
        txtfnt = 'Arial'; % closest to Helvetica
        txtsiz = round(2/3*ppd); % text size is ~2/3 smaller in Windows than MacOSX
    end
    Screen('Preference','TextAlphaBlending',1);
    Screen('Preference','DefaultFontName',txtfnt);
    Screen('Preference','DefaultFontSize',txtsiz);
    Screen('Preference','DefaultFontStyle',0);
    % prepare configuration and open main window
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask','General','UseFastOffscreenWindows');
    PsychImaging('AddTask','General','NormalizedHighresColorRange');
    video.i = iscr;
    video.res = Screen('Resolution',video.i);
    video.h = PsychImaging('OpenWindow',video.i,0);
    [video.x,video.y] = Screen('WindowSize',video.h);
    if syncflip
        video.ifi = Screen('GetFlipInterval',video.h,100,50e-6,10);
    else
        video.ifi = 1/60; % assume 60 Hz
    end
    Screen('BlendFunction',video.h,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    Priority(MaxPriority(video.h));
    Screen('ColorRange',video.h,1);
    Screen('FillRect',video.h,lumibg);
    Screen('Flip',video.h);
    
    % open offscreen window
    video.hoff = Screen('OpenOffscreenWindow',video.h);
    %% TEXTURES
    % load shape textures
    shape_tex = zeros(2,8); % shape_tex(1,i)=black (contour) shape_tex(2,i)=greyish (inside)
    instr_tex = zeros(2,8);
    for i = 1:8
        imgc = double(imread(sprintf('./img/shape%dc.png',i)))/255;
        imgc2 = imresize(imgc,instr_siz/size(imgc,1));
        imgc = imresize(imgc,shape_siz/size(imgc,1));
        shape_tex(1,i) = Screen('MakeTexture',video.h,cat(3,ones(size(imgc)),imgc),[],[],2);
        instr_tex(1,i) = Screen('MakeTexture',video.h,cat(3,ones(size(imgc2)),imgc2),[],[],2);
        imgi = double(imread(sprintf('./img/shape%d.png',i)))/255;
        imgi2 = imresize(imgi,instr_siz/size(imgi,1));
        imgi = imresize(imgi,shape_siz/size(imgi,1));
        shape_tex(2,i) = Screen('MakeTexture',video.h,cat(3,ones(size(imgi)),imgi),[],[],2);
        instr_tex(2,i) = Screen('MakeTexture',video.h,cat(3,ones(size(imgi2)),imgi2),[],[],2);
    end
    shape_rec = zeros(4,4); % shape_rec(1,:)=left shape_rec(2,:)=right shape_rec(3,:)=up shape_rec(4,:)=down
    shape_rec(1,:) = CenterRectOnPoint(Screen('Rect',shape_tex(1)),video.x/2-shape_off,video.y/2);
    shape_rec(2,:) = CenterRectOnPoint(Screen('Rect',shape_tex(1)),video.x/2+shape_off,video.y/2);
    shape_rec(3,:) = CenterRectOnPoint(Screen('Rect',shape_tex(1)),video.x/2,video.y/2+vert_off-shape_off/2);
    shape_rec(4,:) = CenterRectOnPoint(Screen('Rect',shape_tex(1)),video.x/2,video.y/2+vert_off+shape_off/2);
    
    % configure the rectangles around the choice
    shape_box = [[ ...
        video.x/2-shape_off-round(RectWidth(shape_rec(1,:)))*0.5-shape_add; ...
        video.y/2-round(RectWidth(shape_rec(1,:)))*0.5-shape_add; ...
        video.x/2-shape_off+round(RectWidth(shape_rec(1,:)))*0.5+shape_add; ...
        video.y/2+round(RectWidth(shape_rec(1,:)))*0.5+shape_add] ...
        [ ...
        video.x/2+shape_off-round(RectWidth(shape_rec(1,:)))*0.5-shape_add; ...
        video.y/2-round(RectWidth(shape_rec(1,:)))*0.5-shape_add; ...
        video.x/2+shape_off+round(RectWidth(shape_rec(1,:)))*0.5+shape_add; ...
        video.y/2+round(RectWidth(shape_rec(1,:)))*0.5+shape_add]];
    
    % create fixation point
    img = CreateCircularAperture(fix_siz);
    fix_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    fix_rec = CenterRectOnPoint(Screen('Rect',fix_tex(1)),video.x/2,video.y/2);
    
    % create response probe contour
    img = CreateCircle(fix_siz+6*probwdth,probwdth);
    resp_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    resp_rec = CenterRectOnPoint(Screen('Rect',resp_tex(1)),video.x/2,video.y/2);
    
    % create single DOT outcome
    img = CreateCircularAperture(dot_siz);
    outc_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    outc_rec = CenterRectOnPoint(Screen('Rect',outc_tex(1)),video.x/2,video.y/2);
    
    % create bag textures 1:Obs | 2:Agent | 3:Agent (in text)
    [~,~,imgbagi] = imread('./img/bag6.png'); % filled bag
    [~,~,imgbagc] = imread('./img/bagc6.png'); % contour bag
    % big bags
    imgbagi = double(imgbagi)/255;
    imgbagi = imresize(imgbagi,bag_siz/size(imgbagi,1));
    imgbagc = double(imgbagc)/255;
    imgbagc = imresize(imgbagc,bag_siz/size(imgbagc,1));
    bag_tex(1,1) = Screen('MakeTexture',video.h,cat(3,ones(size(imgbagc)),imgbagc),[],[],2);
    bag_tex(1,2) = Screen('MakeTexture',video.h,cat(3,ones(size(imgbagi)),imgbagi),[],[],2);
    % small bags
    imgbagi = imresize(imgbagi,.5*bag_siz/size(imgbagi,1));
    imgbagc = imresize(imgbagc,.5*bag_siz/size(imgbagc,1));
    bag_tex(2,1) = Screen('MakeTexture',video.h,cat(3,ones(size(imgbagc)),imgbagc),[],[],2);
    bag_tex(2,2) = Screen('MakeTexture',video.h,cat(3,ones(size(imgbagi)),imgbagi),[],[],2);
    
    % create instructions rects
    % upper text
    Screen('TextSize',video.h,round(txtsiz*instr_fac));
    rec_agent = CenterRectOnPoint(Screen('TextBounds',video.h,label_agent),video.x/2-.5*bag_siz,video.y/7);
    rec_obs   = CenterRectOnPoint(Screen('TextBounds',video.h,label_obs),video.x/2,video.y/7);
    % bag
    bag_rec(1,:) = CenterRectOnPoint(Screen('Rect',bag_tex(1,1)),video.x/2,video.y/2+vert_off-shape_off/2);
    bag_rec(2,:) = CenterRectOnPoint(Screen('Rect',bag_tex(1,1)),video.x/2,video.y/2+vert_off+shape_off/2);
    bag_rec(3,:) = CenterRectOnPoint(Screen('Rect',bag_tex(2,1)),rec_agent(3)+.5*bag_siz,video.y/7);
    % instructions shapes inside of bag
    instr_rec(1,:) = CenterRectOnPoint(Screen('Rect',instr_tex(1)),video.x/2,bag_rec(1,2)+.65*bag_siz);
    instr_rec(2,:) = CenterRectOnPoint(Screen('Rect',instr_tex(1)),video.x/2,bag_rec(2,2)+.65*bag_siz);
    
    %% WELCOME SCREEN
    % first flip
    t = Screen('Flip',video.h);
    
    Screen('TextStyle',video.h,0);
    Screen('TextSize',video.h,round(txtsiz));
    
    labeltxt = sprintf('appuyez sur [espace] pour démarrer');
    labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,round(1.2*ppd));
    Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
    
    Screen('TextSize',video.h,round(txtsiz*info_fac));
    labeltxt = sprintf('Bienvenue!');
    labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y/2);
    Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
    Screen('DrawingFinished',video.h);
    Screen('Flip',video.h,t+roundfp(1.500,0.500));
    WaitKeyPress(keywait);
    
    % draw fixation point
    Screen('DrawTexture',video.h,fix_tex,[],fix_rec,[],[],[],0);
    Screen('DrawingFinished',video.h);
    t = Screen('Flip',video.h);
    
    nblck      = length(expe.cfg.condtn); % total number of blocks
    nblck_prac = expe.cfg.nprac; % number of practice blocks
    pseudopos  = 3;
    color_prop = flip(calibration.rslt.range);
    
    %% loop on blocks
    for iblck = start_blck:end_blck
        % blck.epimap:    starting color    => 1:color1 or 2:color2
        % blck.color_seq: color profile of outcome => 1:epimap color or 0:other color
        blck = expe.blck(iblck);
        ntrl = blck.ntrl;
        clear pattern;
        pattern = Pattern();
        
        % create stim substructure
        stim = [];
        % shapes for this block shape(1): color1; shape(2):color2
        stim.shape = shape(iblck,:);
        
        %% show instructions screen
        % end of training TODO: à mettre à la fin pour plus de clareté
        if (iblck == (nblck_prac+1)) && (nblck_prac~=0)
            Screen('TextSize',video.h,round(txtsiz*info_fac));
            labeltxt = 'fin de l''entraînement!';
            labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y/2);
            Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
            Screen('TextSize',video.h,round(txtsiz*.7));
            labeltxt = label_esc;
            labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y-round(1.2*ppd));
            Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h,t+roundfp(1));
            WaitKeyPress(keywait);
        end
        
        
        % position of shape following color_seq(_ffb)
        pos = zeros(1,ntrl+1);
        while abs(diff(hist(pos,1:2)))>1
            pos(1:pseudopos) = randi(2,1,pseudopos);
            for i = (pseudopos+1):(ntrl+1)
                ipos = ceil(2*rand(1));
                if sum(pos((i-pseudopos):(i-1))==ipos)>=pseudopos
                    pos(i) = 3-ipos;
                else
                    pos(i) = ipos;
                end
            end
        end
        stim.pos = pos;
        
        % theoretical color outcome => 1: color1 or 2: color2
        if blck.condtn == 1 % simple dot outcome: include false feedbacks
            color_outcome = blck.epimap*blck.colorffb_seq+(3-blck.epimap)*(~blck.colorffb_seq);
        elseif blck.condtn == 2 % bicolor dot outcome: no false feedbacks
            color_outcome = blck.epimap*blck.color_seq+(3-blck.epimap)*(~blck.color_seq);
        end
        
        stim.outcome = color_outcome;
        
        
        % correct: without ffb
        stim.correct = blck.epimap*(blck.reward_seq==blck.p_reward)+(3-blck.epimap)*(~(blck.reward_seq==blck.p_reward));
        
        ntrl = 5;
        
        % regular instructions
        iu = randi(2); % randomize position of upper shape/color
        id = 3-iu;
        draw_instr(iu,id,'wait');
        Screen('DrawingFinished',video.h);
        t = Screen('Flip',video.h,t+roundfp(2));
        
        % bicolor dot pattern
        if blck.condtn == 2
            pattern(1,ntrl) = Pattern();
            tic
            for k = 1:ntrl
                pattern(k) = color_polygon(pattern_raw(iblck,k),color_prop(color_outcome(k)),1e4,1e-4,10);
            end
            toc
            stim.pattern = pattern;
        else
            stim.pattern = [];
        end
        
        draw_instr(iu,id,'esc');
        Screen('DrawingFinished',video.h);
        t = Screen('Flip',video.h,t+roundfp(ntrl*.5));
        if make_scrshots
            imgscreen = Screen('GetImage',video.h);
            if (iblck<=nblck_prac)
                imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/instructions_training%d_%s.png',...
                    subj,iblck,datestr(now,'yymmdd-HHMMSS')));
            else
                imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/instructions_block%d_%s.png',...
                    subj,iblck-nblck_prac,datestr(now,'yymmdd-HHMMSS')));
            end
        end
        
        WaitKeyPress(keywait);
        
        % create results substructure
        rslt        = [];
        rslt.resp   = zeros(1,ntrl+1); % response as stimulus index
        rslt.respkb = zeros(1,ntrl+1); % response as keyboard index
        rslt.shape  = zeros(1,ntrl+1); % response as shape index
        rslt.rt     = zeros(1,ntrl+1); % response time
        
        %% loop on trials
        for itrl = 1:(ntrl+1)
            
            % draw fixation point
            Screen('DrawTexture',video.h,fix_tex,[],fix_rec,[],[],[],0);
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h);
            
            %check if abort key is pressed
            if CheckKeyPress(keyquit)
                aborted = true;
                break
            end
            
            % left/right shape index
            if stim.pos(itrl) == 1
                il = blck.epimap;
                ir = 3-il;
            else
                il = 3-blck.epimap;
                ir = blck.epimap;
            end
            
            Screen('DrawTexture',video.h,fix_tex,[],fix_rec,[],[],[],0);
            draw_stim(il,ir);
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h,t+roundfp(.6,0.200));
            
            if make_scrshots % screenshot for stimulus shown waiting for response
                imgscreen = Screen('GetImage',video.h);
                if (iblck<=nblck_prac)
                    imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/training%d_%s_stim%d.png',...
                        subj,iblck,datestr(now,'yymmdd-HHMMSS'),itrl));
                else
                    imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/block%d_%s_stim%d.png',...
                        subj,iblck-nblck_prac,datestr(now,'yymmdd-HHMMSS'),itrl));
                end
            end
            
            % keyboard input
            [response,tkey] = WaitKeyPress(keyresp,[],false);
            rslt.rt(itrl)   = tkey-t;
            
            if response == 1 % left shape chosen
                rslt.resp(itrl)   = il;
                rslt.respkb(itrl) = 1;
                rslt.shape(itrl)  = stim.shape(il);
            else             % right shape chosen
                rslt.resp(itrl)   = ir;
                rslt.respkb(itrl) = 2;
                rslt.shape(itrl)  = stim.shape(ir);
            end
            Screen('DrawTexture',video.h,fix_tex,[],fix_rec,[],[],[],0);
            draw_stim(il,ir);
            Screen('FrameRect',video.h,color_frame,shape_box(:,response),8); % frame the shape chosen by participant
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            
            if make_scrshots % screenshot for response done
                imgscreen = Screen('GetImage',video.h);
                if (iblck<=nblck_prac)
                    imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/training%d_%s_resp%d.png',...
                        subj,iblck,datestr(now,'yymmdd-HHMMSS'),itrl));
                else
                    imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/block%d_%s_resp%d.png',...
                        subj,iblck-nblck_prac,datestr(now,'yymmdd-HHMMSS'),itrl));
                end
                
            end
            
            WaitSecs(.4);
            
            % time lapse before getting ready for outcome
            if keepshape
                draw_stim(il,ir);
                Screen('FrameRect',video.h,color_frame,shape_box(:,response),8);
            end
            Screen('DrawTexture',video.h,fix_tex,[],fix_rec,[],[],[],0);
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            WaitSecs(.15);
            
            if itrl ~=(ntrl+1)
                % show response prob to focus on outcome
                if keepshape
                    draw_stim(il,ir);
                    Screen('FrameRect',video.h,color_frame,shape_box(:,response),8);
                end
                Screen('DrawTexture',video.h,fix_tex,[],fix_rec,[],[],[],0);
                if respprob == true
                    Screen('DrawTexture',video.h,resp_tex,[],resp_rec,[],[],[],0);
                end
                Screen('DrawingFinished',video.h);
                Screen('Flip',video.h);
                WaitSecs(.35);
                
                % time lapse without response prob before outcome appears
                if keepshape
                    draw_stim(il,ir);
                    Screen('FrameRect',video.h,color_frame,shape_box(:,response),8);
                end
                Screen('DrawTexture',video.h,fix_tex,[],fix_rec,[],[],[],0);
                Screen('DrawingFinished',video.h);
                Screen('Flip',video.h);
                WaitSecs(.25);
                
                %% draw outcome
                if keepshape
                    draw_stim(il,ir);
                    Screen('FrameRect',video.h,color_frame,shape_box(:,response),8);
                end
                
                if blck.condtn == 1
                    Screen('DrawTexture',video.h,outc_tex,[],outc_rec,[],[],[],colors(color_outcome(itrl),:));
                else
                    draw_pat(pattern(itrl));
                end
                
                Screen('DrawingFinished',video.h);
                Screen('Flip',video.h);
                
                if make_scrshots %&&(iblck<=nblck_prac) % screenshot for outcome
                    imgscreen = Screen('GetImage',video.h);
                    if (iblck<=nblck_prac)
                        imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/training%d_%s_outcome%d.png',...
                            subj,iblck,datestr(now,'yymmdd-HHMMSS'),itrl));
                    else
                        imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/block%d_%s_outcome%d.png',...
                            subj,iblck-nblck_prac,datestr(now,'yymmdd-HHMMSS'),itrl));
                    end
                end
                
                WaitSecs(5*video.ifi);
                
                % hide outcome
                Screen('DrawTexture',video.h,fix_tex,[],fix_rec,[],[],[],0);
                Screen('DrawingFinished',video.h);
                Screen('Flip',video.h);
                WaitSecs(.2);
                
                Screen('DrawTexture',video.h,fix_tex,[],fix_rec,[],[],[],0);
                Screen('DrawingFinished',video.h);
                t = Screen('Flip',video.h);
                
            end
            
        end % end of trial loop
        
        % store results and block substructures
        rslt.correct   = stim.correct(1:ntrl);
        rslt.perf      = sum(rslt.resp(1:ntrl)==stim.correct(1:ntrl))/ntrl;
        rslt.perfafter = sum(rslt.resp(2:ntrl+1)==stim.correct(1:ntrl))/ntrl;
        rslt.wsls_perf = run_model(blck,stim, 1000, ntrl, cfg.rangePerc(2));
        
        if  perf_screen == true
            WaitSecs(.9);
            Screen('TextSize',video.h,round(txtsiz));
            labeltxt = sprintf('%d%% de bonnes réponses pour ce bloc, continuez!',round(rslt.perf*100));
            labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y/2);
            Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
            %             if rslt.perf>=rslt.wsls_perf
            %                 labeltxt = 'en bonne voie pour le bonus';
            %                 labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,3*video.y/4);
            %                 Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
            %             end
            Screen('TextSize',video.h,round(txtsiz*.7));
            labeltxt = label_esc;
            labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y-round(1.2*ppd));
            Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h);
            WaitKeyPress(keywait);
        elseif (perf_screen == false) && (iblck>nblck_prac)
            WaitSecs(.9);
            Screen('TextSize',video.h,round(txtsiz));
            labeltxt = sprintf('fin du bloc %d, bravo!',(iblck-nblck_prac));
            labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y/2);
            Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
            Screen('TextSize',video.h,round(txtsiz*.7));
            labeltxt = label_esc;
            labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y-round(1.2*ppd));
            Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h);
            WaitKeyPress(keywait);
        end
        
        expe.stim(iblck) = stim;
        expe.rslt(iblck) = rslt;
        
        % create subblock structure for saving purposes
        expe_blck = struct();
        expe_blck.cfg = cfg;
        expe_blck.stim = stim;
        expe_blck.blck = blck;
        
        expe_blck.rslt = rslt;
        
        % save temporary file (block per block)
        fpath = foldname;
        if (iblck<=nblck_prac)
            fname = sprintf('DOTCAT_S%02d_t%02d_%s',subj,iblck,datestr(now,'yyyymmdd-HHMM'));
        else
            fname = sprintf('DOTCAT_S%02d_b%02d_%s',subj,(iblck-nblck_prac),datestr(now,'yyyymmdd-HHMM'));
        end
        if aborted
            if ~exist([fpath,'/aborted'],'dir')
                mkdir([fpath,'/aborted']);
            end
            fpath = [foldname,'/aborted/'];
            fname = [fname,'_aborted'];
        end
        fname = fullfile(fpath,fname);
        save([fname,'.mat'],'expe_blck');
        
        if aborted
            break
        end
        
    end % end of block loop
    
    if aborted
        Priority(0);
        Screen('CloseAll');
        FlushEvents;
        ListenChar(0);
        ShowCursor;
        return
    end
    
    if end_blck == nblck
        perf_tot = 0;
        perf_wsls_tot = 0;
        for i = (nblck_prac+1):nblck
            perf_tot = perf_tot+expe.rslt(i).perf/(nblck-nblck_prac);
            perf_wsls_tot = perf_wsls_tot+expe.rslt(i).wsls_perf/(nblck-nblck_prac);
        end
        
        fprintf('Win-stay-Loose-switch performance: %d%%\n', round(perf_wsls_tot*100))
        fprintf('Your performance: %d%%', round(perf_tot*100))
        if (perf_tot>=perf_wsls_tot)
            fprintf(' ... 5 euros de bonus!\n\n')
        end
    end
    
    %% draw end screen
    if end_blck == nblck
        Screen('TextStyle',video.h,0);
        Screen('TextSize',video.h,round(txtsiz));
        labeltxt = sprintf('L''expérience est terminée');
        labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,round(1.2*ppd));
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        
        Screen('TextSize',video.h,round(txtsiz*info_fac));
        labeltxt = 'Merci pour votre participation!';
        labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y/2);
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        
        Screen('TextSize',video.h,round(txtsiz*1.5));
    else
        Screen('TextSize',video.h,round(txtsiz*info_fac));
        labeltxt = 'fin de la première partie';
        labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y/2);
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
    end
    Screen('DrawingFinished',video.h);
    Screen('Flip',video.h,t+roundfp(1.500,0.500));
    WaitKeyPress(keywait);
    
    % save temporary file (whole expe in one .mat file)
    expe.hdr = hdr;
    expe = orderfields(expe,{'hdr','cfg','blck','stim','rslt'});
    fpath = foldname;
    fname = sprintf('DOTCAT_S%02d_%s',subj,datestr(now,'yyyymmdd-HHMM'));
    fname = fullfile(fpath,fname);
    save([fname,'.mat'],'expe');
    
    
    % close Psychtoolbox
    Priority(0);
    Screen('CloseAll');
    FlushEvents;
    ListenChar(0);
    ShowCursor;
catch ME
    
    % close Psychtoolbox
    Priority(0);
    Screen('CloseAll');
    FlushEvents;
    ListenChar(0);
    ShowCursor;
    
    % handle error
    if nargout > 2
        errmsg = ME;
        %errmsg = rmfield(errmsg,'stack');
    else
        rethrow(ME);
    end
    
end

%%
    function [t] = roundfp(t,dt)
        % apply duration rounding policy for video flips
        % where t  - desired (input)/rounded (output) duration
        %       dt - desired uniform jitter on duration (default: none)
        n = round(t/video.ifi);
        % apply uniform jitter
        if nargin > 1 && dt > 0
            m = round(dt/video.ifi);
            n = n+ceil((m*2+1)*rand)-(m+1);
        end
        % convert frames to duration
        t = (n-0.5)*video.ifi;
    end

%%
    function draw_stim(il,ir)
        % draw left stimulus
        is = stim.shape(il);
        Screen('DrawTexture',video.h,shape_tex(1,is),[],shape_rec(1,:),[],[],[],0);
        Screen('DrawTexture',video.h,shape_tex(2,is),[],shape_rec(1,:),[],[],[],color_shape);
        % draw right stimulus
        is = stim.shape(ir);
        Screen('DrawTexture',video.h,shape_tex(1,is),[],shape_rec(2,:),[],[],[],0);
        Screen('DrawTexture',video.h,shape_tex(2,is),[],shape_rec(2,:),[],[],[],color_shape);
    end
%%
    function draw_instr(iu,id,txt)
        
        Screen('TextSize',video.h,round(txtsiz*.7));
        if strcmp(txt,'esc')
            rec_esc = CenterRectOnPoint(Screen('TextBounds',video.h,label_esc),video.x/2,video.y-round(1.2*ppd));
            Screen('DrawText',video.h,label_esc,rec_esc(1),rec_esc(2),0);
        elseif strcmp(txt,'wait')
            rec_wait = CenterRectOnPoint(Screen('TextBounds',video.h,label_wait),video.x/2,video.y-round(1.2*ppd));
            Screen('DrawText',video.h,label_wait,rec_wait(1),rec_wait(2),0);
        end
        if iblck>nblck_prac
            Screen('TextSize',video.h,round(txtsiz));
            label_blck = sprintf('Bloc %d/8',iblck-nblck_prac);
            rec_blck = CenterRectOnPoint(Screen('TextBounds',video.h,label_blck),5*video.x/6,video.y-round(1.2*ppd));
            Screen('DrawText',video.h,label_blck,rec_blck(1),rec_blck(2),0);
        end
        
        Screen('TextSize',video.h,round(txtsiz*instr_fac));
        Screen('DrawText',video.h,label_obs,rec_obs(1),rec_obs(2),0);
        Screen('DrawTexture',video.h,bag_tex(1,2),[],bag_rec(1,:),[],[],[],colors(iu,:));
        Screen('DrawTexture',video.h,bag_tex(1,1),[],bag_rec(1,:),[],[],[],0);
        Screen('DrawTexture',video.h,bag_tex(1,2),[],bag_rec(2,:),[],[],[],colors(id,:));
        Screen('DrawTexture',video.h,bag_tex(1,1),[],bag_rec(2,:),[],[],[],0);
        
        % draw up stimulus
        is = stim.shape(iu);
        Screen('DrawTexture',video.h,instr_tex(1,is),[],instr_rec(1,:),[],[],[],0);
        Screen('DrawTexture',video.h,instr_tex(2,is),[],instr_rec(1,:),[],[],[],color_shape);
        
        % draw down stimulus
        is = stim.shape(id);
        Screen('DrawTexture',video.h,shape_tex(1,is),[],instr_rec(2,:),[],[],[],0);
        Screen('DrawTexture',video.h,shape_tex(2,is),[],instr_rec(2,:),[],[],[],color_shape);
    end

    function draw_pat(pat)
        for ipoly = 1: numel(pat.poly_list)
            Screen('FillPoly', video.h,colors(pat.poly_list(ipoly).color_index,:),...
                [pat.poly_list(ipoly).vertex_list(:,1) pat.poly_list(ipoly).vertex_list(:,2)]+[video.x/2 video.y/2],1);
        end
    end
end