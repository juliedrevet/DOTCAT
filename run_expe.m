function [expe,aborted,errmsg] = run_expe(subj,calibration,syncflip,expe_raw,start_blck)

addpath ./Toolboxes/Rand
addpath ./Toolboxes/IO
addpath ./Toolboxes/Stimuli/Visual

% check input arguments
if nargin < 1
    error('Missing subject number!');
elseif nargin == 1
    warning('Will take default calibration!');
    calibration = importdata('./Data/default_calibration.mat');
    syncflip    = true;
    expe_raw    = [];
    start_blck = 1;
elseif nargin == 2
    syncflip = true;
    expe_raw = [];
    start_blck = 1;
elseif nargin == 3
    expe_raw   = [];
    start_blck = 1;
end

if ~isscalar(subj) || mod(subj,1) ~= 0
    error('Invalid subject number!');
end

%% create data folder for interim subject files
foldname = sprintf('./Data/S%02d',subj);
if ~exist(foldname,'dir')
    mkdir(foldname);
end

%% options (customizable)
make_scrshots = false; % make screenshots?
if make_scrshots
    warning('Will make screenshots in ./scrshots/S%02d !',subj);
    mkdir(sprintf('./scrshots/Data_S%02d',subj))
end
perf_screen = true;  % show performance screen at the end of the block?
prob        = false;  % show handful probe contour during handful outcome?
keepshape   = false; % show shapes and choice framed during outcome?
respprob    = true;  % show response probe that disappear just before stim?
pedestal    = true;  % show handful pedestal?

%% generate experiment for subject
if isempty(expe_raw)
    expe  = gen_expe(subj);
    shape = shuffle_shape(expe); % create subjet shape combinations
    expe_raw       = [];
    expe_raw.expe  = expe;
    expe_raw.shape = shape;
    % save raw expe structure in case of premature termination
    fname = sprintf('DOTCAT_S%02d_%s_raw',subj,datestr(now,'yyyymmdd-HHMM'));
    fname = fullfile(foldname,fname);
    save([fname,'.mat'],'expe_raw');
else % in case of premature termination use raw expe struct already generated
    expe  = expe_raw.expe;
    shape = expe_raw.shape;
    nprac = expe(1).cfg.nprac;
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
                expe(b).cfg  = blck2merge.cfg;
                expe(b).blck = blck2merge.blck;
                expe(b).rslt = blck2merge.rslt;
                break;
            end
        end
    end
end

%%
% define output arguments
aborted = false; % aborted prematurely?
errmsg  = []; % error message

% set screen parameters
iscr = 0;  % screen index
res  = []; % screen resolution
fps  = []; % screen refresh rate
ppd  = calibration.setup.ppd; % number of screen pixels per degree of visual angle

% set setup parameters from calibration data
lumibg    = calibration.setup.lumibg;    % background luminance (grey)
fixtn_siz = calibration.setup.fixtn_siz; % fixation point size
dot_siz   = calibration.setup.dot_siz;   % dot size
hand_siz  = calibration.setup.hand_siz;  % handful size (diameter)
probwdth  = calibration.setup.probwdth;  % handful probe contour width

% set stimulation parameters
shape_siz = 6.0*ppd;        % shape size
shape_off = 7*ppd;          % shape offset
shape_add = round(0.1*ppd); % choice rectangle offset

% set instructions parameters
vert_off  = round(1.2*ppd); % vertical offset for instructions screen
info_fac  = 2.5;            % informational text magnification factor
instr_fac = 1.15;           % instructions text magnification factor
txt_off   = 7*ppd;          % instructions text offset left/right from shape (Observer)
b_instr = 2.8*ppd;          % instructions bag size (diameter)
d_instr = .3*ppd;           % instructions dots in bag size
k_instr = 13;               % number of epimap color dots in the instructions bag
n_instr = 19;               % total number of Dots in the instructions bag

% set instructions labels
label_esc    = 'appuyez sur [espace] pour continuer';
label_choice = {'étiquette' 'pour le sac' 'étiquette pour le sac'};

label_obs    = 'Dans quel sac va piocher l''ordinateur?'; % add 'Observer' ?
label_agent  = 'Piochez un maximum de fois dans le sac'; % add 'Agent' ?

% set list of colors (R/G/B values)
color_frame = [96,96,96]/255;
color_shape = [175,175,175]/255;
color_pedestal = calibration.setup.color_pedestal;
colors      = calibration.setup.colors;

% set handful of dots parameters
nDots  = calibration.setup.nDots;
rangeDots = [calibration.setup.maxDots calibration.setup.minDots];

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
    %%
    % load shape textures
    shape_tex = zeros(2,8); % shape_tex(1,i)=black (contour) shape_tex(2,i)=greyish (inside)
    for i = 1:8
        imgc = double(imread(sprintf('./img/shape%dc.png',i)))/255;
        imgc = imresize(imgc,shape_siz/size(imgc,1));
        shape_tex(1,i) = Screen('MakeTexture',video.h,cat(3,ones(size(imgc)),imgc),[],[],2);
        imgi = double(imread(sprintf('./img/shape%d.png',i)))/255;
        imgi = imresize(imgi,shape_siz/size(imgi,1));
        shape_tex(2,i) = Screen('MakeTexture',video.h,cat(3,ones(size(imgi)),imgi),[],[],2);
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
    img = CreateCircularAperture(fixtn_siz);
    fixtn_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    fixtn_rec = CenterRectOnPoint(Screen('Rect',fixtn_tex(1)),video.x/2,video.y/2);
    
    % create pedestal
    img = CreateCircularAperture(hand_siz+1.5*dot_siz);
    ped_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    ped_rec = CenterRectOnPoint(Screen('Rect',ped_tex(1)),video.x/2,video.y/2);
    
    % create response probe contour
    img = CreateCircle(fixtn_siz+6*probwdth,probwdth);
    resp_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    resp_rec = CenterRectOnPoint(Screen('Rect',resp_tex(1)),video.x/2,video.y/2);
    
    % create feedback
    img = CreateCircularAperture(dot_siz);
    fb_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    fb_rec = CenterRectOnPoint(Screen('Rect',fb_tex(1)),video.x/2,video.y/2);
    
    % create handful probe contour
    img = CreateCircle(hand_siz+1.5*dot_siz,probwdth);
    hand_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    hand_rec = CenterRectOnPoint(Screen('Rect',hand_tex(1)),video.x/2,video.y/2);
    
    % create dots bag for instructions
    imbag = CreateCircle(b_instr+2*d_instr+4*probwdth,4*probwdth);
    bag_tex = Screen('MakeTexture',video.h,cat(3,ones(size(imbag)),imbag),[],[],2);
    
    % create instructions rects
    % upper text
    Screen('TextSize',video.h,round(txtsiz*instr_fac));
    rec_agent = CenterRectOnPoint(Screen('TextBounds',video.h,label_agent),video.x/2-b_instr,video.y/7);
    rec_obs   = CenterRectOnPoint(Screen('TextBounds',video.h,label_obs),video.x/2,video.y/7);
    % 'étiquette ... pour le sac' text (Observer only)
    Screen('TextSize',video.h,txtsiz);
    % define text position (x,y) (DrawText options "x" "y" define the text pen start location)
    % workaround because text left and right are not aligned otherwise
    choice_bounds = zeros(3,4); %label_choice = {'étiquette' 'pour le sac' 'étiquette pour le sac'};
    choice_bounds(1,:) = Screen('TextBounds',video.h,label_choice{1});
    choice_bounds(2,:) = Screen('TextBounds',video.h,label_choice{2});
    choice_bounds(3,:) = Screen('TextBounds',video.h,label_choice{3});
    
    choice_xy = zeros(4,2);
    choice_xy(1,1) = video.x/2-txt_off;
    choice_xy(1,2) = video.y/2+vert_off-shape_off/2-round((choice_bounds(2,4))/2);
    choice_xy(2,1) = video.x/2+txt_off-choice_bounds(1,3);
    choice_xy(2,2) = video.y/2+vert_off-shape_off/2-round((choice_bounds(2,4))/2)+choice_bounds(3,4)-choice_bounds(2,4);
    choice_xy(3,1) = video.x/2-txt_off;
    choice_xy(3,2) = video.y/2+vert_off+shape_off/2-round((choice_bounds(2,4))/2);
    choice_xy(4,1) = video.x/2+txt_off-choice_bounds(1,3);
    choice_xy(4,2) = video.y/2+vert_off+shape_off/2-round((choice_bounds(2,4))/2)+choice_bounds(3,4)-choice_bounds(2,4);
    % instruction bag x,y center coordinates
    bag_xy = zeros(3,2);
    bag_xy(1,1) = rec_agent(3)+b_instr;
    bag_xy(1,2) = video.y/7;
    bag_xy(2,1) = choice_xy(2,1)+choice_bounds(2,3)+b_instr;
    bag_xy(2,2) = video.y/2+vert_off-shape_off/2;
    bag_xy(3,1) =choice_xy(2,1)+choice_bounds(2,3)+b_instr;
    bag_xy(3,2) = video.y/2+vert_off+shape_off/2;
    
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
    Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
    Screen('DrawingFinished',video.h);
    t = Screen('Flip',video.h);
    
    nblck      = length(expe); % total number of blocks
    nblck_prac = expe(1).cfg.nprac; % number of practice blocks
    
    %% loop on blocks
    for iblck = start_blck:nblck
        % blck.epimap:    starting|target (O|A) color    => 1:green or 2:blue
        % blck.color_seq: color profile of outcome|shape => 1:epimap color or 0:other color
        
        blck = expe(iblck).blck;
        % shapes for this block shape(1): green; shape(2):blue
        blck.shape = shape(iblck,:);
        % position of shape following color_seq(_ffb)
        blck.pos = randi(2,1,blck.ntrl+1); %???
        
        % theoretical color outcome => 1: green or 2: blue
        if blck.condtn == 1 % single dot outcome: include false feedbacks
            color_outcome = blck.epimap*blck.colorffb_seq+(3-blck.epimap)*(~blck.colorffb_seq);
        elseif blck.condtn == 2 % handful dots outcome: no false feedbacks
            color_outcome = blck.epimap*blck.color_seq+(3-blck.epimap)*(~blck.color_seq);
        end
        
        % practical color outcome => 1: green or 2: blue
        % line 1: left shape ; line 2: right shape chosen (necessary for Agent)
        outcome = zeros(2,blck.ntrl);
        if blck.taskid == 1 % observer: outcome INDEPENDENT of choice
            outcome = repmat(color_outcome,2,1); % both lines same outcome
        else % agent: outcome DEPENDENT of chosen shape position
            for ipos = 1:(length(blck.pos)-1) %???
                outcome(blck.pos(ipos),ipos)   = color_outcome(ipos);
                outcome(3-blck.pos(ipos),ipos) = 3-color_outcome(ipos);
            end
        end
        
        blck.outcome = outcome;
        % correct: without ffb
        blck.correct = blck.epimap*(blck.reward_seq==blck.p_reward)+(3-blck.epimap)*(~(blck.reward_seq==blck.p_reward));
        
        % end of training
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
            Screen('Flip',video.h,t+roundfp(1.500,0.500));
            WaitKeyPress(keywait);
        end
        %% show instructions screen
        iu = randi(2); % randomize position of upper shape/color
        id = 3-iu;
        
        draw_instr(iu,id);
        Screen('DrawingFinished',video.h);
        Screen('Flip',video.h,t+roundfp(1.500,0.500));
        if make_scrshots
            imgscreen = Screen('GetImage',video.h);
            if (iblck<=nblck_prac)
                imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/instructions_training%d_taskid%d_%s.png',...
                    subj,iblck,blck.taskid,datestr(now,'yymmdd-HHMMSS')));
            else
                imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/instructions_block%d_taskid%d_%s.png',...
                    subj,iblck-nblck_prac,blck.taskid,datestr(now,'yymmdd-HHMMSS')));
            end
        end
        WaitKeyPress(keywait);
        
        %% loop on trials
        %ntrl = 3;
        ntrl = blck.ntrl;
        
        % create results substructure
        rslt         = [];
        rslt.resp    = zeros(1,ntrl+1); % response as stimulus index
        rslt.respkb  = zeros(1,ntrl+1); % response as keyboard index
        rslt.shape   = zeros(1,ntrl+1); % response as shape index
        rslt.rt      = zeros(1,ntrl+1); % response time
        
        for itrl = 1:(ntrl+1) %???
            
            % draw fixation point
            if pedestal && (blck.condtn == 2)
                Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
            end
            Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
            if prob && (blck.condtn == 2)
                Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
            end
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h);
            
            %check if abort key is pressed
            if CheckKeyPress(keyquit)
                Screen('CloseAll')
                aborted = true;
                ShowCursor;
                sca;
                break
            end
            
            % left/right shape index
            if blck.pos(itrl) == 1
                il = blck.epimap;
                ir = 3-il;
            else
                il = 3-blck.epimap;
                ir = blck.epimap;
            end
            
            if pedestal && (blck.condtn == 2)
                Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
            end
            Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
            draw_stim(il,ir);
            if prob && (blck.condtn == 2)
                Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
            end
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h,t+roundfp(.6,0.200));
            
            if make_scrshots % screenshot for stimulus shown waiting for response
                imgscreen = Screen('GetImage',video.h);
                if (iblck<=nblck_prac)
                    imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/training%d_taskid%d_%s_stim%d.png',...
                        subj,iblck,blck.taskid,datestr(now,'yymmdd-HHMMSS'),itrl));
                else
                    imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/block%d_taskid%d_%s_stim%d.png',...
                        subj,iblck-nblck_prac,blck.taskid,datestr(now,'yymmdd-HHMMSS'),itrl));
                end
            end
            
            % keyboard input
            [response,tkey] = WaitKeyPress(keyresp,[],false);
            rslt.rt(itrl) = tkey-t;
            
            if response == 1 % left shape chosen
                is = blck.shape(il);
                rslt.resp(itrl) = il;
                rslt.respkb(itrl) = 1;
                rslt.shape(itrl) = is;
                if pedestal && (blck.condtn == 2)
                    Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
                end
                Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
                draw_stim(il,ir);
                Screen('FrameRect',video.h,color_frame,shape_box(:,1),8); % frame the shape chosen by participant
            else             % right shape chosen
                is = blck.shape(ir);
                rslt.resp(itrl) = ir;
                rslt.respkb(itrl) = 2;
                rslt.shape(itrl) = is;
                if pedestal && (blck.condtn == 2)
                    Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
                end
                Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
                draw_stim(il,ir);
                Screen('FrameRect',video.h,color_frame,shape_box(:,2),8); % frame the shape chosen by participant
            end
            if prob && (blck.condtn == 2)
                Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
            end
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            
            if make_scrshots % screenshot for response done
                imgscreen = Screen('GetImage',video.h);
                if (iblck<=nblck_prac)
                    imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/training%d_taskid%d_%s_resp%d.png',...
                        subj,iblck,blck.taskid,datestr(now,'yymmdd-HHMMSS'),itrl));
                else
                    imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/block%d_taskid%d_%s_resp%d.png',...
                        subj,iblck-nblck_prac,blck.taskid,datestr(now,'yymmdd-HHMMSS'),itrl));
                end
                
            end
            
            WaitSecs(.4);
            
            % time lapse before getting ready for outcome
            if keepshape
                draw_stim(il,ir);
                Screen('FrameRect',video.h,color_frame,shape_box(:,response),8);
            end
            if pedestal && (blck.condtn == 2)
                Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
            end
            if prob && (blck.condtn == 2)
                Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
            end
            Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            WaitSecs(.15);
            
            if itrl ~=(ntrl+1)
                % show response prob to focus on outcome
                if keepshape
                    draw_stim(il,ir);
                    Screen('FrameRect',video.h,color_frame,shape_box(:,response),8);
                end
                if pedestal && (blck.condtn == 2)
                    Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
                end
                if prob && (blck.condtn == 2)
                    Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
                    Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
                else
                    Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
                end
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
                if pedestal && (blck.condtn == 2)
                    Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
                end
                if prob && (blck.condtn == 2)
                    Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
                    Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
                else
                    Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
                end
                Screen('DrawingFinished',video.h);
                Screen('Flip',video.h);
                WaitSecs(.25);
                
                %% draw outcome
                if pedestal && (blck.condtn == 2)
                    Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
                end
                if keepshape
                    draw_stim(il,ir);
                    Screen('FrameRect',video.h,color_frame,shape_box(:,response),8);
                end
                
                if blck.condtn == 1
                    Screen('DrawTexture',video.h,fb_tex,[],fb_rec,[],[],[],colors(blck.outcome(response,itrl),:));
                else
                    if prob && (blck.condtn == 2)
                        Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
                    end
                    draw_handful(rangeDots(blck.outcome(response,itrl)),colors(1,:),colors(2,:));
                end
                
                Screen('DrawingFinished',video.h);
                Screen('Flip',video.h);
                
                if make_scrshots %&&(iblck<=nblck_prac) % screenshot for outcome
                    imgscreen = Screen('GetImage',video.h);
                    if (iblck<=nblck_prac)
                        imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/training%d_taskid%d_%s_outcome%d.png',...
                            subj,iblck,blck.taskid,datestr(now,'yymmdd-HHMMSS'),itrl));
                    else
                        imwrite(imgscreen,sprintf('./scrshots/Data_S%02d/block%d_taskid%d_%s_outcome%d.png',...
                            subj,iblck-nblck_prac,blck.taskid,datestr(now,'yymmdd-HHMMSS'),itrl));
                    end
                end
                
                WaitSecs(5*video.ifi);
                
                % hide outcome
                if pedestal && (blck.condtn == 2)
                    Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
                end
                if prob && (blck.condtn == 2)
                    Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
                    Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
                else
                    Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
                end
                Screen('DrawingFinished',video.h);
                Screen('Flip',video.h);
                WaitSecs(.2);
            
            if pedestal && (blck.condtn == 2)
                Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
            end
            if prob && (blck.condtn == 2)
                Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
            end
            Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h);
            
            end
            
        end % end of trial loop
        
        % store results and block substructures
        rslt.correct   = blck.correct(1:ntrl);
        rslt.perf      = sum(rslt.resp(1:ntrl)==blck.correct(1:ntrl))/ntrl;
        rslt.perfafter = sum(rslt.resp(2:ntrl+1)==blck.correct(1:ntrl))/ntrl;
        rslt.wsls_perf = run_model(blck, 1000, ntrl, calibration.setup.maxClass);
        
        if  perf_screen == true
            WaitSecs(.9);
            Screen('TextSize',video.h,round(txtsiz));
            labeltxt = sprintf('%d%% de bonnes réponses pour ce bloc, continuez!',round(rslt.perf*100));
            labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y/2);
            Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
            if rslt.perf>=rslt.wsls_perf
                labeltxt = 'en bonne voie pour le bonus';
                labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,3*video.y/4);
                Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
            end
            Screen('TextSize',video.h,round(txtsiz*.7));
            labeltxt = label_esc;
            labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y-round(1.2*ppd));
            Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h);
            WaitKeyPress(keywait);
        end
        
        expe(iblck).blck = blck;
        expe(iblck).rslt = rslt;
        expe_blck = expe(iblck); % necessary to create variable to save it
        
        % save temporary file (block per block)
        fpath = foldname;
        if (iblck<=nblck_prac)
            fname = sprintf('DOTCAT_S%02d_t%02d_%s',subj,iblck,datestr(now,'yyyymmdd-HHMM'));
        else
            fname = sprintf('DOTCAT_S%02d_b%02d_%s',subj,(iblck-nblck_prac),datestr(now,'yyyymmdd-HHMM'));
        end
        fname = fullfile(fpath,fname);
        save([fname,'.mat'],'expe_blck');
        
    end
    
    perf_tot = 0;
    perf_wsls_tot = 0;
    for i = (nblck_prac+1):nblck
        perf_tot = perf_tot+expe(i).rslt.perf/(nblck-nblck_prac);
        perf_wsls_tot = perf_wsls_tot+expe(i).rslt.wsls_perf/(nblck-nblck_prac);
    end
    
    fprintf('Win-stay-Loose-switch performance: %d%%\n', round(perf_wsls_tot*100))
    fprintf('Your performance: %d%%', round(perf_tot*100))
    if (perf_tot>=perf_wsls_tot)
        fprintf(' ... 5 euros de bonus!\n\n')
    end
    
    %% draw end screen
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
    if perf_tot>=perf_wsls_tot
        labeltxt3 = sprintf('%d%% de bonnes réponses au total, bonus de 5eur!',round(perf_tot*100));
        labelrec3 = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt3),video.x/2,5*video.y/6);
    else
        labeltxt3 = sprintf('%d%% de bonnes réponses au total, bien joué!',round(perf_tot*100));
        labelrec3 = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt3),video.x/2,5*video.y/6);
    end
    Screen('DrawText',video.h,labeltxt3,labelrec3(1),labelrec3(2),0);
    Screen('DrawingFinished',video.h);
    Screen('Flip',video.h,t+roundfp(1.500,0.500));
    WaitKeyPress(keywait);
    
    % save temporary file (whole expe in one .mat file)
    fpath = foldname;
    fname = sprintf('DOTCAT_S%02d_%s',subj,datestr(now,'yyyymmdd-HHMM'));
    fname = fullfile(fpath,fname);
    save([fname,'.mat'],'expe');
    
    if aborted
        Screen('CloseAll');
        return
    end
    
    % close Psychtoolbox
    Priority(0);
    Screen('CloseAll');
    FlushEvents;
    ListenChar(0);
    ShowCursor;
catch
    
    % close Psychtoolbox
    Priority(0);
    Screen('CloseAll');
    FlushEvents;
    ListenChar(0);
    ShowCursor;
    
    % handle error
    if nargout > 2
        errmsg = lasterror;
        errmsg = rmfield(errmsg,'stack');
    else
        rethrow(lasterror);
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
        is = blck.shape(il);
        Screen('DrawTexture',video.h,shape_tex(1,is),[],shape_rec(1,:),[],[],[],0);
        Screen('DrawTexture',video.h,shape_tex(2,is),[],shape_rec(1,:),[],[],[],color_shape);
        % draw right stimulus
        is = blck.shape(ir);
        Screen('DrawTexture',video.h,shape_tex(1,is),[],shape_rec(2,:),[],[],[],0);
        Screen('DrawTexture',video.h,shape_tex(2,is),[],shape_rec(2,:),[],[],[],color_shape);
    end
%%
    function draw_instr(iu,id)
        
        Screen('TextSize',video.h,round(txtsiz*.7));
        rec_esc = CenterRectOnPoint(Screen('TextBounds',video.h,label_esc),video.x/2,video.y-round(1.2*ppd));
        Screen('DrawText',video.h,label_esc,rec_esc(1),rec_esc(2),0);
        if iblck>nblck_prac
            Screen('TextSize',video.h,round(txtsiz));
            label_blck = sprintf('Bloc %d/8',iblck-nblck_prac);
            rec_blck = CenterRectOnPoint(Screen('TextBounds',video.h,label_blck),5*video.x/6,video.y-round(1.2*ppd));
            Screen('DrawText',video.h,label_blck,rec_blck(1),rec_blck(2),0);
        end
        
% % %         if iblck<=nblck_prac
% % %             rec_esc = CenterRectOnPoint(Screen('TextBounds',video.h,label_esc),video.x/2,video.y-round(1.2*ppd));
% % %             Screen('DrawText',video.h,label_esc,rec_esc(1),rec_esc(2),0);
% % %         else
% % %             rec_esc = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y-round(1.2*ppd));
% % %             Screen('DrawText',video.h,label_esc,rec_esc(1),rec_esc(2),0);
% % %             Screen('TextSize',video.h,round(txtsiz));
% % %             label_blck = sprintf('Bloc %d/8',iblck-nblck_prac);
% % %             rec_blck = CenterRectOnPoint(Screen('TextBounds',video.h,label_blck),5*video.x/6,video.y-round(1.2*ppd));
% % %             Screen('DrawText',video.h,label_blck,rec_blck(1),rec_blck(2),0);
% % %         end

        
        Screen('TextSize',video.h,round(txtsiz*instr_fac));
        if blck.taskid ==1 % Observer
            Screen('DrawText',video.h,label_obs,rec_obs(1),rec_obs(2),0);
            Screen('TextSize',video.h,txtsiz);
            Screen('DrawText',video.h,label_choice{1},choice_xy(1,1),choice_xy(1,2),0,[],0);
            Screen('DrawText',video.h,label_choice{2},choice_xy(2,1),choice_xy(2,2),0,[],0);
            Screen('DrawText',video.h,label_choice{1},choice_xy(3,1),choice_xy(3,2),0,[],0);
            Screen('DrawText',video.h,label_choice{2},choice_xy(4,1),choice_xy(4,2),0,[],0);
            draw_handful(k_instr,colors(iu,:),colors((3-iu),:),b_instr,d_instr,n_instr,true,bag_xy(2,1),bag_xy(2,2));
            draw_handful(k_instr,colors(id,:),colors((3-id),:),b_instr,d_instr,n_instr,true,bag_xy(3,1),bag_xy(3,2));
        else % Agent
            Screen('DrawText',video.h,label_agent,rec_agent(1),rec_agent(2),0);
            draw_handful(k_instr,colors(blck.epimap,:),colors((3-blck.epimap),:),b_instr,d_instr,n_instr,true,bag_xy(1,1),bag_xy(1,2));
        end
        
        % draw up stimulus
        is = blck.shape(iu);
        Screen('DrawTexture',video.h,shape_tex(1,is),[],shape_rec(3,:),[],[],[],0);
        Screen('DrawTexture',video.h,shape_tex(2,is),[],shape_rec(3,:),[],[],[],color_shape);
        
        % draw down stimulus
        is = blck.shape(id);
        Screen('DrawTexture',video.h,shape_tex(1,is),[],shape_rec(4,:),[],[],[],0);
        Screen('DrawTexture',video.h,shape_tex(2,is),[],shape_rec(4,:),[],[],[],color_shape);
    end

    function draw_handful(k,color1,color2,hsiz,dsiz,n,bag,xpos,ypos)
        % k dots from color1, nDots-k dots from color2
        % color1 and color2 are RGB triplets
        % hsiz;dsiz;n;bag;xpos;ypos are optional
        
        if nargin < 4
            hsiz = hand_siz;
            dsiz = dot_siz;
            n    = nDots;
            bag  = false;
            xpos = video.x/2;
            ypos = video.y/2;
        end
        
        color1 = reshape(color1,3,1);
        color2 = reshape(color2,3,1);
        
        dim = floor((floor(hsiz/dsiz)-3 )/2) ;
        [x, y] = meshgrid(-dim:1:dim, -dim:1:dim);
        
        pixelScale = hsiz / (dim * 2);
        x = x .* pixelScale;
        y = y .* pixelScale;
        
        dotPositionMatrix = [reshape(x, 1, numel(x)); reshape(y, 1, numel(x))];
        
        dotPositionMatrix(:,sqrt(dotPositionMatrix(1,:).^2+dotPositionMatrix(2,:).^2)>hsiz*.5) = nan;
        dotPositionMatrix = dotPositionMatrix(:,(sum(~isnan(dotPositionMatrix))==2));
        
        if size(dotPositionMatrix,2)<n
            error('too many dots for this handful size, please adapt')
        end
        
        dotPositionMatrix = dotPositionMatrix(:,randperm(size(dotPositionMatrix,2),n));
        
        % rotate handful
        theta = 4/12*pi*rand(1)+pi/12; % a + (b-a).*rand >> interval [pi/12 5*pi/12]
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        dotPositionMatrix = (dotPositionMatrix'*R)';
        dotPositionMatrix = dotPositionMatrix+rand(size(dotPositionMatrix))*6-3;
        
        % set colors
        dotColors = color2.*ones(3,n); % first set all dots from color2
        dotColors(:,randsample(n,k)) = color1.*ones(3,k); % set k dots to color1
        
        Screen('DrawDots', video.h, dotPositionMatrix,...
            dsiz, dotColors, [xpos,ypos], 2);
        
        if bag % for instructions screen, take colored bag contour
            bag_rec = CenterRectOnPoint(Screen('Rect',bag_tex(1)),xpos,ypos);
            if k>n/2
                Screen('DrawTexture',video.h,bag_tex,[],bag_rec,[],[],[],color1);
            else
                Screen('DrawTexture',video.h,bag_tex,[],bag_rec,[],[],[],color2);
            end
        end
    end

end