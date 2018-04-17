function [calibration,aborted,errmsg] = run_calibration(subj,calibration,trainornot)

addpath ./Toolboxes/Rand
addpath ./Toolboxes/IO
addpath ./Toolboxes/Stimuli/Visual
addpath ./Toolboxes/Stimuli/Auditory


%% customizable options
syncflip  = false;
plotornot = false;
respprob  = true;   % show response probe that disappear just before stim?
make_scrshots = false; % make screenshots?
if make_scrshots
    warning('Will make screenshots in ./scrshots!');
end

%%
if nargin<1
    subjstr = '';
    ntrl = 20;
    nblck = 2;
    foldname = sprintf('./Data/old_drafts');
    if ~exist(foldname,'dir')
        mkdir(foldname);
    end
    calibration = gen_calibration(ntrl,nblck);
    trainornot = false;
elseif nargin == 1
    subjstr = sprintf('_S%02d',subj);
    ntrl  = 120;
    nblck = 2;
    % create data folder for interim subject files
    foldname = sprintf('./Data/S%02d',subj);
    if ~exist(foldname,'dir')
        mkdir(foldname);
    end
    fname = fullfile(foldname,sprintf('DOTCAT_calibration_raw_S%02d_%s',subj,datestr(now,'yyyymmdd-HHMM')));
    calibration = gen_calibration(ntrl,nblck,fname);  
    trainornot = false;
elseif nargin == 2
    subjstr = sprintf('_S%02d',subj);
    ntrl = calibration.cfg.ntrl;
    nblck = calibration.cfg.nblck;
    % create data folder for interim subject files
    foldname = sprintf('./Data/S%02d',subj);
    if ~exist(foldname,'dir')
        mkdir(foldname);
    end
    trainornot = false;
elseif nargin == 3
    subjstr = sprintf('_S%02d',subj);
    ntrl = calibration.cfg.ntrl;
    nblck = calibration.cfg.nblck;
    % create data folder for interim subject files
    foldname = sprintf('./Data/S%02d',subj);
    if ~exist(foldname,'dir')
        mkdir(foldname);
    end
end


%%
% define output arguments
aborted = false; % aborted prematurely?
errmsg  = []; % error message


% set screen parameters
iscr = 2;%2  % screen index
res  = []; % screen resolution
fps  = []; % screen refresh rate
ppd  = 40; % number of screen pixels per degree of visual angle

% set stimulation parameters
fix_siz   = 0.2*ppd; % fixation point size
dot_siz   = calibration.cfg.dmtr; % DOT size
probwdth  = 1;	     % response probe contour width
color_off = 7*ppd;   % color offset for instruction screen
info_fac  = 2.5;     % informational text magnification factor

% set colors
lumibg    = 128/255; % background luminance (grey)
% http://www.workwithcolor.com/hsl-color-picker-01.htm? | colorbrewer2.org
color_rgb = [ ...
%     178,223,138; ...
%     166,206,227; ...
    190,190,190; ...
    66, 66, 66 ; ...
    136,236,136; ... % green L 73%, Lum 80%, 120°
    136,220,236; ... % blue  L 73%, Lum 80%, 190°
    255,170,128; ... % orange L:75%; Lum:75%
    128,191,255; ... % bleu   L:75%; Lum:75%
    ]/255;

%set color opacity > depend on screen properties
% color_opa = 2/3; % color opacity (close to 1 no opacity)
% color_rgb = color_rgb*color_opa+lumibg*(1-color_opa);

% add screen parameter to cfg structure
calibration.cfg.ppd      = ppd;
calibration.cfg.lumibg   = lumibg;
calibration.cfg.colors   = color_rgb(1:2,:);
calibration.cfg.fix_siz  = fix_siz;
calibration.cfg.probwdth = probwdth;

pat = calibration.stim.pattern;
color_prop = calibration.stim.color_prop;

calibration.rslt = [];

if trainornot
    train_calib = importdata('./Raw/DOTCAT_training_calibration_1.mat');
    train_idx = Shuffle(1:length(train_calib.stim.pattern));
    train_pat = train_calib.stim.pattern(train_idx);
    train_cprop = train_calib.stim.color_prop(train_idx);
end


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
    
    
    %% load win/lose sounds
    soundsnames = {'bip_win.mp3', 'bip_early.wav'};
    bips = cell(1,2);
    for isound = 1:2
        bips{isound} = audioread(sprintf('./Bips/%s', soundsnames{isound}))';
        if size(bips{isound}, 1) < 2
            bips{isound} = repmat(bips{isound}, [2 1]);
        end
    end
    
    InitializePsychSound(1);
    audio.freq = 48000;
    audio.h= PsychPortAudio('Open', [], 1, 1, 48000, 2);
    PsychPortAudio('RunMode',audio.h,1);
    PsychPortAudio('Volume',audio.h,0.7);
    
    [tonebuf,tonepos] = CreateAudioBuffer( ...
        bips{1}, ...
        bips{2});
    PsychPortAudio('FillBuffer',audio.h,tonebuf);
    
    %% create textures
    
    % create fixation point
    img = CreateCircularAperture(fix_siz);
    fixtn_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    fixtn_rec = CenterRectOnPoint(Screen('Rect',fixtn_tex(1)),video.x/2,video.y/2);
    
    % create response probe contour
    img = CreateCircle(fix_siz+6*probwdth,probwdth);
    resp_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    resp_rec = CenterRectOnPoint(Screen('Rect',resp_tex(1)),video.x/2,video.y/2);
    
    % create instruction dots
    img = CreateCircularAperture(dot_siz);
    instr_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    instr_rec(1,:) = CenterRectOnPoint(Screen('Rect',instr_tex(1)),video.x/2-color_off,video.y/2);
    instr_rec(2,:) = CenterRectOnPoint(Screen('Rect',instr_tex(1)),video.x/2+color_off,video.y/2);
    
    % first flip
    t = Screen('Flip',video.h);
    
    Screen('TextStyle',video.h,0);
    Screen('TextSize',video.h,round(txtsiz*.7));
    labeltxt = sprintf('appuyez sur [espace] pour démarrer');
    labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y-round(1.2*ppd));
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
    
    if trainornot
        
        % draw instruction screen
        draw_instr(2); % >> left = color2
        Screen('Flip',video.h,t+roundfp(1.500,0.500));
        WaitKeyPress(keywait);
        
        Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
        Screen('DrawingFinished',video.h);
        Screen('Flip',video.h);
        WaitSecs(25*video.ifi);
        
        % left = color2 in traiing
        train_correct = (train_cprop>.5)+1;
        
        for itrl = 1:train_calib.cfg.ntrl
            
            %check if abort key is pressed, nothing to save in training, quit
            if CheckKeyPress(keyquit)
                PsychPortAudio('Close');     
                aborted = true;
                Screen('CloseAll');
                FlushEvents;
                ListenChar(0);
                ShowCursor;
                sca;
                return
            end
            
            % draw fixation point before stim
            Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
            if respprob % draw response probe
                Screen('DrawTexture',video.h,resp_tex,[],resp_rec,[],[],[],0);
            end
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            WaitSecs(20*video.ifi);
            
            % response probe disappears
            Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            WaitSecs(15*video.ifi);
            
            % draw stim
            draw_stim(train_pat(itrl));
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            WaitSecs(5*video.ifi); % duration of stim presentation
            
            % disparition of stim awaiting for response
            Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            
            % keyboard input
            [response,tkey] = WaitKeyPress(keyresp,[],false);
            
           
            if train_correct(itrl) == response % win
                PsychPortAudio('SetLoop',audio.h,tonepos(1,1),tonepos(1,2));
                PsychPortAudio('Start',audio.h,1,tkey+.2,0);
            else % loose
                PsychPortAudio('SetLoop',audio.h,tonepos(2,1),tonepos(2,2));
                PsychPortAudio('Start',audio.h,1,tkey+.2,0);
            end
                    
            WaitSecs(40*video.ifi);
                        
        end % end of trial loop
        
        Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
        Screen('DrawingFinished',video.h);
        Screen('Flip',video.h);
        WaitSecs(.5);
        
        % draw end of training screen
        Screen('TextSize',video.h,round(txtsiz*info_fac));
        labeltxt = sprintf('fin de l''entraînement!');
        labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y/2);
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        Screen('TextSize',video.h,round(txtsiz*.7));
        labeltxt = sprintf('appuyez sur [espace] pour continuer');
        labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y-round(1.2*ppd));
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        Screen('DrawingFinished',video.h);
        Screen('Flip',video.h);
        WaitKeyPress(keywait);
    end
    
    for iblck = 1:nblck
        
        % draw instruction screen
        draw_instr(iblck)
        Screen('Flip',video.h,t+roundfp(1.500,0.500));
        WaitKeyPress(keywait);
        
        Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
        Screen('DrawingFinished',video.h);
        Screen('Flip',video.h);
        WaitSecs(25*video.ifi);
        
        resp   = zeros(1,ntrl);
        respkb = zeros(1,ntrl);
        rt     = zeros(1,ntrl);
        
        % left = color(iblck)
        if mod(iblck,2)
            correct = (color_prop(iblck,:)<.5)+1;
        else
            correct = (color_prop(iblck,:)>.5)+1;
        end
        
        for itrl = 1:ntrl
            
            %check if abort key is pressed
            if CheckKeyPress(keyquit)
                aborted = true;
                break
            end
            
            % draw fixation point before stim
            Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
            if respprob % draw response probe
                Screen('DrawTexture',video.h,resp_tex,[],resp_rec,[],[],[],0);
            end
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            WaitSecs(20*video.ifi);
            
            % response probe disappears
            Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            WaitSecs(15*video.ifi);
            
            % draw stim
            draw_stim(pat(iblck,itrl));
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            WaitSecs(5*video.ifi); % duration of stim presentation
            
            % disparition of stim awaiting for response
            Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h);
            
            % keyboard input
            [response,tkey] = WaitKeyPress(keyresp,[],false);
            rt(itrl)  = tkey-t;
            
            % play sound
            if correct(itrl) == response % win
                PsychPortAudio('SetLoop',audio.h,tonepos(1,1),tonepos(1,2));
                PsychPortAudio('Start',audio.h,1,tkey+.2,0);
            else % loose
                PsychPortAudio('SetLoop',audio.h,tonepos(2,1),tonepos(2,2));
                PsychPortAudio('Start',audio.h,1,tkey+.2,0);
            end
            
            if response == 1
                resp(itrl)   = iblck;   % left color majority
            elseif response == 2
                resp(itrl)   = 3-iblck; % right color majority
            end
            respkb(itrl) = response;
            
            WaitSecs(40*video.ifi);
                        
        end % end of trial loop
        
        Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
        Screen('DrawingFinished',video.h);
        Screen('Flip',video.h);
        WaitSecs(.5);
                    
        calibration.rslt.resp(iblck,:)   = resp;
        calibration.rslt.respkb(iblck,:) = respkb;
        calibration.rslt.rt(iblck,:)     = rt;
        
        if aborted
            fprintf('Calibration aborted!\n')
            break
        end
        
        fprintf('B%d: Correct categorization color 1: %d%%\n',iblck,...
                round(100*sum(resp(color_prop(iblck,:)>.5)==1)/sum(color_prop(iblck,:)>.5)));
        fprintf('B%d: Correct categorization color 2: %d%%\n',iblck,...
                round(100*sum(resp(color_prop(iblck,:)<.5)==2)/sum(color_prop(iblck,:)<.5)));
    end
    
    if ~aborted
        Screen('TextSize',video.h,round(txtsiz*info_fac));
        labeltxt = sprintf('Merci!');
        labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y/2);
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        Screen('TextSize',video.h,round(txtsiz*.7));
        labeltxt = sprintf('appuyez sur [espace] pour quitter');
        labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y-round(1.2*ppd));
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        Screen('DrawingFinished',video.h);
        Screen('Flip',video.h);
        WaitKeyPress(keywait);
        
        fprintf('Correct categorization color 1: %d%%\n',round(100*sum(calibration.rslt.resp(color_prop>.5)==1)/sum(color_prop(:)>.5)));
        fprintf('Correct categorization color 2: %d%%\n',round(100*sum(calibration.rslt.resp(color_prop<.5)==2)/sum(color_prop(:)<.5)));
        
        try
            % fit calibration results
            calibration = fit_calib(calibration); % add rtmax argument if necessary
        catch
            fprintf('No fit for psychometric curve found!\n')
        end
    end
    
    % close Psychtoolbox
    Priority(0);
    PsychPortAudio('Stop',audio.h);
    PsychPortAudio('Close');
    Screen('CloseAll');
    FlushEvents;
    ListenChar(0);
    ShowCursor;   
    
    hdr.subj = subj;
    hdr.date = datestr(now,'yyyymmdd-HHMM');
    
    calibration.hdr = hdr;
    try
        calibration = orderfields(calibration,{'hdr','cfg','stim','rslt'});
    catch
        calibration = orderfields(calibration,{'hdr','cfg','stim'});
    end
    
    fpath = foldname;
    fname = sprintf('DOTCAT_calibration%s_%s',subjstr,datestr(now,'yyyymmdd-HHMM'));
    if aborted
        if ~exist([fpath,'/aborted'],'dir')
            mkdir([fpath,'/aborted']);
        end
        fpath = [foldname,'/aborted/'];
        fname = [fname,'_aborted'];
    end
    save([fullfile(fpath,fname),'.mat'],'calibration');
    if plotornot
        savefig(fname);
    end
    
catch ME
    % close Psychtoolbox
    Priority(0);
    Screen('CloseAll');
    FlushEvents;
    ListenChar(0);
    ShowCursor;
    % close audio
    if exist('audio','var') && ~isempty(audio)
        PsychPortAudio('Close');
    end
    
    % handle error
    if nargout > 2
        errmsg = ME;
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

    function draw_stim(pat)
        for ipoly = 1: numel(pat.poly_list)
            Screen('FillPoly', video.h,color_rgb(pat.poly_list(ipoly).color_index,:),...
                [pat.poly_list(ipoly).vertex_list(:,1) pat.poly_list(ipoly).vertex_list(:,2)]+[video.x/2 video.y/2],1);
        end
    end

    function draw_instr(iblck)
        % draw instruction screen
        Screen('TextSize',video.h,round(txtsiz*info_fac));
        labeltxt = sprintf('Je range cette bille dans le sac');
        labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y/7);
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        Screen('DrawTexture',video.h,instr_tex,[],instr_rec(1,:),[],[],[],color_rgb(iblck,:));
        Screen('DrawTexture',video.h,instr_tex,[],instr_rec(2,:),[],[],[],color_rgb((3-iblck),:));
        Screen('TextSize',video.h,round(txtsiz*.7));
        labeltxt = sprintf('appuyez sur [espace] pour continuer');
        labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y-round(1.2*ppd));
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        Screen('DrawingFinished',video.h);
    end
end
