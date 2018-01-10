function [calibration,aborted,errmsg] = run_calibration(subj,xmin,xmax,ntrl)

addpath ./Toolboxes/Rand
addpath ./Toolboxes/IO
addpath ./Toolboxes/Stimuli/Visual

syncflip = true;

if nargin<1
    subjstr = '';
    xmin = .2; % minimal proportion shown for calibration
    xmax = .8; % maximal proportion shown for calibration
    ntrl = 120;
    foldname = sprintf('./Data/old_drafts');
    if ~exist(foldname,'dir')
        mkdir(foldname);
    end
elseif nargin == 1
    subjstr = sprintf('_S%02d',subj);
    xmin = .2; % cut off handfuls with proportion smaller
    xmax = .8; % cut off handfuls with proportion bigger
    ntrl = 120;
    % create data folder for interim subject files
    foldname = sprintf('./Data/S%02d',subj);
    if ~exist(foldname,'dir')
        mkdir(foldname);
    end
elseif nargin == 4
    subjstr = sprintf('_S%02d',subj);
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

%%
make_scrshots = false; % make screenshots?
if make_scrshots
    warning('Will make screenshots in ./scrshots!');
end

% set screen parameters
iscr = 0; % screen index
res  = []; % screen resolution
fps  = []; % screen refresh rate
ppd  = 40; % number of screen pixels per degree of visual angle

% set stimulation parameters
lumibg    = 128/255;        % background luminance (grey)
color_pedestal = 120/255;   % pedestal background luminance
fixtn_siz = 0.2*ppd;        % fixation point size
dot_siz   = 0.4*ppd;        % dot size
hand_siz  = 6.5*ppd;        % handful size
probwdth  = 1;              % handful probe contour width
color_off = 7*ppd;          % color offset for instruction screen
info_fac  = 2.5;            % informational text magnification factor

% set list of color-wise R/G/B values
color_rgb = [ ...
    % http://www.workwithcolor.com/hsl-color-picker-01.htm?
    136,236,136; ... % green L 73%, Lum 80%, 120°
    136,220,236; ... % blue  L 73%, Lum 80%, 190°
    % from colorbrewer2.org
    127,201,127; ... % green 4 > lum 69% (L 64)
    190,174,212; ... % violet 5 > lum 73% >> 69% = 180, 161, 206 (L 72)
    253,192,134; ... % orange
    %
    102,194,165; ... % vert poubelle 7
    252,141,98;  ... % vermillon
    141,160,203; ... % mauve 9
    %
    179,205,227; ... % bleu clair
    222,203,228; ... % mauve
    251,180,174; ... % rouge pale
    % foncé NON
    55,126,184;  ... % bleu fonçé
    152,78,163;  ... % violet
    % http://www.workwithcolor.com/hsl-color-picker-01.htm?
    %154,213,154; ... % green L 72%, Lum 76%, 120°
    %154,203,213; ... % blue  L 72%, Lum 76%, 190°
    ]/255;

% set color opacity
color_opa = 2/3; % color opacity (close to 1 no opacity)
color_rgb = color_rgb*color_opa+lumibg*(1-color_opa); % comment?

% set params stim presentation
%ntrl  = 4;
sig   = 5; % adapt?
nblck = 2;
nDots = 35;
tmin  = 0.2; % minimal correct classification
tmax  = 0.8; % maximal correct classification

% customizable
prob     = false; % show handful probe contour?
respprob = true;  % show response probe that disappear just before stim?
pedestal = true;  % show handful in a pedestal?

colors = [1 2]; % color 1: green; color 2: blue
colorstxt = {'green' 'blue'};

resp = zeros(nblck,ntrl); % 1: green or 2: blue
rt   = zeros(nblck,ntrl); % response time
calibration = struct();
calibration.expe = [];

calibration.setup.ppd       = ppd;
calibration.setup.lumibg    = lumibg;
calibration.setup.color_pedestal = color_pedestal;
calibration.setup.fixtn_siz = fixtn_siz;
calibration.setup.dot_siz   = dot_siz;
calibration.setup.hand_siz  = hand_siz;
calibration.setup.probwdth  = probwdth;
calibration.setup.colors    = color_rgb(colors,:);
calibration.setup.nDots     = nDots;
calibration.setup.minClass  = tmin;
calibration.setup.maxClass  = tmax;
calibration.setup.minCalib  = xmin;
calibration.setup.maxCalib  = xmax;
calibration.setup.ntrl      = ntrl;

pres = 0;
while length(pres)<ntrl*nblck
    pres = round(randn(1,(nblck+1)*ntrl)*sig+nDots/2);
    pres = pres(pres~=(nDots/2));
    pres = pres(pres<=round(xmax*nDots));
    if round(.2*nDots)==0
        pres = pres(pres>=1);
    else
        pres = pres(pres>=round(xmin*nDots));
    end
end

pres = pres(1:nblck*ntrl);
calibration.expe.pres = pres; % contains number of green dots used during calibration
calibration.expe.presfrac = pres/nDots; % contains fraction of green dots used during calibration

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
    
    %% create textures

    % create fixation point
    img = CreateCircularAperture(fixtn_siz);
    fixtn_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    fixtn_rec = CenterRectOnPoint(Screen('Rect',fixtn_tex(1)),video.x/2,video.y/2);
    
    % create handful probe contour
    img = CreateCircle(hand_siz+1.5*dot_siz,probwdth);
    hand_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    hand_rec = CenterRectOnPoint(Screen('Rect',hand_tex(1)),video.x/2,video.y/2);
    
    % create pedestal
    img = CreateCircularAperture(hand_siz+1.5*dot_siz);
    ped_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    ped_rec = CenterRectOnPoint(Screen('Rect',ped_tex(1)),video.x/2,video.y/2);
    
    % create response probe contour
    img = CreateCircle(fixtn_siz+6*probwdth,probwdth);
    resp_tex = Screen('MakeTexture',video.h,cat(3,ones(size(img)),img),[],[],2);
    resp_rec = CenterRectOnPoint(Screen('Rect',resp_tex(1)),video.x/2,video.y/2);
    
    % create instruction dots
    img = CreateCircularAperture(dot_siz*5);
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
    
    for iblck = 1:nblck
        
        % draw instruction screen
        Screen('TextSize',video.h,round(txtsiz*info_fac));
        labeltxt = sprintf('La majorité des billes sont');
        labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y/7);
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        Screen('DrawTexture',video.h,instr_tex,[],instr_rec(1,:),[],[],[],color_rgb(colors(iblck),:));
        Screen('DrawTexture',video.h,instr_tex,[],instr_rec(2,:),[],[],[],color_rgb(colors(3-iblck),:));
        Screen('TextSize',video.h,round(txtsiz*.7));
        labeltxt = sprintf('appuyez sur [espace] pour continuer');
        labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y-round(1.2*ppd));
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        Screen('DrawingFinished',video.h);
        Screen('Flip',video.h,t+roundfp(1.500,0.500));
        WaitKeyPress(keywait);
        
        if pedestal
            Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
        end
        Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
        if prob % draw handful probe contour
            Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
        end
        Screen('DrawingFinished',video.h);
        Screen('Flip',video.h);
        WaitSecs(25*video.ifi);
        
        pres = calibration.expe.pres((iblck-1)*ntrl+(1:ntrl));
        
        for itrl = 1:length(pres)
            
            %check if abort key is pressed
            if CheckKeyPress(keyquit)
                Screen('CloseAll')
                aborted = true;
                ShowCursor;
                sca;
                break
            end
                       
            % draw fixation point before stim
            if pedestal
                Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
            end
            Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
            if prob % draw handful probe contour
                Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
            end
            if respprob % draw response probe
                Screen('DrawTexture',video.h,resp_tex,[],resp_rec,[],[],[],0);
            end
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            WaitSecs(20*video.ifi);
            
            % response probe disappears
            if pedestal
                Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
            end
            Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
            if prob % draw handful probe contour
                Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
            end
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            WaitSecs(15*video.ifi);
            
            % draw stim
            if pedestal
                Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
            end
            draw_handful(pres(itrl),calibration.setup.colors(1,:),calibration.setup.colors(2,:));
            if prob % draw handful probe contour
                Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
            end
            Screen('DrawingFinished',video.h);
            t = Screen('Flip',video.h);
               
            if make_scrshots&&(itrl<=10)&&(iblck ==1)
                imgscreen = Screen('GetImage',video.h);
                imwrite(imgscreen,sprintf('./scrshots/calibration%d.png',itrl));
            end
            WaitSecs(5*video.ifi); % duration of stim presentation
            %WaitKeyPress(keywait);
            
            % disparition of stim awaiting for response
            if pedestal
                Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
            end
            Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
            if prob
                Screen('DrawTexture',video.h,hand_tex,[],hand_rec,[],[],[],0);
            end
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h);
            
            % keyboard input
            [response,tkey] = WaitKeyPress(keyresp,[],false);
            rt(iblck,itrl)  = tkey-t;
            
            if response == 1 
                resp(iblck,itrl) = iblck; % left color majority
            elseif response == 2
                resp(iblck,itrl) = 3-iblck; % right color majority
            end
            
            WaitSecs(40*video.ifi);
                       
        end % end of trial loop      
        if pedestal
            Screen('DrawTexture',video.h,ped_tex,[],ped_rec,[],[],[],color_pedestal);
        end
        Screen('DrawTexture',video.h,fixtn_tex,[],fixtn_rec,[],[],[],0);
        Screen('DrawingFinished',video.h);
        Screen('Flip',video.h);
        WaitSecs(.5);
    end
    
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

    calibration.expe.resp = reshape(resp',1,nblck*ntrl);
    for i = 1:nDots
        calibration.expe.titr(i) = sum(calibration.expe.resp(calibration.expe.pres==i)==1)/sum(calibration.expe.pres==i); % green sigmoid
    end
    calibration.expe.rt = reshape(rt',1,nblck*ntrl);
    
    [b,~,stat] = glmfit(calibration.expe.pres'/nDots,calibration.expe.resp' == 1,'binomial','link','probit');
    
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
    
    figure('Visible','Off');
    subplot(2,1,1)
    histogram(calibration.expe.pres/nDots,20);
    xlim([0 1])
    subplot(2,1,2)
    plot((1:nDots)/nDots,calibration.expe.titr);
    xlim([0 1])
    hold on
    plot(0:0.01:1,normcdf(b(1)+(0:0.01:1)*b(2)));
    xlabel(sprintf('%s Dots (%%)',colorstxt{1}));
    ylabel(sprintf('classified as %s (%%)',colorstxt{1}));
    
    calibration.setup.range   = [fzero(@(xx)normcdf(b(1)+xx*b(2))-tmin,[0,1]) ...
                                 fzero(@(xx)normcdf(b(1)+xx*b(2))-tmax,[0,1])];
    calibration.setup.minDots = round(calibration.setup.range(1)*nDots);
    calibration.setup.maxDots = round(calibration.setup.range(2)*nDots);
    calibration.expe.b        = b;
    calibration.expe.stat     = stat;
    
    fpath = foldname;
    fname = fullfile(fpath,sprintf('DOTCAT_calibration%s_%s',subjstr,datestr(now,'yyyymmdd-HHMM')));
    save([fname,'.mat'],'calibration');
    savefig(fname);
    
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

        dim = floor((floor(hsiz/dsiz)-3)/2) ;
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
        % first set all dots from color2
        dotColors = color2.*ones(3,n);
        % set k dots to color1
        dotColors(:,randsample(n,k)) = color1.*ones(3,k); 
        
        Screen('DrawDots', video.h, dotPositionMatrix,...
            dsiz, dotColors, [xpos,ypos], 2);
        if bag
            % create handful probe contour
            imbag = CreateCircle(hsiz+1.5*dsiz,probwdth);
            bag_tex = Screen('MakeTexture',video.h,cat(3,ones(size(imbag)),imbag),[],[],2);
            bag_rec = CenterRectOnPoint(Screen('Rect',bag_tex(1)),xpos,ypos);
            if k>n/2
                Screen('DrawTexture',video.h,bag_tex,[],bag_rec,[],[],[],color1);
            else
                Screen('DrawTexture',video.h,bag_tex,[],bag_rec,[],[],[],color2);
            end
        end
    end

    function [p] = deg2pix(d,b)
        % convert degrees of visual angle into screen pixels
        % where d - degrees of visual angle
        %       b - rounding factor (default: none)
        p = d*ppd;
        if nargin > 1 && b > 0
            p = max(round(p/b),1)*b;
        end
    end
end
