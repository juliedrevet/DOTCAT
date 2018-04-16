function show_pattern(pat)


addpath ./Toolboxes/Rand
addpath ./Toolboxes/IO
addpath ./Toolboxes/Stimuli/Visual

% customizable options
make_scrshots = true;
draw_prop = true;
dmtr = pat.cfg.dmtr;

% set background color
lumibg    = 128/255;
% set list of color-wise R/G/B values
% http://www.workwithcolor.com/hsl-color-picker-01.htm?
color_rgb = [ ...
    190,190,190; ...
    66, 66, 66 ; ...
    % http://www.workwithcolor.com/hsl-color-picker-01.htm? | colorbrewer2.org
    103,228,103;...
    103,207,228;...
    255,170,128; ... % orange L:75%; Lum:75%
    128,191,255; ... % bleu   L:75%; Lum:75%
    255,153,102; ... % orange L:70%; Lum:71%
    102,179,255; ... % bleu   L:70%; Lum:71%
    255,179,102; ... % orange L:70%; Lum:77%
    102,179,255; ... % bleu   L:70%; Lum:71%
    255,160,65 ; ... % orange AGOBS L:63%; Lum:72%
    65 ,160,255; ... % blue   AGOBS L:63%; Lum:65%
    136,236,136; ... % green L 73%, Lum 80%, 120°
    255,000,000; ... % red
    136,220,236; ... % blue  L 73%, Lum 80%, 190°
    127,201,127; ... % green 4 > lum 69% (L 64)
    190,174,212; ... % violet 5 > lum 73% >> 69% = 180, 161, 206 (L 72)
    253,192,134; ... % orange
    102,194,165; ... % vert poubelle 7
    252,141,98;  ... % vermillon
    141,160,203; ... % mauve 9
    179,205,227; ... % bleu clair
    222,203,228; ... % mauve
    251,180,174; ... % rouge pale
    55,126,184;  ... % bleu foncé
    152,78,163;  ... % violet
    154,213,154; ... % green L 72%, Lum 76%, 120°
    154,203,213; ... % blue  L 72%, Lum 76%, 190°
    ]/255;

% set color opacity
% color_opa = 2/3; % color opacity
% color_rgb = color_rgb*color_opa+lumibg*(1-color_opa);

% set screen parameters
iscr = 0; % screen index

PsychDefaultSetup(2);
video.i = iscr;
video.res = Screen('Resolution',video.i);
video.h = PsychImaging('OpenWindow',video.i,0);
[video.x,video.y] = Screen('WindowSize',video.h);
ifi = Screen('GetFlipInterval', video.h);
Screen('BlendFunction', video.h, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

KbName('UnifyKeyNames');
keywait = KbName('DownArrow');

if nargin==2
    img_rec = CenterRectOnPoint([0 0 dmtr dmtr],video.x/2,video.y/2);
    scr_rec = CenterRectOnPoint([0 0 dmtr*2 dmtr*2],video.x,video.y);
end

% Over-write the screen in grey
Screen('FillRect', video.h, lumibg);

for ipat = 1:length(pat)
    for i = 1: numel(pat(ipat).poly_list)
        Screen('FillPoly', video.h,color_rgb(pat(ipat).poly_list(i).color_index,:),...
            [pat(ipat).poly_list(i).vertex_list(:,1) pat(ipat).poly_list(i).vertex_list(:,2)]+[video.x/2 video.y/2],1);
    end
%     labeltxt = sprintf('%02d',round(pat(ipat).color_prop_true*100));
%     if nargin == 2
%         Screen('DrawText',video.h,labeltxt,img_rec(1)+2,img_rec(2)+2,0);
%     else
%         Screen('DrawText',video.h,labeltxt,video.x/4,video.y/4,0);
%     end
    Screen('DrawingFinished',video.h);
    Screen('Flip',video.h);
    
    if make_scrshots
        if nargin==2
            imgscreen = Screen('GetImage',video.h,scr_rec);
        else
            imgscreen = Screen('GetImage',video.h);
        end
        imwrite(imgscreen,sprintf('./scrshots/fraction%d_%s.png',round(pat(ipat).color_prop_true*100),datestr(now,'yyyymmdd-HHMMSS')));
    end
    
    WaitKeyPress(keywait);
end

KbStrokeWait;
sca;
end