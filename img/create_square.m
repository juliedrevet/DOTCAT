%% create square

%length fixed
large = 360;
start = floor((size(imread('shape1c.png'),2)-large)/2);
start_in = start+31;
last = size(imread('shape1c.png'),2)-floor((size(imread('shape1c.png'),2)-large)/2);
if large==last-start
    last = last+1;
end
last_in = last-31;

% % % % length = circle diameter = 417
circle = imread('shape1c.png');
% % % circle_in = imread('shape1.png');
% % % start = find(circle(size(circle,1)/2,:),1,'first');
% % % start_in = find(circle_in(size(circle_in,1)/2,:),1,'first');
% % % last = find(circle(size(circle,1)/2,:),1,'last');
% % % last_in = find(circle_in(size(circle_in,1)/2,:),1,'last');
% % % large = last-start+1;

% get the corners' structure of the cross (shape6)
im6 = imread('shape6c.png');
im = im6;

for k = 1:150
    idx = find(im6(:,k));
    if any(idx)
        firstidx(k) = idx(1);
        lastidx(k) = idx(end);
    else
        firstidx(k) = 0;
        lastidx(k) = 0;
    end
end

first_line = min(nonzeros(firstidx));
first_col = find(firstidx,1);
corner_length = length(nonzeros(firstidx));

corner1 = im6(first_line:(first_line+corner_length-1),first_col:(first_col+corner_length-1));
corner2 = corner1(flip(1:corner_length),:);
corner3 = corner1(:,flip(1:corner_length));
corner4 = corner1(flip(1:corner_length),flip(1:corner_length));

square = uint8(zeros(size(circle)));
square(start:start+corner_length-1,start:start+corner_length-1) = corner1;
square(last-corner_length+1:last,start:start+corner_length-1) = corner2;
square(start:start+corner_length-1,last-corner_length+1:last) = corner3;
square(last-corner_length+1:last,last-corner_length+1:last) = corner4;

f1 = find(square(1:300,start),1,'last');
f2 = find(square(301:600,start),1,'first');
square(f1:301+f2,start) = 4;
square(start,f1:301+f2) = 4;
square(f1:301+f2,last) = 4;
square(last,f1:301+f2) = 4;
square(f1+1:300+f2,start+1:last-1)=255;
square(start+1:last-1,f1+1:300+f2)=255;

square_in = uint8(zeros(size(circle)));
square_in(start_in:last_in,start_in:last_in) = 255;

imwrite(square,'shape7c.png','PNG')
imwrite(square_in,'shape7.png','PNG')

imshow(square-square_in)