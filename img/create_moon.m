%% create moon contour
lune_raw = imread('./src/moon.png');
if size(lune_raw,3)>1
    lune_raw = lune_raw(:,:,1);
end

lune = 255-lune_raw;
elune_inside = edge(lune,'canny');
lune_inside = imfill(elune_inside,'holes');
lune_inside = lune_inside-elune_inside;
lune_inside = [zeros(size(lune_inside,1),3) lune_inside zeros(size(lune_inside,1),3)];
lune_inside = [zeros(3,size(lune_inside,2));lune_inside; zeros(3,size(lune_inside,2))];

lune = [zeros(size(lune,1),3) lune zeros(size(lune,1),3)];
lune = [zeros(3,size(lune,2));lune; zeros(3,size(lune,2))];
elune_contour = edge(lune,'canny');

lune_contour = imfill(elune_contour, 'holes');

act_height = size(imread('shape1.png'),1);
act_width = size(imread('shape1.png'),2);

leftside = floor((act_width-size(lune_contour,2))/2);
upperside = floor((act_height-size(lune_contour,1))/2);
shape8c = uint8(zeros(size(imread('shape1.png'))));
shape8c(upperside+1:upperside+size(lune_contour,1),leftside+1:leftside+size(lune_contour,2)) = lune_contour*255;
shape8 = uint8(zeros(size(imread('shape1.png'))));
shape8(upperside+1:upperside+size(lune_inside,1),leftside+1:leftside+size(lune_inside,2)) = lune_inside*255;

figure()
imshow(shape8c)
imwrite(shape8c,'shape8c.png','PNG')
figure()
imshow(shape8)
imwrite(shape8,'shape8.png','PNG')
figure()
imshow(shape8c-shape8)