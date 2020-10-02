clear all;
clc;

 HM = imread('C:\Mtech\Mtech thesis\matlab\Implementation\H&MAoutput\newoutputs\26_training_H&MA.tif');
% folder2 = 'C:\Mtech\Mtech thesis\matlab\Implementation\H&MAoutput\newoutputs';
% ds2 = imageDatastore(folder2);
%while hasdata(ds2) 
    
    %HM = read(ds2); 
    HM = im2uint8(HM(:,:,[1 1 1]));
    white_pixel = sum(HM(:));
    imshow(HM);

    % [centers, radii, metric] = imfindcircles(origImg,[1 10])
    % viscircles(centers, radii,'EdgeColor','b');



   number_of_circles = 0;
    %.... you have a loop going on ...
      %.... you find some circles, and their coordinates are in X and Y and size in R)
      %draw the circles
    [centers, radii, metric] = imfindcircles(HM,[1 10]);
    viscircles(centers, radii,'EdgeColor','b');
    number_of_circles = number_of_circles + length(centers); 
   %.... end of loop

%end