
clear all;
clc;

origImg = imread('C:\Mtech\Mtech thesis\matlab\Implementation\H&MAoutput\finaloutput\32_training_H&MA.tif');
origImg = im2uint8(origImg(:,:,[1 1 1]));
white_pixel = sum(origImg(:));
%figure, imshow(origImg);
