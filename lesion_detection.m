clear all;
clc;
%%   input datastore

cleanfilename = 'C:\1_project_implementation\DRIVEdataset\images\*.tif'
ds1 = imageDatastore(cleanfilename); 

maskfilename1 = 'C:\1_project_implementation\BV\*.png';
ds2 = imageDatastore(maskfilename1);

% BV = imcomplement(BV);
% figure, imshow(BV);title('BV');
count = 1;

%%

while hasdata(ds1) 
    
    origImg  = readimage(ds1,count);
    gray = rgb2gray(origImg);
    figure, imshow(gray);title('gray');
    Ig = origImg(:,:,2);
    BV = readimage(ds2,count);
    % his = histeq(gray,5000);
    % figure, imshow(his);title('Bw');

    %% local contrast 
    edgeThreshold = 0.4;
    amount = 0.5;
    C = localcontrast(Ig, edgeThreshold, amount); %localcontrast : Edge-aware local contrast manipulation of images
    % figure, imshow(C);title('Ig local contrast');
    
    %%
    
    OT = otsuthresholding_(C,count);
   
    figure;
            subplot(1,2,1) ;
            imshow(origImg);
            title('Original Image');
            subplot(1,2,2) ;
            imshow(OT);
            title('lesion extraction');


count=count+1;

end
    
