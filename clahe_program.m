clear all;
clc;
%%   input datastore

    folder1 = 'C:\1_project_implementation\DRIVEdataset\images\*.tif'
    ds1 = imageDatastore(folder1); 
%     origImg = imread('C:\1_project_implementation\DRIVEdataset\images\11_test.tif');
    count = 1;
    
%%  CLAHE
    while hasdata(ds1) 
        
        origImg  = readimage(ds1,count);
        figure, imshow(origImg);title('Original image');
%         grayImg = rgb2gray(origImg);
        G = origImg(:,:,2);    
        figure, imshow(G);title('Green channel image');
%         R=origImg(:,:,1);
%         B=origImg(:,:,3);
        claheimg = adapthisteq(G);
        figure, imshow(claheimg);title('Green channel image');
%         imwrite(claheimg,'C:\Users\5559\Downloads\diagrams\clahe\11_clahe.jpg');
%         figure;
%                 subplot(1,2,1) ;
%                 imshow(G);
%                 title('Green channel image');
%                 subplot(1,2,2) ;
%                 imshow(claheimg);
%                 title('CLAHE output');
                
         
%%  save outputs

        folder2 = 'C:\1_project_implementation\programs\CLAHE\';
        if(count<=20)   
            newfilename1 = sprintf('%d_test.tif',count);
        else
            newfilename1 = sprintf('%d_training.tif',count);
        end
        fullfilename1 = fullfile(folder2,newfilename1);
        imwrite(claheimg,fullfilename1);
        
        count = count+1;
        
    end   
%%
