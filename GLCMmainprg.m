clear all;
clc;
%  folder in which fundus images exist
myFolder = 'C:\Mtech\Mtech thesis\matlab\Implementation\inpainting\InputAndOutput\ODinpaintingResults';       

if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
%  Creates a datastore for all images in your folder
ds = imageDatastore(myFolder);     

% image='Filename';
% contrast='contrast';
% correlation='correlation';
% energy='energy';
% entropy='entropy';
% homogenity='homogenity';
% variance='variance';  
% sumaverage='sumaverage';
% sumvariance='sumvariance';
% sumentropy='sumentropy';
% differencevariance='differencevariance';
% differenceentropy='differenceentropy';
% inf1='inf1';
% inf2='inf2';
count = 0;
% worksheetName = 'Results';
% cellReference = sprintf('A%d', count);
% Excelvalues = [image contrast correlation energy entropy homogenity variance sumaverage sumvariance sumentropy differencevariance differenceentropy inf1 inf2]
% xlswrite('C:\Mtech\Mtech thesis\matlab\Implementation\final\test\test.xls', Excelvalues, worksheetName, cellReference);


while hasdata(ds) 
    
    count = count+1;
    origImage = read(ds); 
    filename = ds.Files(count);
    GrayImage = rgb2gray(origImage);
    GLCM = graycomatrix(GrayImage,'Offset',[2 0]);
    out = GLCM_Features(GLCM);

%     [data,info] = read(ds);
     
    Excelvalues = [filename out.contrast out.correlation out.energy out.entropy out.homogenity out.variance  out.sumaverage out.sumvariance out.sumentropy out.differencevariance out.differenceentropy out.inf1 out.inf2];
    % Display image.
    drawnow; % Force display to update immediately.
    worksheetName = 'Sheet1';
    cellReference = sprintf('A%d', count+1);
    xlswrite('C:\Mtech\Mtech thesis\matlab\Implementation\final\test\test.xls', Excelvalues, worksheetName, cellReference);

    
    imshow(origImage);title('result ');
end
