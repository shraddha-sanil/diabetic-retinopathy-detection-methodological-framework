function [result] = otsuthresholding(img,count)
[rows,columns] = size(img);
imgSize = rows*columns;
t=0;
mu=0;
muT=0;
thetaT=0;
max=0;
hist = imhist(img);
px = zeros(1,256);
result = zeros(rows,columns);

%calculate pdf : probability density function
for i=1:256
    px(1,i) = hist(i,1)/imgSize;
end
%calculate mu
for i=1:256
    mu = mu + i * px(1,i);
end

%find biggest threshold
for i=1:256
    thetaT = thetaT + px(1,i);
    muT = muT + i*px(1,i);
    tempMax = ((muT - mu*thetaT)^2)/ (thetaT*(1-thetaT));
    if(tempMax>max)
        max = tempMax;
        t = i;
    end
end
%distinguish image by t 
for i=1:rows
    for j=1:columns
        if (img(i,j)<=t)
            result(i,j)=255;
        end
    end
end

%[x y d] = size(result);
result = uint8(result);
end

