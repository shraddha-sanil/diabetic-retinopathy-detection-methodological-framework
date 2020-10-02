%% GLCM_Features Computes a subset of GLCM features in a vectorized fashion.
function [out] = GLCM_Features(glcm)
% Input:
%   glcm - Ng x Ng x Ndir matrix (Ng - number of gray levels, Ndir - number
%       of directions for which GLCMs were computed)
%
% Output:
%   out - a structure containing values of Haralick's features in the
%       respective fields
%
% -- Features computed --
% Contrast: matlab/[1,2]                    (out.contr) [f4]
% Correlation: [1,2]                        (out.corrp) [f6]
% Energy (ASM): matlab / [1,2]              (out.energ) [f1]
% Entropy: [2]                              (out.entro) [f2]
% Homogeneity (Inv. Diff. Moment): [2]      (out.homop) [f7]
% Sum of sqaures: Variance [1]              (out.sosvh) [f12]
% Sum average [1]                           (out.savgh) [f13]
% Sum variance [1]                          (out.svarh) [f14]
% Sum entropy [1]                           (out.senth) [f15]
% Difference variance [1]                   (out.dvarh) [f16]
% Difference entropy [1]                    (out.denth) [f17]
% Information measure of correlation1 [1]   (out.inf1h) [f18]
% Informaiton measure of correlation2 [1]   (out.inf2h) [f19]

% if ((nargin > 2) || (nargin == 0))
%    error('Too many or too few input arguments. Enter GLCM and pairs.');
% else
%    if ((size(glcm,1) <= 1) || (size(glcm,2) <= 1))
%        error('The GLCM should be a 2-D or 3-D matrix.');
%     elseif ( size(glcm,1) ~= size(glcm,2) )
%        error('Each GLCM should be square with NumLevels rows and NumLevels cols');
%     end    
% end

%%
format long e
Ng = size(glcm,1);  % Number of gray levels=Ng
size_glcm_1 = size(glcm,1);
size_glcm_2 = size(glcm,2);
size_glcm_3 = size(glcm,3);
glcm_mean = zeros(size_glcm_3,1);
    % glcm_var  = zeros(size_glcm_3,1);
    % checked p_x p_y p_xplusy p_xminusy
p_x = zeros(size_glcm_1,size_glcm_3);   % Ng x #glcms[1]  
p_y = zeros(size_glcm_2,size_glcm_3);   % Ng x #glcms[1]
p_xminusy = zeros(size_glcm_1,1,size_glcm_3); %[1]
[Mj, Mi] = meshgrid(1:size_glcm_1,1:size_glcm_2);
Mi3d = repmat(Mi, [1 1 size_glcm_3]);
Mj3d = repmat(Mj, [1 1 size_glcm_3]);
T3d = repmat(abs(Mi - Mj) .^ 2, [1 1 size_glcm_3]);
contr = sum(sum(T3d .* glcm, 2), 1);
out.contrast = permute(contr, [2 3 1]);
glcm_sum = permute(sum(sum(glcm, 2) ,1), [2 3 1]);
for k = 1:size_glcm_3 % number glcms
    glcm(:,:,k) = glcm(:,:,k)./glcm_sum(k);    % Normalize each glcm
    glcm_mean(k) = mean2(glcm(:,:,k));     % compute mean after norm
end

%%
% Compute glcm_var
% n = size_glcm_1*size_glcm_2;
% mean3d = repmat(reshape(glcm_mean, [1 1 size_glcm_3]), [size_glcm_1 size_glcm_2 1]);
% d23d = (glcm - mean3d).^2;
% s = 1/(n-1) * sum(sum(d23d,2),1);
% glcm_var = permute(s, [2 3 1]);
    

%% Compute energ
t = glcm .^ 2;
energ3d = sum(sum(t,2),1);
out.energy = permute(energ3d, [2 3 1]);

%% Compute sosvh
uX2 = sum(sum(Mi3d .* glcm));
uX2_3d = repmat(uX2, [size_glcm_1 size_glcm_2 1]);
t = glcm .* (Mi3d - uX2_3d) .^ 2;
sosvh3d = sum(sum(t,2),1);
out.variance = permute(sosvh3d, [2 3 1]);

%% Compute entro
t = -glcm .* log(glcm + eps);
entro3d = sum(sum(t,2),1);
out.entropy = permute(entro3d, [2 3 1]);

%% Compute homop
t = glcm ./ repmat(1 + (Mi - Mj) .^ 2, [1 1 size_glcm_3]);
homop3d = sum(sum(t,2),1);
out.homogenity = permute(homop3d, [2 3 1]);

%% Compute p_x and p_y
for k = 1:size_glcm_3
    glcm_k = glcm(:,:,k);
    p_x(:,k) = sum(glcm_k, 2);
    p_y(:,k) = sum(glcm_k, 1);
end

% seq1 = 2:2*size_glcm_1;
% seq2 = 0:(size_glcm_1-1);
% Compute p_xplusy and p_xminusy
% for i = 1:size_glcm_1
%     for j = 1:size_glcm_2
%         NOTE: No need to check this condition - i+j ALWAYS falls in the
%         range of 2:2*size_glcm_1, as size_glcm_1 == size_glcm_2
%         if any(seq1 == i + j)
%             p_xplusy((i+j)-1,1,:) = p_xplusy((i+j)-1,1,:) + glcm(i,j,:);
%         end
        
%         NOTE: No need to check this condition - |i-j| ALWAYS falls in the
%         range of 0:(size_glcm_1-1), as size_glcm_1 == size_glcm_2
%         if any(seq2 == abs(i-j))
%             p_xminusy((abs(i-j))+1,1,:) = p_xminusy((abs(i-j))+1,1,:) + glcm(i,j,:);
%         end
%     end
% end
% p_xplusy = permute(p_xplusy, [1 3 2]);
% p_xminusy = permute(p_xminusy, [1 3 2]);
% -----

dim1 = size_glcm_1;
dim1m1 = dim1 - 1;
indexes = cell(2 * dim1 - 1,1);

% for k = 1:(dim1*2 - 1)
%     indexes{k} = k + (0:(k-1)) * dim1m1;
% end

indexes{1} = 1;
for k = 2:dim1
    indexes{k} = [indexes{k-1} + 1, k + (k-1) * dim1m1];
end
for k = (dim1+1):size(indexes,1)
    indexes{k} = indexes{k-1}(2:end) + 1;
end
for k = 1:size_glcm_3
    glcmk = glcm(:,:,k);
    for m = 1:length(indexes)
        p_xplusy(m,k) = sum(glcmk(indexes{m}));
    end
end

% dim1 = size_glcm_1;
% dim1m1 = dim1 - 1;
% 
% indexes = zeros(2 * dim1 - 1,dim1);
% indexes(1) = 1;
% for k = 2:dim1
%     indexes(k,1:k) = [indexes(k-1, 1:(k-1)) + 1, k + (k-1) * dim1m1];
% end
% for k = (dim1+1):size(indexes,1)
%     indexes(k,1:(size(indexes,1)-k+1)) = indexes(k-1,2:size(indexes,1)-k+2) + 1;
% end
% 
% indexes(indexes == 0) = 1;
% 
% modifiers = [(dim1 - 1):-1:0 1:(dim1-1)]';
% 
% p_xplusy = zeros(2 * dim1 - 1, size_glcm_3);
% for k = 1:size_glcm_3
%     glcmk = glcm(:,:,k);
%     p_xplusy(:,k) = sum(glcmk(indexes),2) - modifiers * glcm(1,1,k);
% end

p_xplusy = permute(p_xplusy, [1 3 2]);

dim1p1 = dim1 + 1;
indexesL = cell(dim1,1);
indexesU = cell(dim1,1);
indexesL{1} = 1:dim1p1:dim1^2;
indexesU{1} = 1:dim1p1:dim1^2;
for k = 2:dim1
    indexesL{k} = indexesL{k-1}(1:end-1) + 1;
    indexesU{k} = indexesU{k-1}(2:end) - 1;
end
indexesU{1} = [];
indexes = cellfun(@(c1, c2) [c1 c2], indexesL, indexesU, 'UniformOutput',false);
for k = (dim1+1):size(indexes,1)
    indexes{k} = indexes{k-1}(2:end) + 1;
end
for k = 1:size_glcm_3
    glcmk = glcm(:,:,k);
    for m = 1:length(indexes)
        p_xminusy(m,k) = sum(glcmk(indexes{m}));
    end
end
% -----
p_xplusy2d = permute(p_xplusy, [1 3 2]);

%% Compute sum average
out.sumaverage = sum(repmat((2:(2*size_glcm_1))', [1 size_glcm_3]) .* p_xplusy2d);

%% Compute sum entropy
out.sumentropy = -sum(p_xplusy2d .* log(p_xplusy2d + eps));

%% Compute sum variance with the help of sum entropy
t = repmat((0:(2*size_glcm_1-2))', [1 size_glcm_3]) - repmat(out.sumentropy, [size(p_xplusy,1) 1]);
out.sumvariance = sum(t.^2 .* p_xplusy2d);

%% Compute difference entropy
out.differenceentropy = -sum(p_xminusy .* log(p_xminusy + eps));
out.differenceentropy = squeeze(out.differenceentropy)';

%% Compute difference variance
p_xminusy2d = permute(p_xminusy, [1 3 2]);
% t = repmat(((0:(size_glcm_1-1)) .^ 2)', [1 size_glcm_3]);
% out.dvarh = sum(t .* p_xminusy2d);
%
% (use formula implemented in WND-CHARM - below)
ssq = sum(repmat((1:Ng)'.^2, [1 size_glcm_3]) .* p_xminusy2d);
s = sum((repmat((1:Ng)', [1 size_glcm_3]) .* p_xminusy2d) .^2);
out.differencevariance = ssq - s;
    
hx = -sum(p_x .* log(p_x + eps));
hy = -sum(p_y .* log(p_y + eps));
hxy = out.entropy;
p_x3d = reshape(p_x, [size_glcm_1 1 size_glcm_3]);
p_xM = repmat(p_x3d, [1 size_glcm_2 1]);
p_y3d = reshape(p_y, [1 size_glcm_2 size_glcm_3]);
p_yM = repmat(p_y3d, [size_glcm_1 1 1]);
p_xyM = p_xM .* p_yM;
p_xylogM = log(p_xyM + eps);
t = -(glcm .* p_xylogM);
hxy13d = sum(sum(t,2),1);
hxy1 = squeeze(hxy13d)';
t = -(p_xyM .* p_xylogM);
hxy23d = sum(sum(t,2),1);
hxy2 = squeeze(hxy23d)';
out.inf1 = ( hxy - hxy1 ) ./ ( max([hx;hy]) );
out.inf2 = ( 1 - exp( -2*( hxy2 - hxy ) ) ) .^ 0.5;
% Compute u_x
t = Mi3d .* glcm;
u_x3d = sum(sum(t,2),1);
u_x = squeeze(u_x3d)';
% Compute u_y
t = Mj3d .* glcm;
u_y3d = sum(sum(t,2),1);
u_y = squeeze(u_y3d)';
% Compute s_x
u_xM = repmat(u_x3d, [size_glcm_1 size_glcm_2]);
t = (repmat(Mi, [1 1 size_glcm_3]) - u_xM) .^ 2 .* glcm;
s_x3d = sum(sum(t,2),1);
s_x = squeeze(s_x3d) .^ 0.5;
% Compute s_y
u_yM = repmat(u_y3d, [size_glcm_1 size_glcm_2]);
t = (repmat(Mj, [1 1 size_glcm_3]) - u_yM) .^ 2 .* glcm;
s_y3d = sum(sum(t,2),1);
s_y = squeeze(s_y3d) .^ 0.5;
% Compute corp
t = repmat(Mi .* Mj, [1 1 size_glcm_3]) .* glcm;
corp3d = sum(sum(t,2),1);
corp = squeeze(corp3d);
% Compute corm
% t = (Mi3d - u_xM) .* (Mj3d - u_yM) .* glcm;
% corm3d = sum(sum(t,2),1);
% corm = permute(corm3d, [2 3 1])';
out.correlation = (corp' - u_x .* u_y) ./ (s_x .* s_y)';

