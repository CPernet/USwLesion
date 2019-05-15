function SSIM = SSI(varargin)

% computes the structural similarity index
% Wang et al 2004 IEEE Trans Image Process 13, 600-612
%
% FORMAT SSIM = SSI(image1,image2,roi)
%        SSIM = SSI(image1,image2,roi,c1,c2)
%
% INPUT image1 and image2 are the image to process (char or matrix)
%       roi is the region use to mask imgages 1 and 2 (ie restrict SSI to that ROI,
%       could be ICV mask for the whole brain or a smaller region)
%       C1 and C2 (optional) are the mean and variance constant stabilizer
%       (useful to use the same C1 and C2 when compararing multiple SSIM -
%       otherwise this is proprotional to the images intensity range)
%
% OUTPUT SSIM is the Structural SIMilarity value
%
% Cyril Pernet, University of Edinburgh, May 2019

%% inputs
if nargin < 2
    errror('not enough arguments in')
else
    for in = 1:3
        if ischar(varargin{1})
            image1 = spm_read_vols(spm_vol(varargin{1}));
        else
            image1 = varargin{1};
        end
        if ischar(varargin{2})
            image2 = spm_read_vols(spm_vol(varargin{2}));
        else
            image2 = varargin{2};
        end
        if ischar(varargin{3})
            roi = spm_read_vols(spm_vol(varargin{3}));
        else
            roi= varargin{3};
        end
    end
end

if ~any(size(image1)==size(image2))
    error('imput images are of different dimensions')
end

if ~any(size(image1)==size(roi))
    error('the ROI has a different dimension than images 1 and 2')
else
    roi = roi>0; % make sure it's binary
end

%% compute

% quality index = ((2*mean(image1)*mean(image2))*(2*S) / (mean(image1)^2+mean(image2)^2)*(var(image1)+var(image2)
% if means or variances are close to 0, the quality index is unstable, we thus need to define usefull constant

if nargin == 4
    C1 = varargin{4};
    if isnumeric(C1)
        C2 = 3*C1;
    else
        error('the constant input is not a numerical value, see help')
    end
elseif nargin == 5
    C1 = varargin{4};
    C2 = varargin{5};
    if any([isnumeric(C1) isnumeric(C2)])
        error('the constant input are not a numerical values, see help')
    end
else
    C1 = (0.01*max([range(image1(:)) range(image2(:))])).^2; % 0.01 as in the matlab ssim function
    C2 = (0.03*max([range(image1(:)) range(image2(:))])).^2;
end

S = cov(image1(roi),image2(roi));
SSIM = ((2*nanmean(image1(roi))*nanmean(image2(roi))+C1)*(2*S(1,2)+C2)) / ...
    (nanmean(image1(roi))^2+nanmean(image2(roi))^2+C1)*(nanvar(image1(roi))+nanvar(image2(roi))+C2);
