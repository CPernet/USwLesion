function [mJ,mHd,overlap,Dice] = compare2gdtruth(varargin)

% FORMAT: [mJ,mHd,overlap,Dice] = compare2gdtruth(c1,c2,c3,c4,ground_truth,save_option);
%
% INPUT c1, c2, c3, c4 
%          ground_truth
%        save_option 'on' or 'off' 
%

% OUTPUT mJ,mHd,overlap,Dice
%        
 
%% default
if nargin < 5
    error('not enough inputs')
elseif nargin ==5
    save_option = 'off';
elseif nargin ==6
    save_option = varargin{6};
end

%% 1st load images and check all dimension match

V1 = spm_vol(varargin{1});
V2 = spm_vol(varargin{2});
V3 = spm_vol(varargin{3});
V4 = spm_vol(varargin{4});

if isequal(V1.dim,V2.dim,V3.dim,V4.dim);
    disp('dimensions do match');
else
    error('image dimensions do not match');
end

c1_vol = spm_read_vols(V1);
c2_vol = spm_read_vols(V2);
c3_vol = spm_read_vols(V3);
c4_vol = spm_read_vols(V4);

%% create binary images of c3

% tumour prob > sum of all other classes 
% normal tissue is c1+c2+c4

normal_tissue = c1_vol + c2_vol + c4_vol;
c3a = c3_vol > normal_tissue;
c3a = extent_thresold(c3a);
if strcmp(save_option,'on')
    V3a = V3;
    [root,filename,ext]=fileparts(V3.fname);
    V3a.fname = [root filesep filename '_bigger_than_sum_of_others' ext];
    V3a.descrip = 'c3 > c1 + c2 + c4;';
    spm_write_vol(V3a,c3a);
end

% tumour prob > any of other classes

c3b = (c3_vol > c1_vol) | (c3_vol > c2_vol) | (c3_vol > c4_vol);
c3b = extent_thresold(c3b);
if strcmp(save_option,'on')
    V3b = V3;
    [root,filename,ext]=fileparts(V3.fname);
    V3b.fname = [root filesep filename '_bigger_than_at_least_one' ext];
    V3b.descrip = 'c3 > c1 | c3 > c2 | c3 > c4;';
    spm_write_vol(V3b,c3b);
end

c3c = (c3_vol > c1_vol) & (c3_vol > c2_vol) & (c3_vol > c4_vol);
c3c = extent_thresold(c3c);
if strcmp(save_option,'on')
    V3c = V3;
    [root,filename,ext]=fileparts(V3.fname);
    V3c.fname = [root filesep filename '_bigger_than_each' ext];
    V3c.descrip = 'c3 > c1 & c3 > c2 & c3 > c4';
    spm_write_vol(V3c,c3c);
end

% tumour > 99%
c3d = c3_vol > 0.99;
c3d= extent_thresold(c3d);
if strcmp(save_option,'on')
    V3d = V3;
    [root,filename,ext]=fileparts(V3.fname);
    V3d.fname = [root filesep filename '_bigger_than_99' ext];
    V3d.descrip = 'c3 > 99';
    spm_write_vol(V3d,c3d);
end

% tumour > 50%
c3e = c3_vol > 0.5;
c3e = extent_thresold(c3e);
if strcmp(save_option,'on')
    V3e = V3;
    [root,filename,ext]=fileparts(V3.fname);
    V3e.fname = [root filesep filename '_bigger_than_50' ext];
    V3e.descrip = 'c3 > 50';
    spm_write_vol(V3e,c3e);
end
 

clear c1 c2 c3 c4

%% now compute overlap

GT = spm_read_vols(spm_vol(varargin{5}));

[mJ(1),mHd(1),overlap(1)] = image_overlap(c3a,GT);
Dice(1) = 2*overlap(1).voxel.tp / (2*overlap(1).voxel.tp+2*overlap(1).voxel.fp+2*overlap(1).voxel.fn);

[mJ(2),mHd(2),overlap(2)] = image_overlap(c3b,GT);
Dice(2) = 2*overlap(2).voxel.tp / (2*overlap(2).voxel.tp+2*overlap(2).voxel.fp+2*overlap(2).voxel.fn);

[mJ(3),mHd(3),overlap(3)] = image_overlap(c3c,GT);
Dice(3) = 2*overlap(3).voxel.tp / (2*overlap(3).voxel.tp+2*overlap(3).voxel.fp+2*overlap(3).voxel.fn);

[mJ(4),mHd(4),overlap(4)] = image_overlap(c3d,GT);
Dice(4) = 2*overlap(4).voxel.tp / (2*overlap(4).voxel.tp+2*overlap(4).voxel.fp+2*overlap(4).voxel.fn);

[mJ(5),mHd(5),overlap(5)] = image_overlap(c3e,GT);
Dice(5) = 2*overlap(5).voxel.tp / (2*overlap(5).voxel.tp+2*overlap(5).voxel.fp+2*overlap(5).voxel.fn);


end


function mapout= extent_thresold(mapin)
%% remove all clusters but the biggest

CC = bwconncomp(mapin,18);
cluster_size = NaN(CC.NumObjects,1);
for n=1:CC.NumObjects
    cluster_size(n) = length(CC.PixelIdxList{n});
end
[~,position]=max(cluster_size);
mapout = zeros(size(mapin));
mapout(CC.PixelIdxList{position}) = 1;

end





