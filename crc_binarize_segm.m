function [fn_out,fn_nc] = crc_binarize_segm(fn_in,fn_msk,opt)
%
% Simple routine to binarize 2D or 3D probability images.
% For example, the posterior probability images (ci or wci) generated by 
% segmenting structural images and return binary (bci or bwci)images.
%
% The idea is to assign each voxel to one specific tissue class according
% to one or more criteria. A voxel is said to be of class i (bci = 1)
% if the probability of being i is higher than
% - a fixed threshold, i.e. a classic binarization;
% - that of being any other class, i.e. a "most likely" (ML) binarization;
% - the sum of probabilities to be from the other classes, i.e. a "more
%   likely than the rest" (mltr) binarization.
%  An extra image (nc) is generated containing all the voxels that are not
%  assigned to any class (no-class voxels).
%
% FORMAT
%   [fn_out,fn_nc] = crc_binarize_segm(fn_in,fn_msk,opt)
% or
%   [vx_out,vx_nc] = crc_binarize_segm(vx_in,vx_msk,opt)
%
% INPUT:
% - fn_in/vx_in  : file name of the posterior probability images (e.g. ci) 
%                  or 4D array with the pre-loaded probability images.
% - fn_msk/vx_msk: binary mask defining which bit of the image is to be
%                  considered or 3D array with the preloaded mask.
%                  By default [], i.e. no masking.
% - opt    : set of options
%   .bin_classic : classic binarization [1, def.]
%   .bin_ML      : "most likely" (ML) binarization [1, def.]
%   .bin_mltr    : "more likely than the rest" (mltr) binarization [0, def.]
%   .thr         : fixed threshold for classic binarization [.2, def]
%   .singleImg   : force consider there is only one single image for input
%                  (this is to avoid some dimension issue when handling the
%                  images) [0, def].
%   .ListOut     : index of input images that need to be written out, only
%                  when images are written on disk. If 0, then all images 
%                  are written [0, def.].
%
% OUTPUT:
% - fn_out/vx_out : file name of the binarized probability images (e.g. 
%                   binarized tissue classes, bci) or a 4D array with the 
%                   voxel values
% - fn_nc/vx_nc   : file name of the no-class voxel image (nc) or 3D array 
%                   with the voxel values (within the mask if provided).
%
% NOTE:
% - the 'ML' and 'mltr' binarization are mutually exclusive, i.e. only one 
%   of these 2 can be applied. ML will have priority if both are selected.
% - when classic and ML/mltr binarization are used together, then voxels
%   must meet both criteria (like a 'AND' operation).
% - all input images and mask should have excactly the same number of
%   voxels and orientation!
% - for a single 3D image, you should really set the option 'singleImg' to 
%   'true', otherwise the routine will crash. This will also desable any 
%   ML/mltr binarization, as they require at least 2 images!
%
%_______________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% Deal with options
opt_def = struct(... % Default options
    'bin_classic', true, ...
    'bin_ML', true, ...
    'bin_mltr', false, ...
    'singleImg', false, ...
    'ListOut', 0, ...
    'thr', .2);
opt = crc_check_flag(opt_def,opt);

if opt.singleImg
    opt.bin_ML = false;   %\_ No comparison across multiple maps possible
    opt.bin_mltr = false; %/
end

if opt.bin_ML && opt.bin_mltr
    opt.bin_ML = true;    %\_ Only one is possible -> priority to ML
    opt.bin_mltr = false; %/
    warning('Cannot do both ML and mltr binarization -> keeping ML!');
end

%% Deal with input files/data
if nargin<2, fn_msk = []; end % if no mask
if nargin<1
    warning('Wrong input format, see help here under.')
    help crc_binarize_segm;
    return
end

% input images
if ischar(fn_in)
    V_in = spm_vol(fn_in);
    spm_check_orientations(V_in);
    vx_in = spm_read_vols(V_in);
    save_img = true;
    if any(opt.ListOut)
        save_list = opt.ListOut;
    else
        save_list = 1:numel(V_in);
    end
elseif isnumeric(fn_in)
    vx_in = fn_in;
    save_img = false;
else
    error('Wrong input format!');
end

% get sizes and vectorize images
SZ_in = size(vx_in); % could be a 4D (N 3D images) or 2D (N vectorized images) array
if opt.singleImg
    Nb_in = 1; % case of a single image...
    SZ_in = [SZ_in 1];
    vx_in = vx_in(:);
else
    Nb_in = SZ_in(end); % Last size should be number of images
    vx_in = reshape(vx_in,[prod(SZ_in(1:end-1)) Nb_in]);
end

% mask image
if isempty(fn_msk)
    apply_mask = false; % No masking
else
    apply_mask = true; % use mask but should do further checks!
    if ischar(fn_msk) && ischar(fn_in)
        V_msk = spm_vol(fn_msk);
        spm_check_orientations([V_in(1) ; V_msk])
        vx_msk = spm_read_vols(V_msk);
        vx_msk = vx_msk(:);
    elseif isnumerical(fn_msk) && isnumeric(fn_in)
        vx_msk = fn_msk;
        if any(size(vx_msk~=SZ_in(1:end-1)))
            warning(['Mismatch between mask and input volume.' ...
                'No mask applied!'])
            apply_mask = false;
        else
            vx_msk = vx_msk(:);
        end
    else
        warning(['Wrong mask volume format or incompatibility with ' ...
            'input volumes. No mask applied!'])
        apply_mask = false;
    end
end

%% Apply the binarization
% One at a time and mask at the end.

% initialize vectorized binary volume(s)
bi = zeros(prod(SZ_in(1:end-1)),Nb_in);

if opt.bin_classic
    bi = vx_in>opt.thr;
    bi_classic = bi;
end

if opt.bin_ML % should have a p higher than all the others individually
    list_c = 1:Nb_in;
    u_vec = ones(1,Nb_in-1);
    for ii=list_c
        list_c_ii = list_c; list_c_ii(ii) = []; % list of others
        bi(:,ii) = sum((vx_in(:,ii)*u_vec)>vx_in(:,list_c_ii),2)==(Nb_in-1);
    end
    if opt.bin_classic
        bi = bi.*bi_classic;
    end
end

if opt.bin_mltr % should have a p higher than the sum of the others
    list_c = 1:Nb_in;
    u_vec = zeros(1,Nb_in-1);
    for ii=1:list_c %#ok<*BDSCI>
        list_c_ii = list_c; list_c_ii(ii) = []; % list of others
        bi(:,ii) = (vx_in(:,ii)*u_vec)>sum(vx_in(:,list_c_ii),2);
    end
    if opt.bin_classic
        bi = bi.*bi_classic;
    end
end

% left out voxels, i.e. not-classified (nc)
nc = ~sum(bi,2);

% apply mask, if requested
if apply_mask
    for ii=1:Nb_in
        bi(:,ii) = bi(:,ii) .* vx_msk;
    end
    nc = nc .* vx_msk;
end

%% write down images or pass values
if save_img
    % create images from bi and nc
    Nb_out = numel(save_list);
    fn_out = spm_file(fn_in(save_list,:), 'prefix', 'b', 'number', '');
    V_out = V_in(1:Nb_out);
    for ii=1:Nb_out
        V_out(ii).fname = deblank(fn_out(ii,:));
        V_out(ii).dt(1) = 2;
        V_out(ii) = spm_write_vol(V_out(ii), ...
            reshape(bi(:,save_list(ii)),SZ_in(1:end-1)));
    end
    fn_nc = spm_file(fn_in(1,:), 'prefix', 'nc', 'number', '');
    V_nc = V_in(1);
    V_nc.fname = deblank(fn_nc);
    V_nc.dt(1) = 2;
    V_nc = spm_write_vol(V_nc,reshape(nc,SZ_in(1:end-1))); %#ok<*NASGU>
else
    % just change the name of variables for the output
    fn_out = reshape(bi,SZ_in);
    fn_nc = reshape(nc,SZ_in(1:end-1));
end

end


