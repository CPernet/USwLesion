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
