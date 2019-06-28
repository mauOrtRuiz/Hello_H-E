ab=IIIne4; 
nrows = size(ab,1); 
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,1); 



nColors = 3; % repeat the clustering 3 times to avoid local minima 
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...                                       
    'Replicates',1);
pixel_labels = reshape(cluster_idx,nrows,ncols); 
%imshow(pixel_labels,[]), title('image labeled by cluster index');

segmented_images = cell(1,3); 
rgb_label = repmat(pixel_labels,[1 1 3]);  
for k = 1:nColors     
    color = he;     
    color(rgb_label ~= k) = 0;     
    segmented_images{k} = color; 
end
umb=min(cluster_center);
umbpi=mean([median(cluster_center) min(cluster_center)]);



