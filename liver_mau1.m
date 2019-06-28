he=imread('Test1.jpg');
lab_he = rgb2lab(he);
ab = lab_he(:,:,2:3); 
nrows = size(ab,1); 
ncols = size(ab,2); 
ab = reshape(ab,nrows*ncols,2);  
nColors = 4;
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...                                       
    'Replicates',4);
pixel_labels = reshape(cluster_idx,nrows,ncols); 
%imshow(pixel_labels,[]), title('image labeled by cluster index');
segmented_images = cell(1,3); 
rgb_label = repmat(pixel_labels,[1 1 3]);  
for k = 1:nColors     
    color = he;     
    color(rgb_label ~= k) = 0;     
    segmented_images{k} = color; 
end
I1=segmented_images{1};
I2=segmented_images{2};
I3=segmented_images{3};
I4=segmented_images{4};
