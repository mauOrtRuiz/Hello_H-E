%for loops to process automated the full SET
nRange = 1:3;
for n = nRange
    filename = ['ch99788_',num2str(n),'.tif'];
    he = imread(filename);
    imshow(he)
    n
    pause
end