% Read in original RGB image.
global rgbImage;
rgbImage = (im2double(imread('image1.bmp')));
[height, width, numberofChannels] = size(rgbImage);

% get the size of the original image
fid = fopen('image1.bmp');
fseek(fid, 0, 'eof');
originalSize = ftell(fid);
fclose(fid);

% Extract color channels.
redChannel = rgbImage(:,:,1); 
greenChannel = rgbImage(:,:,2); 
blueChannel = rgbImage(:,:,3); 

% Display each channel
figure('Name','red Channed','NumberTitle','off');
imshow(redChannel);
figure('Name','green Channel','NumberTitle','off');
imshow(greenChannel);
figure('Name','blue Channel','NumberTitle','off');
imshow(blueChannel);

% get the dct matrix of dimension 8*8
global T;
T = dctmtx(8);

max = 4; % the maximum m to be used
psnr = zeros(max,1);  % to store psnr of the different m used
compressionRatio = zeros(max, 1);   % to store the compression ratio

for i=1:max
    compressedSize = DCTCompress(i, height, width, redChannel, greenChannel, blueChannel);
    compressionRatio(i) = originalSize / compressedSize;
end

for i=1:max
    thePSNR = DCTDecompress(i, height, width);
    psnr(i) = thePSNR;
end

% display the psnr of the different m values
x = 1:1:max;
figure('Name','PSNR','NumberTitle','off');
bar(x, psnr);
xlabel('m');
ylabel('PSNR');

% display the compression ratio of the different m values
figure('Name','Compression Ratio with different m','NumberTitle','off');
bar(x, compressionRatio);
xlabel('m');
ylabel('Compression Ratio');


function compressedSize = DCTCompress(m, height, width, redChannel, greenChannel, blueChannel)
    
    % make three matrices to hold the r, g and b dct coefficients
    newWidth = width / (8 / m);
    newHeight = height / (8 / m);
    newRedChannel = zeros(newHeight, newWidth);
    newBlueChannel = zeros(newHeight, newWidth);
    newGreenChannel = zeros(newHeight, newWidth);

    % dct compress and mask coefficients around the m*m block of the 8*8
    % block (keeping only m*m coefficient) of each 8*8 coefficient
    for i=1:height/8 
        for j=1:width/8
           startX = (i-1) * m + 1;
           endX =  i * m;
           startY = (j-1) * m + 1;
           endY = j * m;
           
           block = compressBlock(blueChannel, i, j, m);
           newBlueChannel(startX:endX, startY:endY) = block;
           
           block = compressBlock(redChannel, i, j, m);
           newRedChannel(startX:endX, startY:endY) = block;
           
           block = compressBlock(greenChannel, i, j, m);
           newGreenChannel(startX:endX, startY:endY) = block;
        end
    end
    
    % save the coefficients to a file
    compressedRgbImage = (cat(3,newRedChannel,newGreenChannel,newBlueChannel));
    compressedRgbImage = half(compressedRgbImage);
    s1 = mfilename('fullpath');
    s2 = 'encodedm';
    s3 = '.mat';
    save([s1 s2 int2str(m) s3], 'compressedRgbImage');
   
    % compressed size is the size needed to store the coefficients and
    % because they are float minimum size in matlab is 2 bytes per
    % coefficient as minimum size float is (half) in matlab
    fid = fopen([s1 s2 int2str(m) s3]);
    fseek(fid, 0, 'eof');
    compressedSize = ftell(fid);
    fclose(fid);
    
end

function thePSNR = DCTDecompress(m, height, width)
    % make matrices to hold the channels of the decompressed image
    newRedChannel = zeros(height, width);
    newBlueChannel = zeros(height, width);
    newGreenChannel = zeros(height, width);
    
    % load the coefficients from the stored comrpessed file
    s1 = mfilename('fullpath');
    s2 = 'encodedm';
    s3 = '.mat';
    rgbImageCompressed = load([s1 s2 int2str(m) s3]);
    rgbImageCompressed = rgbImageCompressed.compressedRgbImage;
    
    % Extract color channels.
    encodedRedChannel = rgbImageCompressed(:,:,1); 
    encodedGreenChannel = rgbImageCompressed(:,:,2); 
    encodedBlueChannel = rgbImageCompressed(:,:,3);
    
    % inverse dct for each block 8*8 which consists of m*m coefficient and
    % zeros around that m*m block 
    for i=1:height/8 
        for j=1:width/8
           block = DecompressBlock(encodedBlueChannel, i, j, m);
           newBlueChannel((i-1) * 8 + 1: i * 8, (j-1) * 8 + 1: j * 8) = block;

           block = DecompressBlock(encodedRedChannel, i, j, m);
           newRedChannel((i-1) * 8 + 1: i * 8, (j-1) * 8 + 1: j * 8) = block;

           block = DecompressBlock(encodedGreenChannel, i, j, m);
           newGreenChannel((i-1) * 8 + 1: i * 8, (j-1) * 8 + 1: j * 8) = block;
        end
    end
    
    % save the decompressed image to file
    DecompressedRgbImage = (cat(3,newRedChannel,newGreenChannel,newBlueChannel));
    s1 = mfilename('fullpath');
    s2 = '.bmp';
    imwrite(DecompressedRgbImage, [s1 int2str(m) s2]);
    
    % calculate the psnr and compare the output with the original visually
    global rgbImage;
    thePSNR = psnr(DecompressedRgbImage, rgbImage);
    s1 = 'comparison between original and compressed image with m = ';
    figure('Name',[s1 int2str(m)],'NumberTitle','off');
    imshowpair(DecompressedRgbImage, rgbImage, 'montage');
end


function block = compressBlock(Channel, start_x, start_y, m)
  global T;
  temp = Channel((start_x-1) * 8 + 1: start_x * 8, (start_y-1) * 8 + 1: start_y * 8);
  temp = T * temp * T';
  block = temp(1:m, 1:m);
  
end

function block = DecompressBlock(Channel, start_x, start_y, m)
  global T;
  block = zeros(8, 8);
  block(1:m, 1:m) = Channel((start_x-1) * m + 1: start_x * m, (start_y-1) * m + 1: start_y * m);
  block = T' * block * T;
end