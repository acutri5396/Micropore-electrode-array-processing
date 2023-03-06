 close all
clear variables
clc

%Varun Mannam and Vignesh Sundaresan
%date: 6th Febraruy 2021
%extract the pore mean value (by subtracting the background)
format long
font = 14;
linewidth = 2;
%results:
tic
%load the tiff file first frame

addpath('E:\Wide Field\2023\February\02132023\100pc succinate 24hr\MPS6');

file_name = '100pc succinate laser 20mW MPS6.tiff';

a1 = imread(file_name,1);
a11 = size(a1);

numimgs = size(imfinfo(file_name),1);
ipx = zeros(a11(1),a11(2),numimgs);

for i=1:numimgs
    ip = imread(file_name,i); %read ith time frame
    ip = double(ip); %covert to double
    ipx(:,:,i) = ip; %store to matrix
end

for f=1:1:numimgs
        val(:,:,1)=imread(file_name,i);   %read in each pixel at frame f
        localmax=max(val,[],3); %compare the value for each pixel at frame f to the max value so far
        val(:,:,2)=localmax;    %keep the highest pixel intensity and store
       end

% new_image = zeros(a11(1),a11(2));
% for i=1:a11(1)
%     for j=1:a11(2)
%         temp1 = ip(i,j,:);
%         new_image(i,j) = max(temp1);
%     end
% end
        
%get the first image
% a1 = ipx(:,:,1);
a1 = localmax;

%bianry mask
a1_mask = segmentImage(a1); %segmenting image
a1_mask = double(a1_mask);
a1_mask1 = uint8(a1_mask);

%get all the labels in this image
labeledImage1 = bwlabel(a1_mask1, 8);
%coloredLabels1 = label2rgb (labeledImage1, 'hsv', 'k', 'shuffle'); % pseudo random color labels

%get the measurements here
blobMeasurements1 = regionprops(labeledImage1, a1, 'all');
numberOfBlobs1 = size(blobMeasurements1, 1); %number of pours with atleast 1 pixel
boundaries1 = bwboundaries(a1_mask); %what are boundaries of these regions

a2 = [blobMeasurements1.Area].'; %get the number of pixels in each pour
blobMeasurements2 = blobMeasurements1(a2>1); %atleast 2x2 image %select pours atleast 4 pixels, 
%less than that it is not a pore, but it is the single pixel -> misleads
%our analysis

a2_max = [blobMeasurements2.MaxIntensity].'; %get the max intensity of each pore
a2_min = [blobMeasurements2.MinIntensity].'; %get the min intensity of each pore
blobMeasurements3 = blobMeasurements2(a2_max-a2_min>150); %atleast 100 intesnity counts difference, 
%less than that it is not a pore, but it is mislocalization for our
%analysis


%get the centroids of the pore
a3 = [blobMeasurements3.Centroid].';
n1 = size(a3,1)/2; %number of pours with valid pixels
n1 =  round(n1);
x_c = zeros(n1,1);
y_c = zeros(n1,1);
for i = 1:n1 %ge the centroid values -> of each pore
    x_c(i) = round(a3((i-1)*2+1)); 
    y_c(i) = round(a3((i)*2)); 
end


%% plot image with centroids
figure,
imagesc(a1,[min(a1(:)) max(a1(:))]);
colormap('gray'), colorbar;
%axis image; % Make sure image is not artificially stretched because of screen's aspect ratio.
hold on;
for i = 1:n1
	plot(x_c(i), y_c(i), 'm-*', 'LineWidth', 2);
end
hold off;
title('The centroids of the pores');
set(gca,'FontSize',font);

%% get the mean of the pore (after subtracting the background intensities)
n_pixels_side = 4;
proc_im = zeros(n1,numimgs);
for i=1:n1 %ith pours
    %background image cordinates
    b_xc = x_c(i)-(n_pixels_side*2+1);
    if b_xc < 1 %%safety boundary
        b_xc = 1;
    end
    b_yc = y_c(i)-(n_pixels_side*2+1);
    if b_yc < 1 %safety boundary
        b_yc = 1;
    end
    
    for j=1:numimgs %each image here %over time 
        a_temp = ipx(:,:,j); %jth image
         b_im = a_temp(b_yc-n_pixels_side:b_yc+n_pixels_side,b_xc-n_pixels_side:b_xc+n_pixels_side);
         b_im_mean = mean(b_im(:));
         
         b_im_mean_1(i,j)= mean(b_im(:)); %store the background data into a new variable for data analysis later 
%b_im_mean = 0;

        im_prc = a_temp(y_c(i)-n_pixels_side:y_c(i)+n_pixels_side,x_c(i)-n_pixels_side:x_c(i)+n_pixels_side);
        im_prc_mean = mean(im_prc(:));
        
        %save the value here
        proc_im(i,j) = im_prc_mean - b_im_mean; %pour 3x3 mean - background mean
    end
    
end


%% desired index of the pour
fprintf('Max number of pours (start =1) is %d in this image.\n',n1); %index2 shows the pulse
index = input('Enter the index of the pore to visualization? ', 's'); 
index = str2double(index);

x1= [1:numimgs]'; %multiply the frame rate to convert to seconds
figure,
plot(x1,proc_im(index,:),'b-*','Linewidth',2);
title(string('Mean pixel value over time of Pore:') + num2str(index));
xlabel('#Frame number');
ylabel('Mean pixel value');
set(gca,'FontSize',font);


%%
%numbered image
 h =  imagesc(localmax);
    colormap 'gray'
    hold on
  for ll = 1:n1
    rectangle('position',[(x_c(ll)-(n_pixels_side)) (y_c(ll)-(n_pixels_side)) 9 9],'EdgeColor','yellow')
    axis([0 a11(1) 0 a11(2)])
    txt1 = num2str(ll);
    text((x_c(ll)-((n_pixels_side)+8)),(y_c(ll)-((n_pixels_side)+8)),txt1,'Color','yellow','FontSize',14)    
  end
   
    saveas(h,sprintf('roi_numbered.png')); 

%%
%export data as a single txt file
x1= [1:numimgs]'; %multiply the frame rate to convert to seconds
data(:,1) = x1*0.1; %convert  frame # to time 0.1 is the integration time
data(:,2:n1+1) = proc_im.';
save([file_name  '_data'  '.txt'],'data','-ascii');

%%
%export bakcground as a single txt file
x1= [1:numimgs]'; %multiply the frame rate to convert to seconds
data1(:,1) = x1*0.1; %convert  frame # to time 0.1 is the integration time
data1(:,2:n1+1) = b_im_mean_1.';
save([file_name  '_bkgd'  '.txt'],'data','-ascii');

%%
 %Plotting and saving individual timetrace
 for kk = 1:n1
 g = figure ();
        
  plot(data(:,1), data(:,kk+1),'LineWidth',1.5,'Color','k');
  xlabel('Time (s)','FontSize',12)
  xlim([0 max(data(:,1))])       
  ylabel('Intensity (a. u.)', 'FontSize',12)
  set(gca,'FontSize',12)
  %axis tight
  saveas(g,sprintf('pore-%d.png',kk));    
  close (g);
  end
%%
toc

%end of code