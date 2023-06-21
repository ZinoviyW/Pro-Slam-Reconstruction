close force all;
clc;
clear variables;

%%Reconstruction Folder

% The anatomical 512x512 image
anatImg = imread("p.png");

load("segmentation.mat")
L = segmentation;
numLabels = 5;
% Load ELtons' separated images
recon_data = imread("crappy_image_i.png");

%%Go for the reconstruction
tic

orig_img = recon_data;
orig_ksp = fftshift(fft2(orig_img));
% The reconstruction function
recon_image = recon_algo(5,L,numLabels,0.5,orig_img,orig_ksp);
cd('/home/zinoviy/Desktop/ZInoviy/PRO-SLAM/PRO-SLAM/ForZinoviy/Stefan/Eltons_Data/Test')
imagesc(recon_image)
figure
imagesc(orig_img)

toc