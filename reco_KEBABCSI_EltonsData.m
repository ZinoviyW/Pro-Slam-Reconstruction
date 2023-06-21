%% init and dataset listing
close force all;
warning off all;
clc;
clear variables;
tic
%path to separated data folder
invivoCSIdataFolder=['/home/zinoviy/Desktop/ZInoviy/PRO-SLAM/' ...
    'PRO-SLAM/ForZinoviy/Stefan/Eltons_Data/Separated_Data'];
%path to anatomic images folder
anatomicImagesFolder=['/home/zinoviy/Desktop/ZInoviy/PRO-SLAM/' ...
    'PRO-SLAM/ForZinoviy/Stefan/Eltons_Data/Anatomical_Data'];
cd(invivoCSIdataFolder);
%% reconstruction parameters
ProSLAM.DoSLAM =true;
ProSLAM.num_iter=3; % 5 iterations worked best for phantoms and 3 for in-vivo
ProSLAM.InteractiveContours = false; % interactive contour drawing
ProSLAM.interactiveCS = true; % interactive spectral window
ProSLAM.forceZeroPhase = false;
ProSLAM.ShowMovie = false;
ProSLAM.calcSNRfactor = false;
ProSLAM.sigmaGaussian = 1; % sigma is usually 1. use sigma=8 for the glucose thermal phantom

%spatial zero-filling
% For SLAM in phantoms use ZF of ~46, for in-vivo use 30
if (ProSLAM.DoSLAM)
%     spatialZeroFillingPts=40; % 30 works for Stefan's data % 36 is nice on some of Anne's data
    spatialZeroFillingPts = 30;
end
%0=no change
%x=x voxels added (MUST BE A MULTIPLE OF 2)

%spatial apodization
spatialApodization=1;
%0=no spatial apodization
%1=gaussian k-space apodization

%spectral zero-filling
spectralZeroFillingPts=0;
%0=no change
%x=x zero time points added
%spectral apodization

%apodization
spectralDampingFactor=80; %was 20Hz

%quantification map size
if (ProSLAM.DoSLAM)
    nVoxelsQuantificationMaps=52;  % I've used 52 for the phantom, 32 for in-vivo
else
    nVoxelsQuantificationMaps=32; % was 20 for FT
end
%x=number of voxels of the maps

%how to normalize quantification maps?
normalizeQuantificationMaps=3;
%0=do nothing
%1=normalize - repetitions (at least one red voxel on each repetitions) !!! use ONLY for display !!!
%2=normalize - slices (at least one red voxel for each slice)
%3=normalize - all the data (at least one red voxel)

%how to normalize time curves?
normalizeTimeCurves=0;
%0=do nothing
%1=normalize each time curve to the maximum

%select the number of repetion to process
NRprocessOnly=0;
%0=process all the repetitions
%x=process only the x first repetitions

%average the repetitions into one image
NRaverage=0;
%0=do not average
%1=average all the repetitions or only NRprocessOnly (if >0)
%"a la radial" reconstruction will not be applied!

%interpolation "a la radial"
moreRepetitions=0;
%x=x more images artificially created

%% ELton's data
% GF - hyperpolarized 13C pyruvate + LDH in-vitro 10x10
iScan=8;
iScanAnatomy={9};
anatomicImageFilenames={'Anatomical.png'}; %0715c.png
zeroTime=0; %s
%cd ./LDH/
cd ./expt01/
cd(num2str(iScan));
ProSLAM.DisplayProtonSlice = 15; 
ProSLAM.NumCompartments = -1; 
%% display multislice anatomy
if(exist('iScanAnatomy','var'))
    scanFolder=pwd;
    %check 13C FOV
    lirePVParam('/home/tangir','truc');
    %% 
    PVM_Fov=lirePVParam('.','PVM_Fov');
    ACQ_trim=lirePVParam('.','ACQ_trim');
    PVM_GradCalConst=lirePVParam('.','PVM_GradCalConst');
    ACQ_O1_list=lirePVParam('.','ACQ_O1_list');
    PVM_SliceThick=lirePVParam('.','PVM_SliceThick');
    slicePos=ACQ_O1_list./(ACQ_trim(1,1)/100.0*PVM_GradCalConst);
    
    for i=1:length(iScanAnatomy)
        %read 1H anatomic images
        cd([scanFolder,'/../',num2str(iScanAnatomy{i})]);
        lirePVParam('/home/tangir','truc');
        PVM_Matrix_proton=lirePVParam('.','PVM_Matrix');
        PVM_Fov_proton=lirePVParam('.','PVM_Fov');
        PVM_SliceThick_proton=lirePVParam('.','PVM_SliceThick');
        PVM_SPackArrSliceOffset=lirePVParam('.','PVM_SPackArrSliceOffset');
        
        %rebuild 1H slice positionning
        nSlices_proton=lirePVParam('.','PVM_SPackArrNSlices');
        slicePos_proton=PVM_SliceThick_proton.*( (1:nSlices_proton) - (nSlices_proton/2) )+PVM_SPackArrSliceOffset-PVM_SliceThick_proton/2;
        sliceIndex_proton=1:nSlices_proton;
        
        %keep only the 1H slices included in the 13C slice
        sliceIndex_proton_sel=[];
        for s=1:length(slicePos)
            tmp=sliceIndex_proton(slicePos_proton<(slicePos(s)+PVM_SliceThick/2) & slicePos_proton>(slicePos(s)-PVM_SliceThick/2));
            sliceIndex_proton_sel=[sliceIndex_proton_sel,tmp];
        end
        
        %read 2dseq data
        cd('pdata/1');
        fid=fopen('2dseq','rb','ieee-le');
        data_proton=fread(fid,'int16');
        fclose(fid);
        
        %reshape data
        data_proton=reshape(data_proton,[PVM_Matrix_proton,nSlices_proton]);
        
        %rotate/flip/smooth
        psf=fspecial('gaussian',5,0.5);
        for s=1:nSlices_proton
            data_proton(:,:,s)=rot90(data_proton(:,:,s),-1);
            data_proton(:,:,s)=data_proton(:,end:-1:1,s);
            data_proton(:,:,s)=imfilter(data_proton(:,:,s),psf,'conv');
        end
        
        %check if everything is square
        if( size(data_proton,1)~=size(data_proton,2) || ...
                PVM_Fov(1)~=PVM_Fov(2) || ...
                PVM_Fov_proton(1)~=PVM_Fov_proton(2))
            error('Mmmh... Something (Matrix or FOV) is not square!');
        end
        
        %crop?
        if(~isequal(PVM_Fov,PVM_Fov_proton))
            newsize_proton=round(PVM_Fov(1)/PVM_Fov_proton(1)*size(data_proton,1));
            edge_proton=round((size(data_proton,1)-newsize_proton)/2);
            data_proton_cropped=data_proton(edge_proton:(edge_proton+newsize_proton),edge_proton:(edge_proton+newsize_proton),:);
        else
            data_proton_cropped=data_proton;
        end
        
        %resize
        % GF GF GF
        % UNCOMMENT WHEN POSSIBLE
        %data_proton_cropped=imresize(data_proton_cropped,3);
        
        %normalize RGB 0-255
        data_proton_cropped_uint=data_proton_cropped-min(data_proton_cropped(:));
        data_proton_cropped_uint=255.*data_proton_cropped_uint./max(data_proton_cropped_uint(:));
        data_proton_cropped_uint=uint8(data_proton_cropped_uint);
        data_proton_cropped_uint=repmat(data_proton_cropped_uint,[1 1 1 3]);
        data_proton_cropped_uint=permute(data_proton_cropped_uint,[1 2 4 3]);
        
        %1H slices not included in 13C slice --> red
        sliceIndex_proton_nonsel=setdiff(1:nSlices_proton,sliceIndex_proton_sel);
        data_proton_cropped_uint(:,:,2,sliceIndex_proton_nonsel)=0;
        data_proton_cropped_uint(:,:,3,sliceIndex_proton_nonsel)=0;
        
        %display stack
        ShowProtonStack = ProSLAM.ShowMovie;
        if ShowProtonStack
            implay(uint8(data_proton_cropped_uint));
        end
    end
    
    cd(scanFolder);
end

%% read data and start reconstruction

%read acq parameters
lirePVParam('/home/tangir','truc');
ACQ_trim=lirePVParam('.','ACQ_trim');
BF1=lirePVParam('.','BF1');
CSIMatrix=lirePVParam('.','CSIMatrix');
NR=lirePVParam('.','NR');
nSlices=lirePVParam('.','PVM_SPackArrNSlices');
PVM_DigShift=lirePVParam('.','PVM_DigShift');
PVM_ExcPulseAngle=lirePVParam('.','PVM_ExcPulseAngle');
PVM_GradCalConst=lirePVParam('.','PVM_GradCalConst'); %Hz/mm
PVM_RepetitionTime=lirePVParam('.','PVM_RepetitionTime');
PVM_SpecAcquisitionTime=lirePVParam('.','PVM_SpecAcquisitionTime');
PVM_SpecSWH=lirePVParam('.','PVM_SpecSWH');
PVM_SPackArrSliceOrient=lirePVParam('.','PVM_SPackArrSliceOrient');
snailYesNo=lirePVParam('.','snailYesNo');
rampT1=lirePVParam('.','rampT1');

%display some acq parameters for info
disp(['CSIMatrix=',num2str(CSIMatrix)]);
disp(['NR=',num2str(NR)]);
disp(['nSlices=',num2str(nSlices)]);
disp(['PVM_ExcPulseAngle=',num2str(PVM_ExcPulseAngle)]);
disp(['PVM_SpecAcquisitionTime=',num2str(PVM_SpecAcquisitionTime)]);
disp(['PVM_SliceThick=',num2str(PVM_SliceThick)]);
disp(['snailYesNo=',snailYesNo]);
disp('-----------------------------------------------------------------------------');

%warning the user
warning on all;
if(spatialZeroFillingPts>0)
    warning('spatialZeroFillingPts>0! We are zero-filling the spatial dimensions!');
end
if(spectralZeroFillingPts>0)
    warning('spectralZeroFillingPts>0! We are zero-filling the spectral dimension!');
end
if(NRprocessOnly>0)
    warning(['NRprocessOnly>0! We are going to process only the ',num2str(NRprocessOnly),' first repetitions!']);
end
if(moreRepetitions>1)
    warning('moreRepetitions>1! We are doing a la radiale processing!');
    if(snailYesNo(1)~='y')
        error('but this is not a good idea: the scan is not a centric k-space fill!');
    end
end

%read data
fid=fopen('fid','rb','ieee-le');
if(fid==-1)
    warning('No fid file! fid.orig?');
    fid=fopen('fid.orig','rb','ieee-le');
end
data=fread(fid,'int32');
fclose(fid);

disp('-----------------------------------------------------------------------------');
warning off all;

%reshape data
nVoxels=CSIMatrix(1);
nPts=length(data)/(2*(nVoxels^2)*NR*nSlices);
data=reshape(data,[2,nPts,nSlices,nVoxels,nVoxels,NR]);
data=data(1,:,:,:,:,:)+1i.*data(2,:,:,:,:,:); % GF - this is the CSI data [1   fid     1    x    y    time-repetition]

%reduce the number of repetitions to process
if(NRprocessOnly>0 && NRprocessOnly<NR)
    disp(['Only ',num2str(NRprocessOnly),' repetitions out of ',num2str(NR),' will be processed!']);
    data=data(:,:,:,:,:,1:NRprocessOnly);
    NR=NRprocessOnly;
end

%average the repetitions
if(NRaverage)
    data=mean(data,6);
    NR=1;
end

%about scan time line (s)
scantime=PVM_RepetitionTime/1000*nVoxels.^2;
timeline=linspace(0,NR*scantime,NR)+scantime/2;

%"à la radiale" reconstruction
if(moreRepetitions>=2)
    %the idea is to recreated some fake acquired kspaces between the existing
    %one by "sliding a window", NR will be changed
    
    %new scan time line
    floatingIndex=1:(1/moreRepetitions):NR;
    NR2=length(floatingIndex);
    timeline=linspace(0,NR*scantime,NR2)+scantime/2;
    
    %init
    data2=data(:,:,:,:,:,1);
    dataOnOneLinePerRepetition=reshape(data,[nPts,nSlices,nVoxels*nVoxels,NR]);
    
    for s=1:nSlices
        for r2=1:NR2
            %first, the base image
            baseImageIndex=floor(floatingIndex(r2));
            tmpImage=dataOnOneLinePerRepetition(:,s,:,baseImageIndex);
            %then let's overwrite some kspace points with the next image
            perctNextImage=floatingIndex(r2)-baseImageIndex;
            nVoxelsFromNextImage=round((nVoxels*nVoxels)*perctNextImage);
            if(nVoxelsFromNextImage>0)
                tmpImage(:,1:nVoxelsFromNextImage)=dataOnOneLinePerRepetition(:,s,1:nVoxelsFromNextImage,baseImageIndex+1);
            end
            %reshaping the crap
            data2(1,:,s,:,:,r2)=reshape(tmpImage,[nPts,nVoxels,nVoxels]);
        end
    end
    disp([num2str(NR2),' images reconstructed from ',num2str(NR),' originally acquired using "a la radial" !']);
    
    %substitute the real data with "time interpolated" data
    data=data2;
    NR=NR2;
end
disp('-----------------------------------------------------------------------------');

%% snail reorganization, spatial FFT

%snail trajectory
nSnailPts=nVoxels*nVoxels-1;
x=zeros(nSnailPts+1,1);
y=zeros(nSnailPts+1,1);
k=2;
for i=1:nSnailPts
    Shell=floor((sqrt(i)+1)/2);
    Leg=floor((i-(2.*Shell-1).^2)./(2.*Shell));
    Element=(i-(2.*Shell-1).^2)-2.*Shell.*Leg-Shell+1;
    
    if(Leg==0)
        y(k)=Shell;
    else
        if(Leg==1)
            y(k)=-Element;
        else
            if(Leg==2)
                y(k)=-Shell;
            else
                y(k)=Element;
            end
        end
    end
    
    if(Leg==0)
        x(k)=Element;
    else
        if(Leg==1)
            x(k)=Shell;
        else
            if(Leg==2)
                x(k)=-Element;
            else
                x(k)=-Shell;
            end
        end
    end
    
    k=k+1;
end

x=x+floor(nVoxels/2)-1;
y=y+floor(nVoxels/2)-1;

%data reorganizing
data_reorg=zeros(nPts,nSlices,nVoxels,nVoxels,NR);
for s=1:nSlices
    for r=1:NR
        k=1;
        for j=1:nVoxels
            for i=1:nVoxels
                data_reorg(:,s,x(k)+1,y(k)+1,r)=squeeze(data(1,:,s,i,j,r));
                k=k+1;
            end
        end
    end
end

% DC bias subtraction
doBiasCorrection = false; % inform Lucio that this was either not helpful or increased the artifact in some cases
if doBiasCorrection
    DC_bias = mean(squeeze(data_reorg(end-20:end,1,:,:)));
    figure,imagesc(abs(squeeze(DC_bias))),colormap jet

    for n=1:size(data_reorg,1)
        DC_for_sub(1,1,:,:) = squeeze(DC_bias);
        data_reorg(n,:,:,:) = data_reorg(n,:,:,:) - DC_for_sub ;
    end
end
    
%spatial zero-filling
nVoxels2=nVoxels+spatialZeroFillingPts;
data_reorg_zf=zeros(nPts,nSlices,nVoxels2,nVoxels2,NR);
data_reorg_zf(:,:,(spatialZeroFillingPts/2+1):(spatialZeroFillingPts/2+nVoxels),(spatialZeroFillingPts/2+1):(spatialZeroFillingPts/2+nVoxels),:)=data_reorg;

%display of the kspace - max of signal
kspaceMax=squeeze(max(abs(data_reorg_zf),[],1));
figure(1),subplot(3,1,1);
imagesc(reshape(permute(kspaceMax,[2,1,3,4]),[nSlices*nVoxels2,nVoxels2*NR]));
axis equal;
axis off;
title(['Kspaces ',num2str(nVoxels2),'x',num2str(nVoxels2),' voxels, ',num2str(nSlices),' slices, ',num2str(NR), ' repetitions - Max of signal']);

%k-space apodization
data_reorg_apo=data_reorg_zf;
spatialApodization = true;
if(spatialApodization)
    % Zinoviy - change this after installing Signal Processing toolbox
    gw=gausswin(nVoxels2);
%     gw = fspecial('gaussian',[nVoxels2 1],3.5);
%     load('C:\Users\Gilf\Dropbox\Reconstruction\Stefan\my_gausswin.mat'); % NOMINAL
%     gw = my_gausswin;
    [maskr,maskc]=meshgrid(gw,gw);
    gw2d=maskr.*maskc;
    for s=1:nSlices
        for r=1:NR
            for t=1:nPts
                data_reorg_apo(t,s,:,:,r)=squeeze(data_reorg_zf(t,s,:,:,r)).*gw2d;
            end
        end
    end
end

%display of the kspace - max of signal
kspaceMax=squeeze(max(abs(data_reorg_apo),[],1));
figure(1),subplot(3,1,2);
imagesc(reshape(permute(kspaceMax,[2,1,3,4]),[nSlices*nVoxels2,nVoxels2*NR]));
axis equal;
axis off;
title(['Apodized Kspaces ',num2str(nVoxels2),'x',num2str(nVoxels2),' voxels, ',num2str(nSlices),' slices, ',num2str(NR), ' repetitions - Max of signal']);

%% spatial 2D FFT
data_reorg_fft2=ifftshift(ifftshift(ifft(ifft(  fftshift(fftshift(  data_reorg_apo  ,3) ,4)  ,nVoxels2,3),nVoxels2,4),3),4);

%display of the image - max of signal
imagMax=squeeze(max(abs(data_reorg_fft2),[],1));
figure(1),subplot(3,1,3);
imagesc(reshape(permute(imagMax,[2,1,3,4]),[nSlices*nVoxels2,nVoxels2*NR]));
axis equal;
axis off;
title(['Images ',num2str(nVoxels2),'x',num2str(nVoxels2),' voxels, ',num2str(nSlices),' slices, ',num2str(NR), ' repetitions - Max of signal']);

%% spectral zero-filling, apodization, spectral FFT
%truncating the first points (digital filter delay)
data_reorg_fft2_trunc=data_reorg_fft2(PVM_DigShift:end,:,:,:,:);

%spectral zero-filling and apodization
nPts2=size(data_reorg_fft2_trunc,1)+spectralZeroFillingPts;
data_reorg_fft2_zf=zeros(nPts2,nSlices,nVoxels2,nVoxels2,NR);
data_reorg_fft2_zf_apo=zeros(size(data_reorg_fft2_zf));
t=linspace(0,nPts2/PVM_SpecSWH,nPts2);

%display before apodization
figure(2),subplot(2,1,1);
%plot(t,squeeze(data_reorg_fft2_trunc(:,:)));
plot(t,squeeze(data_reorg_fft2_zf(:,:)));
hold on;
plot(t,max(abs(squeeze(data_reorg_fft2_trunc(:)))).*exp(-spectralDampingFactor.*t).','k--')
hold off;
grid on;
xlabel('time (s)');
title(['Before spectral apodization - exponential ',num2str(spectralDampingFactor),' Hz']);

data_reorg_fft2_zf(1:size(data_reorg_fft2_trunc,1),:,:,:,:)=data_reorg_fft2_trunc;
for s=1:nSlices
    for r=1:NR
        for j=1:nVoxels2
            for i=1:nVoxels2
                data_reorg_fft2_zf_apo(:,s,i,j,r)=squeeze(data_reorg_fft2_zf(:,s,i,j,r)).*exp(-spectralDampingFactor.*t).';
            end
        end
    end
end

%display after apodization
figure(2),subplot(2,1,2);
plot(t,squeeze(data_reorg_fft2_zf_apo(:,:)));
grid on;
xlabel('time (s)');
title(['After spectral apodization - exponential ',num2str(spectralDampingFactor),' Hz']);

%% spectral FFT
data_reorg_fft3_abs=abs(fftshift(fft(data_reorg_fft2_zf_apo,[],1),1));

% 1. display the sum image of the first frame
sum_img = squeeze(sum(squeeze(data_reorg_fft3_abs(:,:,:,:,1)),1));
SEG_IMAGE = abs(data_proton_cropped(:,:,ProSLAM.DisplayProtonSlice));
figure(18),imagesc(SEG_IMAGE),axis square,colormap gray
axis off
% 2. define the compartmental mask on the 1st time frame from the reconstructed CSI data
NUM = ProSLAM.NumCompartments; % NUMBER OF PRO-SLAM COMPONENTS
N = size(data_proton_cropped,1);
Ni = length(sum_img);

h = fspecial('gaussian',[7 7],0.001); % was 0.001. use stronger smoothing to reduce partial volume when down-sampling
S = size(data_reorg_apo,3);

% in this case the compartment ROIs will be segmented automatically
if (NUM==-1)
    %BW = (SEG_IMAGE) > max(SEG_IMAGE(:))*0.5;
    BW = (SEG_IMAGE) > max(SEG_IMAGE(:))*0.15;
    BW = imfill(BW,'holes');

    [L,NUM] = bwlabeln(BW);
    figure(7),imagesc(L),axis equal
    NUM_C=NUM+1;
    for n = 1:NUM_C
        if (n<=NUM)
            mask = (L==n);
%             se = strel('disk',5);
%             mask = imdilate(mask, se);
        else
            % add the background as another component
            mask = (L==0);
%             se = strel('disk',5);
%             mask = imerode(mask, se);
        end
        
        [xg yg]=meshgrid(linspace(-1,1,N),linspace(-1,1,N));
        [xi yi]=meshgrid(linspace(-1,1,Ni),linspace(-1,1,Ni));
        mask_interp = interp2(xg,yg,single(mask),xi,yi,'linear');
        
        new_L(:,:,n) = mask_interp;
    end
    
    L = new_L;
else
    
    for n = 1:NUM
        % here the user draws the compartment ROI's
        if (ProSLAM.InteractiveContours)
            [x,y] = ginput;
            %[x,y] = smooth_contours(x,y);
            hold on
            save(['c:\temp\contours' num2str(n) '.mat'],'x','y');
        else % or alternatively we can load the ROI's
            load(['c:\temp\contours' num2str(n) '.mat'],'x','y');
        end
        figure(18),hold on,plot([x;x(1)],[y;y(1)],'linewidth',2,'Color','g')
        
        mask=poly2mask(x,y,N,N);
        mask = conv2(single(mask),h,'same');
        
        % interpolate mask to the zero filled C13 matrix size
        [xg yg]=meshgrid(linspace(-1,1,N),linspace(-1,1,N));
        [xi yi]=meshgrid(linspace(-1,1,Ni),linspace(-1,1,Ni));
        mask_interp = interp2(xg,yg,single(mask),xi,yi,'linear');
        if (n==1)
            display_ROI = mask_interp;
        else
            display_ROI = or(display_ROI,mask_interp);
        end
        L(:,:,n) = mask_interp;

        
    end
    if(ProSLAM.InteractiveContours)
        hold off
    end
    NUM_C = NUM+1;
    % add the background as another component
    L(:,:,NUM_C) = 1 - display_ROI;
    L(:,:,NUM_C) = min(L(:,:,NUM_C),1);
    L(:,:,NUM_C) = max(L(:,:,NUM_C),0);
    % display selected regions

end
if(ProSLAM.InteractiveContours)
    hold off
    figure(19),imagesc(display_ROI),axis equal
end
ProSLAM.CSISlice = 1;
SliceNum = ProSLAM.CSISlice;

if ProSLAM.DoSLAM
    
    % spatial 2D FFT
    data_reorg_fft2=ifftshift(ifftshift(ifft(ifft(  fftshift(fftshift(  data_reorg_apo  ,3) ,4)  ,nVoxels2,3),nVoxels2,4),3),4);
    
    % spectral zero-filling, apodization, spectral FFT
    % truncating the first points (digital filter delay)
    if true
        data_reorg_fft2_trunc=data_reorg_fft2(PVM_DigShift:end,:,:,:,:);
    else
        data_reorg_fft2_trunc=data_reorg_fft2(1:end,:,:,:,:);
    end
    
    %spectral zero-filling and apodization
    if true
        nPts2=size(data_reorg_fft2_trunc,1)+spectralZeroFillingPts;
        data_reorg_fft2_zf=zeros(nPts2,nSlices,nVoxels2,nVoxels2,NR);
        data_reorg_fft2_zf_apo=zeros(size(data_reorg_fft2_zf));
        t=linspace(0,nPts2/PVM_SpecSWH,nPts2);
    end
    
    if true
        data_reorg_fft2_zf(1:size(data_reorg_fft2_trunc,1),:,:,:,:)=data_reorg_fft2_trunc;
        for s=1:nSlices
            for r=1:NR
                for j=1:nVoxels2
                    for i=1:nVoxels2
                        data_reorg_fft2_zf_apo(:,s,i,j,r)=squeeze(data_reorg_fft2_zf(:,s,i,j,r)).*exp(-spectralDampingFactor.*t).';
                    end
                end
            end
        end
    else
        data_reorg_fft2_zf_apo=data_reorg_fft2_zf;
    end
    
    if true
        % Spectral FFT
        % If requested, enforce zero phase in the images
        if (ProSLAM.forceZeroPhase)
            data_reorg_fft3_abs = abs(fftshift(fft(data_reorg_fft2_zf_apo,[],1),1));
        else
            data_reorg_fft3_abs = fftshift(fft(data_reorg_fft2_zf_apo,[],1),1);
        end
        
        data_reorg_fft3_abs_no_slam = data_reorg_fft3_abs;
        
    end
    
    %% Main loop of the iterative recon
    for iter = 1:ProSLAM.num_iter
        
        M = 7;
        m = floor(M/2);

        for f=1:size(data_reorg_fft3_abs,5)
            
            if (iter<ProSLAM.num_iter)

                % Replace every pixel with a gaussian weighted average of its peers
                new_data_reorg_fft3_abs = padarray(data_reorg_fft3_abs,[0 0 m m 0]);
                filtered_data_reorg_fft3_abs = new_data_reorg_fft3_abs;
                K = size(new_data_reorg_fft3_abs,1);
                for n = 1:NUM_C
                    
                    % W = fspecial('gaussian',[M M],2); % sigma was 1
                    W = fspecial('gaussian',[M M],ProSLAM.sigmaGaussian);
                    if (n==NUM_C)
                        W = fspecial('gaussian',[M M],ProSLAM.sigmaGaussian); % use sigma=8 for the thermal phantom
                    end
                    W = repmat(W,[1 1 K]);
                    W = permute(W,[3,1,2]);
                    
                    mask = single(squeeze(L(:,:,n)));
                    mask = padarray(mask,[m m]);
                    mask = repmat(mask,[1 1 K]);
                    mask=permute(mask,[3,1,2]);
                    
                    % filter (smooth) the magnitude image
                    for x = 1+m:size(new_data_reorg_fft3_abs,4)-m
                        for y = 1+m:size(new_data_reorg_fft3_abs,3)-m
                            
                            if (mask(1,y,x)>0)
                               MyPeerAbs = squeeze(abs(new_data_reorg_fft3_abs(:,SliceNum,y-m:y+m,x-m:x+m,f))).* mask(:,y-m:y+m,x-m:x+m) .*W;

                               peers =  (squeeze(sum(squeeze(sum(MyPeerAbs,2)),2))) ./ squeeze(sum(squeeze(sum(mask(:,y-m:y+m,x-m:x+m) .* W,2)),2));
                               
                                filtered_data_reorg_fft3_abs(:,SliceNum,y,x,f) = peers ;

                            end
                        end
                    end
                end
                data_reorg_fft3_abs=filtered_data_reorg_fft3_abs(:,:,1+m:end-m,1+m:end-m,:);
            else % last iteration
                figure(9),subtightplot(1,10,f),imagesc(squeeze(abs(sum(data_reorg_fft3_abs(:,SliceNum,:,:,f),1)))),colormap jet,axis square,axis off
            end
            
            % set central k-space back to the experimental measurements
            for s=1:size(data_reorg_fft3_abs,1)

                orig_img = squeeze(data_reorg_fft3_abs_no_slam(s,SliceNum,:,:,f));
                orig_ksp = fftshift(fft2(orig_img));
                
                this_img = squeeze(data_reorg_fft3_abs(s,SliceNum,:,:,f));
                this_img = abs(this_img) .*exp(1i*angle(orig_img)); % we take only magnitude changes in the image
                this_ksp = fftshift(fft2(this_img));

                doTapering = true;
                if doTapering
                    % W = padarray(my_hanning(nVoxels)',floor((nVoxels2-nVoxels)/2)); % use hanning(nVoxels).^2 if you really trust the filter
                    % GF changed on 04/01/2023
                    W = padarray(hann(nVoxels)',floor((nVoxels2-nVoxels)/2));
                    W = repmat(W,[1 length(this_ksp)]);
                    W = W.*W';
                    combined_ksp = W.*orig_ksp + (1-W).*this_ksp;
                end
                
                data_reorg_fft3_abs(s,SliceNum,:,:,f) = ifft2(combined_ksp);
            end

        end
        disp(['iter ' num2str(iter)])      
    end
    
end

%rotate and flip if needed
if(isequal(PVM_SPackArrSliceOrient,'axial'))
    data_reorg_fft3_abs_rotated=data_reorg_fft3_abs;
    for t=1:nPts2
        for s=1:nSlices
            for r=1:NR
                data_reorg_fft3_abs(t,s,:,:,r)=rot90(squeeze(data_reorg_fft3_abs_rotated(t,s,:,:,r)),-1);
            end
        end
    end
    data_reorg_fft3_abs=data_reorg_fft3_abs(:,:,:,end:-1:1,:);
end

%% quantification maps

if (ProSLAM.DoSLAM)

    data_reorg_fft3_abs = abs(data_reorg_fft3_abs);
                        
end

if (ProSLAM.interactiveCS)
    % MyVec = real(data_reorg_fft3_abs(:,1:1:end));
    MyVec = real(data_reorg_fft3_abs(:,1:5:end));
    figure,plot(MyVec)
    title('Select signal range! (2 clicks)');
    [x,y]=ginput(2);
    x=sort(x);
    x(1)=round(x(1));
    x(2)=round(x(2));
else
    load('c:\temp\spectral_win.mat','x','b');
end

%test CSdisplacement
[~,ipeak]=max(max(data_reorg_fft3_abs(x(1):x(2),:),[],2),[],1);
f=linspace(-PVM_SpecSWH/2,+PVM_SpecSWH/2,nPts2);
BF1_offset=f(ipeak+x(1));%Hz
Gss=ACQ_trim(1,1);
GssHzmm=Gss/100*PVM_GradCalConst;
CSdisplacement=BF1_offset/GssHzmm;
disp(['INFO: CS displacement for a ',num2str(BF1_offset),' Hz offset in frequency = ',num2str(CSdisplacement),' mm']);

if (ProSLAM.interactiveCS)
    title('Select noise range! (2 clicks)');
    [b,y]=ginput(2);
    close;
    b=sort(b);
    b(1)=round(b(1));
    b(2)=round(b(2));
end

if (ProSLAM.interactiveCS)
    save('C:\Users\Owner\Desktop\PRO-SLAM\PRO-SLAM\ForZinoviy\Stefan\temp\temp.mat','x','b');
end
%% Measurement of metabolite kinetics
if (false)
    kinetics = zeros([1 length(timeline)]);
    % GF - plot the kinetics of the metabolite
    for f=1:size(data_reorg_fft3_abs,5)
        TmpDataset = data_reorg_fft3_abs(x(1):x(2),:,:,:,f);
        kinetics(f) = abs(sum(TmpDataset(:)));
    end
    figure,plot(timeline,kinetics,'LineWidth',2.5)
else
    kinetics = zeros([NUM_C-1 length(timeline)]);
    TmpImage = zeros(size(data_reorg_fft3_abs));
    % claculate kinetics curve per ROI
    for n = 1:NUM_C-1
        mask_BW = squeeze(L(:,:,n));
        K = size(data_reorg_fft3_abs,1);
        mask = repmat(mask_BW,[1 1 K]);
        mask=permute(mask,[3,1,2]);
        % iterate over time-reps
        for f=1:size(data_reorg_fft3_abs,5)
            
            for xx = 1:size(data_reorg_fft3_abs,4)
                for yy = 1:size(data_reorg_fft3_abs,3)
                    
                    if (mask(1,yy,xx))
                        removeBackground = true;
                        if (removeBackground)
                            baseLine = mean(squeeze((data_reorg_fft3_abs(b(1):b(2),SliceNum,yy,xx,f))));
                         
                            TmpVec = (squeeze(data_reorg_fft3_abs(x(1):x(2),SliceNum,yy,xx,f)) - baseLine).* mask(x(1):x(2),yy,xx);
                            TmpVec(real(TmpVec)<0)=0;
                            TmpImage(x(1):x(2),SliceNum,yy,xx,f) = TmpVec ;
    
                        else
                            TmpImage(x(1):x(2),SliceNum,yy,xx,f) = squeeze((data_reorg_fft3_abs(x(1):x(2),SliceNum,yy,xx,f))).* mask(x(1):x(2),yy,xx) ;
                        end
                        
                    else
                        TmpImage(:,SliceNum,yy,xx,f) = 0;
                    end
                end
            end
            
            % normalize to the number of pixels in the compartment
            kinetics(n,f) = sum(sum(sum(TmpImage(:,SliceNum,:,:,f)))) / sum(mask_BW(:));    
            % calculate the max pixel in this component
            max_kinetics(n,f) = max(max(max(TmpImage(:,SliceNum,:,:,f)))) ;  
            % calculate signal for the COM pixel
            
        end % loop over repetitions
    end % loop over ROI's
    
    % calculate noise maps
    calcSNRfactor = ProSLAM.calcSNRfactor; % SET DoSLAM to true if you want to see the noise maps
    if calcSNRfactor
        
        for f=1:1
            for xx = 1:size(data_reorg_fft3_abs,4)
                for yy = 1:size(data_reorg_fft3_abs,3)
                    noise_map(yy,xx) = std(squeeze(abs(data_reorg_fft3_abs(1:50,SliceNum,yy,xx)))); % was 30
                    noise_map_ft(yy,xx) = std(squeeze(abs(data_reorg_fft3_abs_no_slam(1:50,SliceNum,yy,xx))));
                    
                    sig_map(yy,xx) = mean(squeeze(abs(data_reorg_fft3_abs(x(1):x(2),SliceNum,yy,xx))));
                    sig_map_ft(yy,xx) = mean(squeeze(abs(data_reorg_fft3_abs_no_slam(x(1):x(2),SliceNum,yy,xx))));
                end
            end
            if true
                max_noise_FT = max(noise_map_ft(:));
                min_noise_SLAM = min(noise_map(:));
                figure(11),imagesc(squeeze(noise_map(:,:)),[min_noise_SLAM floor(max_noise_FT)]),colormap jet,axis square,axis off,colorbar('FontSize',16)
                figure(12),imagesc(squeeze(sig_map(:,:))),colormap jet,axis square,axis off,colorbar('FontSize',16)
                max_sigmap = max(sig_map(:));
                figure(13),imagesc(squeeze(noise_map_ft(:,:)),[min_noise_SLAM floor(max_noise_FT)]),colormap jet,axis square,axis off,colorbar('FontSize',16)
                figure(14),imagesc(squeeze(sig_map_ft(:,:)),[0 max_sigmap]),colormap jet,axis square,axis off,colorbar('FontSize',16)
                figure(15),imagesc(squeeze(sig_map(:,:) .* noise_map_ft(:,:) ./ (noise_map(:,:) .* sig_map_ft(:,:))),[0 4]),colormap jet,axis square
                set(gca,'YTick','','XTick','');
                figure(16),plot(squeeze(abs(data_reorg_fft3_abs(:,SliceNum,34,15))),'LineWidth',3,'Color','k'),set(gca,'FontSize',16)
                set(gca,'YTick','','XTick','');

                %figure(10),plot(squeeze(abs(data_reorg_fft3_abs_no_slam(:,SliceNum,20:25,20,f))),'LineWidth',3,'Color','k'),set(gca,'FontSize',16)
                figure(10),plot(squeeze(abs(data_reorg_fft3_abs_no_slam(:,SliceNum,1:25,20,f))),'LineWidth',2),set(gca,'FontSize',16)
                set(gca,'YTick','','XTick','');

            end

        end

    end
end

ProSLAM.ExcelFile = ['c:\temp\kinetics_' date '_' num2str(floor(rand()*100))];
if ProSLAM.DoSLAM
    xlswrite(ProSLAM.ExcelFile,[timeline;kinetics],'SLAM');
else
    xlswrite(ProSLAM.ExcelFile,[timeline;kinetics],'FT');
end
%xlswrite('c:\temp\kinetics',cluster_count,'cluster');
xlabel('time [s]','FontSize',16), ylabel('Integrated Signal','FontSize',16)
set(gca,'FontSize',16)

%snr
[s,~]=max(abs(squeeze(data(1,:,1,1))));
spectrumToAnalyse=squeeze(data_reorg_fft3_abs(b(1):b(2),:));
moyB=mean(spectrumToAnalyse(:));
disp(['INFO: SNR (1st k-space pt) = ',num2str(s/moyB)]);

peakArea=zeros(nSlices,nVoxels2,nVoxels2,NR);
peakArea2=zeros(nSlices,nVoxelsQuantificationMaps,nVoxelsQuantificationMaps,NR);
for s=1:nSlices
    for r=1:NR
        for j=1:nVoxels2
            for i=1:nVoxels2
                %selecting spectrum to quantify
                spectrumToAnalyse=squeeze(data_reorg_fft3_abs(:,s,i,j,r));
                %calculating noise level (mean) of the noise region
                moyB=mean(spectrumToAnalyse(b(1):b(2)));
%                 moyB = 0; % GF CHECK THAT Jan 2023
                if(mean(spectrumToAnalyse(x(1):x(2)))>moyB)
                    %selecting signal region and remove noise level
                    spectrumToAnalyse=spectrumToAnalyse(x(1):x(2))-moyB;
                    spectrumToAnalyse(spectrumToAnalyse<0)=0;
                    %integrate
                    peakArea(s,i,j,r)=trapz(spectrumToAnalyse);
                else
                    %if no signal, put to zero (removing noise)
                    peakArea(s,i,j,r)=0;
                end
            end
        end
        %smoothing/interpolation
        peakArea2(s,:,:,r)=imresize(squeeze(peakArea(s,:,:,r)),[nVoxelsQuantificationMaps nVoxelsQuantificationMaps],'bilinear');
    end
end

%display quantification maps
figure(3),subplot(2,1,1);
imagesc(reshape(permute(peakArea,[2,1,3,4]),[nSlices*nVoxels2,nVoxels2*NR]));
axis equal;
axis off;
title(['Quantification maps (peak area) ',num2str(nVoxels2),'x',num2str(nVoxels2),' voxels (raw resolution), ',num2str(nSlices),' slices, ',num2str(NR), ' repetitions']);

figure(3),subplot(2,1,2);
imagesc(reshape(permute(peakArea2,[2,1,3,4]),[nSlices*nVoxelsQuantificationMaps,nVoxelsQuantificationMaps*NR]));
axis equal;
axis off;
title(['Interpolated quantification maps (peak area) ',num2str(nVoxelsQuantificationMaps),'x',num2str(nVoxelsQuantificationMaps),' voxels, ',num2str(nSlices),' slices, ',num2str(NR), ' repetitions']);

%% fitting time curves
disp('INFO: fitting the time curves...');
diffusionCoeff=zeros(nSlices,nVoxelsQuantificationMaps,nVoxelsQuantificationMaps);
for s=1:nSlices
    mslice=max(peakArea(s,:));
    for j=1:nVoxelsQuantificationMaps
        for i=1:nVoxelsQuantificationMaps
            tc=squeeze(peakArea2(s,i,j,:));
            [p,S]=polyfit(1:NR,tc.',1);
            diffusionCoeff(s,i,j)=p(1);
        end
    end
end

%display diffusion maps
figure(4);
imagesc(reshape(permute(diffusionCoeff,[2,1,3]),[nSlices*nVoxelsQuantificationMaps,nVoxelsQuantificationMaps]));
axis equal;
axis off;
colorbar;
title(['Linear fit of time domain curves ',num2str(nVoxelsQuantificationMaps),'x',num2str(nVoxelsQuantificationMaps),' voxels, ',num2str(nSlices),' slices']);

%% overlaying with anatomic images

%checking
if(length(anatomicImageFilenames)~=nSlices)
    error('Mmmh... The number of slices acquired in 13C does not match the number of 1H images provided!');
end
%cd to images folder
cd(anatomicImagesFolder);

%max for this dataset (normalization)
maxForAll=max(squeeze(peakArea2(:)));
maxForAllRaw=max(squeeze(peakArea(:)));
for s=1:nSlices
    %read anatomic image
    protonImg=mean(double(importdata(anatomicImageFilenames{s})),3);
    %making it really square
    nVoxelsProtonImg=size(protonImg,1);
    protonImg=imresize(protonImg,[nVoxelsProtonImg nVoxelsProtonImg]);
    %normalize it and convert it to a full rgb vector
    protonImg=protonImg-min(protonImg(:));
    protonImg=protonImg./max(protonImg(:));
    protonImg=repmat(protonImg,[1,1,3]);
    %max for this slice (normalization)
    maxForThisSlice=max(squeeze(peakArea2(s,:)));
    maxForThisSliceRaw=max(squeeze(peakArea(s,:)));
    for r=1:NR
        %normalize the carbon image
        %0=do nothing
        %1=normalize - repetitions (at least one red voxel on each repetitions)
        %2=normalize - slices (at least one red voxel for each slice)
        %3=normalize - all the data (at least one red voxel)
        switch(normalizeQuantificationMaps)
            case 0
                carbonImg=squeeze(peakArea2(s,:,:,r));
                carbonImgRaw=squeeze(peakArea(s,:,:,r));
            case 1
                carbonImg=squeeze(peakArea2(s,:,:,r))./max(max(squeeze(peakArea2(s,:,:,r))));
                carbonImgRaw=squeeze(peakArea(s,:,:,r))./max(max(squeeze(peakArea(s,:,:,r))));
            case 2
                carbonImg=squeeze(peakArea2(s,:,:,r))./maxForThisSliceRaw;
                carbonImgRaw=squeeze(peakArea(s,:,:,r))./maxForThisSliceRaw;
            case 3
                carbonImg=squeeze(peakArea2(s,:,:,r))./maxForAll;
                carbonImgRaw=squeeze(peakArea(s,:,:,r))./maxForAllRaw;
        end
        %no negative value (!)
        carbonImg(carbonImg<0)=0;
        %resize to proton image size without smoothing effect
        carbonImg=imresize(carbonImg,[nVoxelsProtonImg nVoxelsProtonImg],'nearest');
        carbonImgRaw=imresize(carbonImgRaw,[nVoxelsProtonImg nVoxelsProtonImg],'nearest');
        %store
        protonImgs(s,:,:,r,:)=protonImg;
        carbonImgs(s,:,:,r)=carbonImg;
        carbonImgsRaw(s,:,:,r)=carbonImgRaw;
    end
end
cd(scanFolder);

%spatial shift calculation for 13C/1H overlaying
shiftXY=-0.5*nVoxelsProtonImg/nVoxels2;

figure(5);
imagesc(reshape(permute(protonImgs,[2,1,3,4,5]),[nSlices*nVoxelsProtonImg,nVoxelsProtonImg*NR,3]));
hold on;
h=imagesc((1:(nVoxelsProtonImg*NR))+shiftXY,(1:(nSlices*nVoxelsProtonImg))+shiftXY,reshape(permute(carbonImgsRaw,[2,1,3,4]),[nSlices*nVoxelsProtonImg,nVoxelsProtonImg*NR]));
hold off
adata=reshape(permute(carbonImgsRaw,[2,1,3,4]),[nSlices*nVoxelsProtonImg,nVoxelsProtonImg*NR]);
set(h,'AlphaData',adata);
colormap jet
axis equal;
axis off;
%colorbar;
title(['Maps overlayed on anatomic images ',num2str(nVoxelsProtonImg),'x',num2str(nVoxelsProtonImg),' voxels, ',num2str(nSlices),' slices, ',num2str(NR), ' repetitions']);

figure(6);
subtightplot(2,1,1),imagesc(reshape(permute(protonImgs,[2,1,3,4,5]),[nSlices*nVoxelsProtonImg,nVoxelsProtonImg*NR,3]));
colormap jet
axis equal;
axis off;

subtightplot(2,1,2)
imagesc(reshape(permute(protonImgs,[2,1,3,4,5]),[nSlices*nVoxelsProtonImg,nVoxelsProtonImg*NR,3]));
hold on;
h=imagesc((1:(nVoxelsProtonImg*NR))+shiftXY,(1:(nSlices*nVoxelsProtonImg))+shiftXY,reshape(permute(carbonImgs,[2,1,3,4]),[nSlices*nVoxelsProtonImg,nVoxelsProtonImg*NR]));
hold off
adata=reshape(permute(carbonImgs,[2,1,3,4]),[nSlices*nVoxelsProtonImg,nVoxelsProtonImg*NR]);
set(h,'AlphaData',adata);
axis equal;
axis off;
FigHandle = figure(6)
%set(h, 'Position', [100, 100, 3*390, 3*160]);

% GF - my part ends here
toc
return;

