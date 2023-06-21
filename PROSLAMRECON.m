close force all;
clc;
clear variables;


%% The data
% Anatomic 1H images (batch7_mouse3_anatomic_1H.mat - "anaimage" is collected before injection,
% "anaimageEND" is after injection):
% RARE sequence 512 x 512, 10 slices 40 mm x 40 mm, 0.8mm thickness each
% 
% the separated 2H images (batch7_mouse3_SepImages.mat - "SepImages" is 64 x 64 x 26 x 5):
% 64 x 64, 40 mm x 40 mm, thickness whole animal
% 26 slices (same spatial position in function of 26 point in time)
% dimension 4: 1 - water; 2 - glucose; 3 - lactate; 4 - b0 map; 5 - T2 map.



%%Reconstruction Folder
% The metabolic data folder
recon_data = ['/home/zinoviy/Desktop/ZInoviy/PRO-SLAM/' ...
    'PRO-SLAM/ForZinoviy/Stefan/Eltons_Data/Recon_Folder'];
%T he anatomical data folder
anatData = ['/home/zinoviy/Desktop/ZInoviy/PRO-SLAM/PRO-SLAM/' ...
    'ForZinoviy/Stefan/Eltons_Data/Recon_Folder/Anatomical_Data'];
cd(anatData)
load('batch7_mouse3_anatomic_1H.mat')
% The anatomical 512x512 image
anatImg = anaimage(:,:,5);

%Superpixel segmentation
[L,numLabels] = superpixels(anatImg,70,"NumIterations",20,"Compactness",0.5,"IsInputLab",false);
%Display the original segmentation
figure
BW = boundarymask(L);
imshow(imoverlay(anatImg,BW,'cyan'),[])

% Fixing the superpixel segmentation
means = double.empty(numLabels,0);
for i=1:numLabels
    means(i) = mean(nonzeros(double((L==i)).*anatImg));
end


for j=1:length(means)
    if(means(j)<=1)
        L (L==j) = 1;
    end
end

n = length(nonzeros(means<=1));
newL = length(means) - n+1;
i = 2;
while(max(max(L))~=newL)
    m = max(max(L));
    L(L==m)=i;
    i=i+1;
end

numLabels = newL;

%Display the new segmentation
figure
BW = boundarymask(L);
imshow(imoverlay(anatImg,BW,'cyan'),'InitialMagnification',67)
cd(recon_data)

%Create binary segmentation masks
new_L=zeros([64 64 numLabels]);
    for n = 1:numLabels
        if (n<=numLabels)
            mask = (L==n);
        else
            % add the background as another component
            mask = (L==0);
        end
        
        [xg, yg]=meshgrid(linspace(-1,1,512),linspace(-1,1,512));
        [xi, yi]=meshgrid(linspace(-1,1,64),linspace(-1,1,64));
        mask_interp = interp2(xg,yg,single(mask),xi,yi,'linear');
        
        new_L(:,:,n) = mask_interp;
    end

L = new_L;
cd(recon_data)
% Load ELtons' separated images
load("mouse3_sep_images.mat")
recon_data = SepImages;% data [64 64 26 5] d*d of 25 time pts for 5 metabolites

%%Go for the reconstruction
tic
for d=1:5
    for s=1:26
        orig_img = recon_data(:,:,s,d);
        orig_img_abs=abs(orig_img);
        orig_ksp = fftshift(fft2(orig_img_abs));
        % The reconstruction function
        recon_image = recon_algo(5,L,numLabels,1.5,orig_img,orig_ksp);
        cd(['/home/zinoviy/Desktop/ZInoviy/PRO-SLAM/PRO-SLAM/' ...
            'ForZinoviy/Stefan/Eltons_Data/Recon_Folder/Recon_Data'])
        save(convertStringsToChars(string(d)+"_"+string(s)),'recon_image');
    end
end


cd("/home/zinoviy/Desktop/ZInoviy/PRO-SLAM/PRO-SLAM/ForZinoviy/Stefan/Eltons_Data")
load("SepImages.mat")
cd('/home/zinoviy/Desktop/ZInoviy/PRO-SLAM/PRO-SLAM/ForZinoviy/Stefan/Eltons_Data/Recon_Folder/Recon_Data')


%% Displays the results
water_figure1 = figure("Name","Water - After Reconstraction");
title("Water - After Reconstraction")
w = tiledlayout(6,9,"TileSpacing","Compact");


for i=1:26
    % Display after recon
    w_pr = tiledlayout(w,1,1);
    w_pr.Layout.Tile = i;
    w_pr.Layout.TileSpan = [1 1];
    nexttile(w_pr);
    load("1_"+string(i)+".mat");
    imagesc(recon_image)
    axis off
    title("t="+string(i))
    
    
end

water_figure2 = figure("Name","Water - Before Reconstraction");
title("Water - Before Reconstraction")
w2 = tiledlayout(6,9,"TileSpacing","Compact");
for i=1:26
    % Display before recon
    w_unpr = tiledlayout(w2,1,1);
    w_unpr.Layout.Tile = i;
    w_unpr.Layout.TileSpan = [1 1];
    nexttile(w_unpr);
    imagesc(abs(SepImages(:,:,i,1)))
    axis off
    title("t="+string(i))

end

glucose_figure = figure("Name","Glucose");

g = tiledlayout(6,9);


for i=1:26
    % Display after recon
    
    w_pr = tiledlayout(g,1,1);
    w_pr.Layout.Tile = i;
    w_pr.Layout.TileSpan = [1 1];
    nexttile(w_pr);
    load("2_"+string(i)+".mat");
    imagesc(recon_image)
    axis off
    title("t="+string(i))
    
    
end


for i=1:26
    % Display before recon
    w_unpr = tiledlayout(g,1,1);
    w_unpr.Layout.Tile = i+26;
    w_unpr.Layout.TileSpan = [1 1];
    nexttile(w_unpr);
    imagesc(abs(SepImages(:,:,i,2)))
    axis off
    title("t="+string(i))

end

lactate_figure = figure("Name","Lactate");

l = tiledlayout(6,9);


for i=1:26
    % Display after recon
    
    w_pr = tiledlayout(l,1,1);
    w_pr.Layout.Tile = i;
    w_pr.Layout.TileSpan = [1 1];
    nexttile(w_pr);
    load("3_"+string(i)+".mat");
    imagesc(recon_image)
    axis off
    title("t="+string(i))
    
    
end


for i=1:26
    % Display before recon
    w_unpr = tiledlayout(l,1,1);
    w_unpr.Layout.Tile = i+26;
    w_unpr.Layout.TileSpan = [1 1];
    nexttile(w_unpr);
    imagesc(abs(SepImages(:,:,i,3)))
    axis off
    title("t="+string(i))

end

b0_figure = figure("Name","B0");

b = tiledlayout(6,9);


for i=1:26
    % Display after recon
    
    w_pr = tiledlayout(b,1,1);
    w_pr.Layout.Tile = i;
    w_pr.Layout.TileSpan = [1 1];
    nexttile(w_pr);
    load("4_"+string(i)+".mat");
    imagesc(recon_image)
    axis off
    title("t="+string(i))
    
    
end


for i=1:26
    % Display before recon
    w_unpr = tiledlayout(b,1,1);
    w_unpr.Layout.Tile = i+26;
    w_unpr.Layout.TileSpan = [1 1];
    nexttile(w_unpr);
    imagesc(abs(SepImages(:,:,i,4)))
    axis off
    title("t="+string(i))

end

t2_figure = figure("Name","T2");

s = tiledlayout(6,9);


for i=1:26
    % Display after recon
    
    w_pr = tiledlayout(s,1,1);
    w_pr.Layout.Tile = i;
    w_pr.Layout.TileSpan = [1 1];
    nexttile(w_pr);
    load("5_"+string(i)+".mat");
    imagesc(recon_image)
    axis off
    title("t="+string(i))
    
    
end


for i=1:26
    % Display before recon
    w_unpr = tiledlayout(s,1,1);
    w_unpr.Layout.Tile = i+26;
    w_unpr.Layout.TileSpan = [1 1];
    nexttile(w_unpr);
    imagesc(abs(SepImages(:,:,i,5)))
    axis off
    title("t="+string(i))

end


% % size of test image
% sz = [64 64];
% t = 1:26;
% setnames = {'apple','banana','cherry','pear','orange'};
% 
% tiling = [5 5]; % number of plots per page [rows cols]
% fontsize = 10; 
% 
% htl = tiledlayout(tiling(1),tiling(2));
% ht1.TileSpacing = 'none'; % set the tile spacing as tight as you want
% for tidx = 1:tiling(1)
%     % index through the images in each set
%     for setidx = 1:tiling(2)
%         % index between sets
%         nexttile
%         thisdata = testdata(:,:,tidx,setidx);
%         imagesc(thisdata)
%         axis off
%         
%         % put labels only on the sides
%         if setidx == 1
%             ht = text(0,0,sprintf('t = %d',t(tidx)));
%             set(ht,'horizontalalignment','center','rotation',90,'units','normalized', ...
%                 'position',[-0.2 0.5],'fontweight','bold','fontsize',fontsize);
%         end
%         if tidx == 1
%             ht = text(0,0,setnames{setidx});
%             set(ht,'horizontalalignment','center','rotation',0,'units','normalized', ...
%                 'position',[0.5 1.15],'fontweight','bold','fontsize',fontsize);
%         end
%     end
% end


toc