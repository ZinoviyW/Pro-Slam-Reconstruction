function [proSLAM_IMG] = recon_algo(numIter,anatL,ncomp,sigmaGaussian,orig_img,orig_ksp)


%RECON_ALGO PRO SLAM reconstruction
%   
    for iter = 1:numIter
        %gaussian compartmental smoothing 
        for j=1:ncomp
            W = fspecial('gaussian',[5 5],sigmaGaussian);
            smooth_img=roifilt2(W,abs(orig_img),squeeze(anatL(:,:,j)));
        end
        
        %restore k-space data points to experimental values
        smooth_img_ksp = abs(smooth_img) .*exp(1i*angle(orig_img));
        smooth_img_ksp = fftshift(fft2(smooth_img_ksp));
        W = padarray(hann(32)',16);
        W = repmat(W,[1 length(smooth_img_ksp)]);
        W = W.*W';
        combined_ksp = W.*orig_ksp + (1-W).*smooth_img_ksp;  
        I = ifft2(combined_ksp);    
    end
proSLAM_IMG=abs(I);
