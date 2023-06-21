p = phantom(512);
freq = fft2(p);
absfreq = abs(freq);
subplot(2,1,1)
figure
plot(abs(freq))
figure
imagesc(p)

