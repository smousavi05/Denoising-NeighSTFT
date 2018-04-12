
clc
clear all
close all 

%% Read the original 
load d7
data.x=d7.l1400;
data.dt=0.005;
data.t = linspace(0,length(data.x)*data.dt,length(data.x));
data.x = data.x/max(data.x);

tic
% block thresholding 
[denoised] = neigBlock(data);
% plot(real(denoised))
toc

Xnoisy = data.x/max(data.x);
denoised_norm = denoised/max(denoised); 
% denoised_norm = denoised_norm*-1; 
denoised_norm = real(denoised_norm);


figure(1)
subplot(2,3,1);plot(data.t,data.x);
axis tight
title('Recorded Signal','Rotation',0,'FontSize',14);
xlabel({'Time (s)'},'FontSize',12); 
ylabel({'Normalized Amplitude'},'FontSize',12); 
xlim([min(data.t) max(data.t)]);
grid on, grid minor 
ax = gca;
ax.TitleFontSizeMultiplier = 2.1;
ax.LabelFontSizeMultiplier=2.1;
ax.FontWeight='bold';
hold off

subplot(2,3,4);
time_win = 1500;
factor_redund = 1;
stftCoef = STFT(data.x, time_win, factor_redund, 1/data.dt);
[a b] = size(stftCoef);
imagesc(abs(stftCoef(1:60,:)));
% colormap(jet)
clim=get(gca,'clim');
% clim=[0 35];
set(gca,'clim',clim)

title('STFT','Rotation',0,'FontSize',14);
xlabel({'Time (s)'},'FontSize',12); 
ylabel({'Frequency (Hz)'},'FontSize',12); 
ax.TitleFontSizeMultiplier = 2.1;
ax.LabelFontSizeMultiplier=2.1;
ax.FontWeight='bold';
hold off


subplot(2,3,2);
waveform = data.x;
sampling_rate = 200; 
highpass_freq = 2; 
lowpass_freq = 15; 
order = 2;
nyquist_freq = sampling_rate/2;                   % Nyquist frequency
Wn = [highpass_freq, lowpass_freq]/nyquist_freq;  % non-dimensional frequency
[b,a] = butter(order/2, Wn, 'bandpass');  % construct the filter
xxf = filtfilt(b, a, waveform);          % filter the data with zero phase

plot(data.t,xxf)
axis tight
grid on; grid minor
title('Band-Pass Filtered','Rotation',0,'FontSize',14);
xlabel({'Time (s)'},'FontSize',12);
ylabel({'Normalized Amplitude'},'FontSize',12);
xlim([min(data.t) max(data.t)]);
ax.LabelFontSizeMultiplier=2.1;
ax.FontWeight='bold';
hold off

subplot(2,3,5);
stftCoef = STFT(xxf, time_win, factor_redund, 1/data.dt);
[a b] = size(stftCoef);
for i=1:2;
stftCoef(i,:)=0;
end
for i=20:a;
stftCoef(i,:)=0;
end
imagesc(abs(stftCoef(1:60,:)));
% colormap(jet)
clim=get(gca,'clim');
% clim=[0 10];
set(gca,'clim',clim)

title('STFT','Rotation',0,'FontSize',14);
xlabel({'Time (s)'},'FontSize',12); 
ylabel({'Frequency (Hz)'},'FontSize',12); 
ax.TitleFontSizeMultiplier = 2.1;
ax.LabelFontSizeMultiplier=2.1;
ax.FontWeight='bold';
hold off


subplot(2,3,3);
plot(data.t,denoised_norm(1:length(data.t)))
axis tight
grid on; grid minor
title('Denoised Signal','Rotation',0,'FontSize',14);
xlabel({'Time (s)'},'FontSize',12);
ylabel({'Normalized Amplitude'},'FontSize',12);
xlim([min(data.t) max(data.t)]);
ax.LabelFontSizeMultiplier=2.1;
ax.FontWeight='bold';
hold off

subplot(2,3,6);
stftCoef = STFT(denoised_norm, time_win, factor_redund, 1/data.dt);
[a b] = size(stftCoef);
imagesc(abs(stftCoef(1:60,:)));
% colormap(jet)
clim=get(gca,'clim');
% clim=[0 10];
set(gca,'clim',clim)

title('STFT','Rotation',0,'FontSize',14);
xlabel({'Time (s)'},'FontSize',12); 
ylabel({'Frequency (Hz)'},'FontSize',12); 
ax.TitleFontSizeMultiplier = 2.1;
ax.LabelFontSizeMultiplier=2.1;
ax.FontWeight='bold';
hold off




 
