%{
Title: Assignment 2
Course Code: AESM1511
Name: Karla Echeverry Caro
Studentnumber: 4537270
Date Created: 04 October 2022
Date modified: 04 October 2022
Mail: K.V.EcheverryCaro@student.tudelft.nl
Assignment 2:  Fourier Transformation and Filtering
%}
%% Initial Settings

clc;clear;close all;
%% Initial Inputs

%Velocities
Velp= [1200, 1800, 3000];% Velocities of the P waves
Vels= [500, 1000, 1800]; % Velocities of the S waves

%Traces information
Recevspace=2.5; %Spacing of the receivers
FRecev= 510; %Start Location of the receiver
ERecev= 1510; %End Location of the receiver
samples=401;% Amount of traces
Distance=FRecev:Recevspace:ERecev; %Horizontal Distance
wavenumber=1:1:samples;

%Time information
Timesam=1e-3; % 1 microsecond sampling time;
timesamples=1001; %Total sampled time

% calculate time frequency sampling
tf_s = 1/(timesamples*Timesam);

% calculate space frequency sampling
sf_s = 1/(samples*Recevspace);

%Densities are ignored for this assignment

%Filter
maxfreq=200;
minfreq=0;

%Low Cut Filter Parameters
lcf1max=60; %MAX Frequency
lcf1min=30; %MIN Frequency
lcf1d=185; %Displacement 

%Path and name of the file
path=['C:\Users\khato\OneDrive\Desktop\Signal & Analysis\Matlab\AESM1511_2022' ...
    '_A2_EcheverryCaro_4537270'];%path
name='refl_3layers_fp50_dx0p5_500_rvz.bin';% name of the file
parameter=append(path,'\',name);%combine path and name

%Additional Information
%{
This Code uses a function created by me and it is called plotcreator.m

The code can be found below:

function plotcreator(data,cmap, xlab, ylab, ti, xvalue)
% plotcreator plots imagesc with specified colormap, x-label, y-label and title.


imagesc(data,'xData',xvalue);    % plots data
set(gca,'linewidth',2); % sets linewidth
colormap(cmap);         % set colormap to parameter
colorbar;               % plot colorbar
t = title(ti);          % set title to parameter
lx = xlabel(xlab);      % set xlabel to parameter
ly = ylabel(ylab);      % set ylabel to parameter
t.FontSize = 14;        % set font sizes
lx.FontSize = 14;
ly.FontSize = 14;
end
%}
%% Background Information
%{
Ricker Wavelet:
Negative normalized second derivative of a Gaussian gunction up to
scale and normalisation the second hermite function.
%}

%% Read Bin File
fid=fopen(parameter,'r');
temp=fread(fid,samples*timesamples,'float32',0,'ieee-le');
fclose(fid);
data_refl(:,:)=reshape(temp,timesamples,samples);


%% Bin File Visualization

figure();
plotcreator(data_refl,gray,'Horizontal Distance [m]','Two-way traveltime [s]','Common Source-Gather',Distance)

%% Transformation from Time to Frequency

% This is done using Fast Fourier Transform 
% fft(data in time-space domain,[],?)*(time sampling)
fftdata1=fft(data_refl,[],1)*(Timesam); % FFT from all the columns
fftdata2=fft(data_refl,[],2)*(Timesam); % FFT from all the rows
%% Fourier Transform Visualization

figure();
hold on
subplot(1,2,1);
plotcreator(real(fftdata1),colormap,'Horizontal Distance [m]','Frequency [Hz]','Fast Fourier Transform per column',Distance);
%Imagesc does not use imaginary values so we use the real part of the values.

subplot(1,2,2);
plotcreator(real(fftdata2),colormap,'Horizontal Distance [m]','Frequency [Hz]','Fast Fourier Transform per column',Distance); 
%Imagesc does not use imaginary values so we use the real part of the values.
hold off

fprintf(['Note that we get the following error ''complex values are not supported. Specify the color' ...
    '\ndata as numeric or logical values.''So to solve this problem we only use the real part of the data.\n'])

%% Filtering Data

fftdata1f=fftdata1(1:200,:);

%% Fourier Transform Visualization

figure();
hold on
subplot(1,2,1);
plotcreator(abs(fftdata1),colormap,'Horizontal Distance [m]','Frequency [Hz]','Fourier Transformed Data',Distance); 
%Imagesc does not use imaginary values so we use the real part of the values.

subplot(1,2,2);
plotcreator(abs(fftdata1f),colormap,'Horizontal Distance [m]','Frequency [Hz]','Filtered Fourier Transformed Data',Distance);  
%Imagesc does not use imaginary values so we use the real part of the values.
hold off

%% Direct Surface Wave

%{
Comparing the commonsource gather in the time-space and frequency-space domain, 
to what events in the time-space domain do you think does that energy correspond?
%}
fprintf(['Comparing the common source gather in the time-space and frequency-space domain,\n'... 
        'to what events in the time-space domain do you think does that energy correspond?\n\n'])

fprintf(['The energy we see is source generated waves in this assignment we treat it as noise\n' ...
    'since we are only interested in reflected waves and not direct surface waves. This signal\n' ...
    'is the wave starting at the source going straight to the receivers.\n\n'])

%% Low Cut Filter
fftdata1f_lcf=fftdata1f;
fftdata1f_lcf(1:lcf1max,1:lcf1d)= 0;

%Inverse Fourier Transform:
ifftdata1f_lcf=fftshift(fft(fftdata1f_lcf,[],1)*(Recevspace),1);

figure();
hold on
subplot(1,3,1);
plotcreator(abs(fftdata1f),colormap,'Horizontal Distance [m]','Frequency [Hz]','Fourier Transformed Data',Distance); 

subplot(1,3,2);
plotcreator(abs(fftdata1f_lcf),colormap,'Horizontal Distance [m]','Frequency [Hz]','Low cut filtered data',Distance)

subplot(1,3,3);
plotcreator(abs(ifftdata1f_lcf),colormap,'Horizontal Distance [m]','Time [s]','Time-Sapce filtered data',Distance)

hold off

%% Low Cut Linearity Filter

fftdata1f_lcfl=fftdata1f;
linearity=linspace(0,1,31);
fftdata1f_lcfl(lcf1min:lcf1max,1:lcf1d)= transpose(linearity).*fftdata1f(lcf1min:lcf1max,1:lcf1d);
fftdata1f_lcfl(1:lcf1min,1:lcf1d)= 0;

%Inverse Fourier Transform:
ifftdata1f_lcfl=fftshift(fft(fftdata1f_lcfl,[],1)*(Recevspace),1);

figure();
hold on
subplot(1,3,1);
plotcreator(abs(fftdata1f),colormap,'Horizontal Distance [m]','Frequency [Hz]','Fourier Transformed Data',Distance); 
%Imagesc does not use imaginary values so we use the real part of the values.

subplot(1,3,2);
plotcreator(abs(fftdata1f_lcfl),colormap,'Horizontal Distance [m]','Frequency [Hz]','Linear Low Cut Filtered Data',Distance);  
%Imagesc does not use imaginary values so we use the real part of the values.

subplot(1,3,3);
plotcreator(abs(ifftdata1f_lcfl),colormap,'Horizontal Distance [m]','Two-way Travel Time [s]','TIme Space Linear filtered',Distance);  

hold off

%% Low Cut Sinusoidal Filter

fftdata1f_lcfs=fftdata1f;

%Create a vector using linspace that will filter the data in a sinusoidal
%manner
sinusoidal=linspace(0,(sin(pi()/4)),((lcf1max-lcf1min)+1));
fftdata1f_lcfs(lcf1min:lcf1max,1:lcf1d)= transpose(sinusoidal).*fftdata1f(lcf1min:lcf1max,1:lcf1d);
fftdata1f_lcfs(1:lcf1min,1:lcf1d)= 0;

%Inverse Fourier Transform:
ifftdata1f_lcfs=fftshift(fft(fftdata1f_lcfs,[],1)*(Recevspace),1);

figure();
hold on
subplot(1,3,1);
plotcreator(abs(fftdata1f),colormap,'Horizontal Distance [m]','Frequency [Hz]','Fourier Transformed Data',Distance); 
%Imagesc does not use imaginary values so we use the real part of the values.

subplot(1,3,2);
plotcreator(abs(fftdata1f_lcfs),colormap,'Horizontal Distance [m]','Frequency [Hz]','Sinusoidal Low Cut Filtered Data',Distance);  
%Imagesc does not use imaginary values so we use the real part of the values.

subplot(1,3,3);
plotcreator(abs(ifftdata1f_lcfs),colormap,'Horizontal Distance [m]','Two-way Travel-time [s]','Time-Space Sinusoidal Filteried',Distance);
%Imagesc does not use imaginary values so we use the real part of the values.
hold off


%% Extra Filter images

%Filtered Data in Frequency-Space Domain
figure();
hold on
subplot(1,4,1);
plotcreator(abs(fftdata1f),colormap,'Horizontal Distance [m]','Frequency [Hz]','Fourier Transformed Data',Distance); 
%Imagesc does not use imaginary values so we use the real part of the values.

subplot(1,4,4);
plotcreator(abs(fftdata1f_lcfs),colormap,'Horizontal Distance [m]','Frequency [Hz]','Sinusoidal Low cut filtered data',Distance);  
%Imagesc does not use imaginary values so we use the real part of the values.

subplot(1,4,3);
plotcreator(abs(fftdata1f_lcfl),colormap,'Horizontal Distance [m]','Frequency [Hz]','Linear Low cut filtered data',Distance);  
%Imagesc does not use imaginary values so we use the real part of the values.

subplot(1,4,2);
plotcreator(abs(fftdata1f_lcf),colormap,'Horizontal Distance [m]','Frequency [Hz]','Low cut filtered data',Distance);
%Imagesc does not use imaginary values so we use the real part of the values.
hold off

%Filtered Data in Time-Space Domain
n=figure();
hold on
subplot(1,4,1);
plotcreator(abs(fftdata1),gray,'Horizontal Distance [m]','Time [s]','Fourier Transformed Data',Distance); 
%Imagesc does not use imaginary values so we use the real part of the values.

subplot(1,4,4);
plotcreator(abs(ifftdata1f_lcfs),gray,'Horizontal Distance [m]','Time [s]','Sinusoidal Filter',Distance);  
%Imagesc does not use imaginary values so we use the real part of the values.

subplot(1,4,3);
plotcreator(abs(ifftdata1f_lcfl),gray,'Horizontal Distance [m]','Time [s]','Linear Filter',Distance);  
%Imagesc does not use imaginary values so we use the real part of the values.

subplot(1,4,2);
plotcreator(abs(ifftdata1f_lcf),gray,'Horizontal Distance [m]','Time [s]','Filter',Distance);
%Imagesc does not use imaginary values so we use the real part of the values.
hold off

fprintf(['What can be observed about the events – what events are preserved and what is their preservation \n' ...
    'quality? Which filter gave the clearest result?\n\n'])

fprintf(['As we can see from figure %g the linear filtered gave us clearer results than the other two filters\n' ...
    'applied to our data. Giving us sharper results and visualizating the reflected waves by removind the\n' ...
    'direct surface waves.\n\n'],n)

%% Frequency-Space Domain to Frequency-Wavenumber Domain

freqtowave_lcf=fftshift(fft(fftdata1,[],1)*(Recevspace),1);

figure()

hold on
subplot(1,1,1);
plotcreator(abs(freqtowave_lcf),colormap,'Wavenumber [-]','Frequency [Hz]','Frequency to Wavenumber Domain',wavenumber);
hold off

fprintf(['What can be observed about the energy?\n'])
