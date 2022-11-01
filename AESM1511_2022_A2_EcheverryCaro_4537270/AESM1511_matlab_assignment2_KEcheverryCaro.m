%% AESM1511 MATLAB ASSIGNMENT 2:Fourier Transformation and Filtering
%{
Title: Assignment 2
Course Code: AESM1511
Author: Karla Echeverry Caro
Studentnumber: 4537270
Date Created: 04 October 2022
Date modified: 09 October 2022
Mail: K.V.EcheverryCaro@student.tudelft.nl
Assignment 2:  Fourier Transformation and Filtering

Status: COMPLETED

Other m-files required:
function    plorcreator.m
function    plorcreatorb.m
MAT-files required: none

Parameter Dictionary:
path       [-]          Filepath to seismic measurement data
name       [-]          Name of the seismic measurement data
no_traces  [-]          Number of sampling times
SR         [-]          Location of the starting receiver
LR         [-]          Location of the last receiver
xdt        [-]          Traces Vector 
dtt        [minisec]    Total sampled time
dt         [minisec]    Sampling time
ds         [meters]     Sampling space of the receivers
xd         [meters]     Horizontal displacement array
minfreq    [Hz]         Minimal frequency used in the frequency domain
maxfreq    [Hz]         Maximal frequency used in the frequency domain
minffreq   [Hz]         Minimal frequency used to filter the data in the
                        frequency domain caused by Direct Surface Waves in
                        this assignment
maxffreq   [Hz]         Maximal frequency used to filter the data in the
                        frequency domain caused by Direct Surface Waves in
                        this assignment
maxffdis      [meters]     Minimal displacement used to filter the data in the
                        frequency domain
lcf1max    [meters]     Maximal displacement used to filter the data in the
                        frequency domain
tf_s       [minisec]    Time frequency sampling vector
sf_s       [Hz]         Space frequency sampling
linearity  [-]          Linear Filter for the Frequency domain
sinusoidal [-]          Sinusoidal Filter for the Frequency domain
%}
%% Initial Settings
clc;clear;close all;
%% Path and name of the file
path=['C:\Users\khato\OneDrive\Desktop\Signal & Analysis\Matlab\AESM1511_2022' ...
    '_A2_EcheverryCaro_4537270'];%path
name='refl_3layers_fp50_dx0p5_500_rvz.bin';% name of the file
parameter=append(path,'\',name);%combine path and name
%% Initial Input
%Time Information
dt=1e-3; 
dtt=1001; 

%Traces Information
ds=2.5; 
SR= 510; 
LR= 1510; 
no_traces=401;

%Frequency Domain Information
minfreq=1;
maxfreq=200;

%Filter Parameters
maxffreq=60; 
minffreq=30; 
minffdis=0; 
maxffdis=185;

%Densities are ignored for this assignment!
%% Parameters
xd=linspace(0,LR-SR,no_traces);
xdt=1:no_traces;
freq=minfreq:maxfreq;
tf_s = 1/(dtt*dt);
sf_s = 1/(no_traces*ds);
%% Read Bin File
fid=fopen(parameter,'r');
temp=fread(fid,no_traces*dtt,'float32',0,'ieee-le');
fclose(fid);
data_refl(:,:)=reshape(temp,dtt,no_traces);
%% Background Information
%{
Ricker Wavelet:
Negative normalized second derivative of a Gaussian gunction up to
scale and normalisation the second hermite function.
%}
%% Transformation From Time Domain To Frequency Domain
% This is done using Fast Fourier Transform 
% fft(data in time-space domain,[],?)*(time sampling)
fftdata1=fft(data_refl,[],1)*(dt); % FFT from all the columns
%% Cropping The Data
fftdata1f=fftdata1(minfreq:maxfreq,:);
%% Low Cut Filter
fftdata1f_lcf=fftdata1f;
fftdata1f_lcf(1:maxffreq,1:maxffdis)= 0;
%% Linear Low Cut Filter
fftdata1f_lcfl=fftdata1f;
linearity=linspace(0,1,31);
fftdata1f_lcfl(minffreq:maxffreq,1:maxffdis)= transpose(linearity).*fftdata1f(minffreq:maxffreq,1:maxffdis);
fftdata1f_lcfl(1:minffreq,1:maxffdis)= 0;
%% Sinusoidal Low Cut Filter
fftdata1f_lcfs=fftdata1f;
sinusoidal=linspace(0,(sin(pi()/2)),((maxffreq-minffreq)+1));
fftdata1f_lcfs(minffreq:maxffreq,1:maxffdis)= transpose(sinusoidal).*fftdata1f(minffreq:maxffreq,1:maxffdis);
fftdata1f_lcfs(1:minffreq,1:maxffdis)= 0;
%% Inverse Fourier Transform: Transform the data from Frequency-Space Domain to Time-Space Domain
%Linearity
ifftdata1f_lcfl=real(ifft(fftdata1f_lcfl,[],1))*(maxfreq)*2*tf_s;
%Sinusoidal
ifftdata1f_lcfs=real(ifft(fftdata1f_lcfs,[],1))*(maxfreq)*2*tf_s;
%Low Cut Filter
ifftdata1f_lcf=real(ifft(fftdata1f_lcf,[],1))*(maxfreq)*2*tf_s;
%% Frequency-Space Domain to Frequency-Wavenumber Domain
%Original Data
original_data=fftshift(fft(fftdata1f(:,1:1:end),[],2),2)*ds;
%Filter
freqtowave_lcf=fftshift(fft(fftdata1f_lcf(:,1:1:end),[],2),2)*ds;
%Linear Filter
freqtowave_lcfl=fftshift(fft(fftdata1f_lcfl(:,1:1:end),[],2),2)*ds;
%Sinusoidal Filter
freqtowave_lcfs=fftshift(fft(fftdata1f_lcfs(:,1:1:end),[],2),2)*ds;
%% Frequency-Wavenumber Domain every knd Trace
%Sinusoidal Filter
%k=2
data_wave_2 = fftshift(fft(fftdata1f(:,1:2:end),[],2),2)*ds;
%k=4
data_wave_4 = fftshift(fft(fftdata1f(:,1:4:end),[],2),2)*ds;
%k=8
data_wave_8 = fftshift(fft(fftdata1f(:,1:8:end),[],2),2)*ds;
%% Bin File Visualization
figure();
plotcreator(data_refl,'gray','Horizontal Displacement [m]','Two-Way Travel Time  [ms]','Common Source-Gather',xd)
%% Fourier Transform Visualization
figure();
plotcreator(abs(fftdata1f),'colormap','Horizontal xd [m]','Frequency [Hz]','Filtered Fourier Transformed Data',xd); 
%% Visualization: Filtered Data In Frequency-Space Domain
figure();
hold on
subplot(2,2,1);
plotcreatorb(abs(fftdata1f),'colormap','Horizontal Displacement [m]','Frequency [Hz]','Fourier Transformed Data',xd,[freq(1), freq(2)]); 
subplot(2,2,2);
plotcreatorb(abs(fftdata1f_lcf),'colormap','Horizontal Displacement [m]','Frequency [Hz]','Low cut filtered data',xd,[freq(1), freq(2)]);
subplot(2,2,3);
plotcreatorb(abs(fftdata1f_lcfl),'colormap','Horizontal Displacement [m]','Frequency [Hz]','Linear Low cut filtered data',xd,[freq(1), freq(2)]);  
subplot(2,2,4);
plotcreatorb(abs(fftdata1f_lcfs),'colormap','Horizontal Displacement [m]','Frequency [Hz]','Sinusoidal Low cut filtered data',xd,[freq(1), freq(2)]);  
hold off
%% Visualization: Filtered Data in Time-Space Domain
n=figure();
hold on
subplot(2,2,1);
plotcreator(data_refl,'gray','Horizontal Displacement [m]','Two-Way Travel Time [ms]','Original Data',xd); 
subplot(2,2,2);
plotcreator(ifftdata1f_lcf,'gray','Horizontal Displacement [m]','Two-Way Travel Time [ms]','Filter',xd);
subplot(2,2,3);
plotcreator(ifftdata1f_lcfl,'gray','Horizontal Displacement [m]','Two-Way Travel Time [ms]','Linear Filter',xd);  
subplot(2,2,4);
plotcreator(ifftdata1f_lcfs,'gray','Horizontal Displacement [m]','Two-Way Travel Time [ms]','Sinusoidal Filter',xd);  
hold off
%% Frequency-Wavenumber Domain Visualization
s1=figure();
hold on
subplot(2,2,1)
plotcreatorb(abs(original_data),'colormap','Wavenumber [1/m]','Frequency [Hz]',{'Traces In Frequency-Wavenumber Domain ','Original Data'}, [-0.5*sf_s*no_traces, 0.5*sf_s*no_traces],[freq(1), freq(2)]);
subplot(2,2,2)
plotcreatorb(abs(freqtowave_lcf),'colormap','Wavenumber [1/m]','Frequency [Hz]',{'Traces In Frequency-Wavenumber Domain ',''}, [-0.5*sf_s*no_traces, 0.5*sf_s*no_traces],[freq(1), freq(2)]);
subplot(2,2,3)
plotcreatorb(abs(freqtowave_lcfl),'colormap','Wavenumber [1/m]','Frequency [Hz]',{'Traces In Frequency-Wavenumber Domain ','Linear'}, [-0.5*sf_s*no_traces, 0.5*sf_s*no_traces],[freq(1), freq(2)]);
subplot(2,2,4)
plotcreatorb(abs(freqtowave_lcfs),'colormap','Wavenumber [1/m]','Frequency [Hz]',{'Traces In Frequency-Wavenumber Domain ','Sinusoidal'}, [-0.5*sf_s*no_traces, 0.5*sf_s*no_traces],[freq(1), freq(2)]);
hold off
%% Frequency-Wavenumber Domain every knd Trace Visualization 
s2=figure();
hold on
subplot(2,2,1)
plotcreatorb(abs(freqtowave_lcfs),'colormap','Wavenumber [1/m]','Frequency [Hz]',{'Traces In Frequency-Wavenumber Domain ','Every Trace'}, [-0.5*sf_s*no_traces, 0.5*sf_s*no_traces],[freq(1), freq(2)]);
subplot(2,2,2)
plotcreatorb(abs(data_wave_2),'colormap','Wavenumber [1/m]','Frequency [Hz]',{'Traces In Frequency-Wavenumber Domain ','Every 2nd Trace'}, [-0.5*sf_s*no_traces, 0.5*sf_s*no_traces],[freq(1), freq(2)]);
subplot(2,2,3)
plotcreatorb(abs(data_wave_4),'colormap','Wavenumber [1/m]','Frequency [Hz]',{'Traces In Frequency-Wavenumber Domain ','Every 4nd Trace'}, [-0.5*sf_s*no_traces, 0.5*sf_s*no_traces],[freq(1), freq(2)]);
subplot(2,2,4)
plotcreatorb(abs(data_wave_8),'colormap','Wavenumber [1/m]','Frequency [Hz]',{'Traces In Frequency-Wavenumber Domain ','Every 8nd Trace'}, [-0.5*sf_s*no_traces, 0.5*sf_s*no_traces],[freq(1), freq(2)]);
hold off
%% Answering Question
fprintf(['Question 2: Note that we get the following error ''complex values are not supported. Specify the color' ...
    '\ndata as numeric or logical values.''So to solve this problem we only use the real part of the data.\n\n'])
fprintf(['Question 3: Comparing the common source gather in the time-space and frequency-space domain,\n'... 
        'to what event in the time-space domain do you think does that energy correspond?\n\n'])
fprintf(['The energy we see is source generated waves in this assignment we treat it as noise\n' ...
    'since we are only interested in reflected waves and not direct surface waves. This signal\n' ...
    'is the wave starting at the source going straight to the receivers.\n\n'])
fprintf(['Question 4: What can be observed about the event – what evendtt are preserved and what is their preservation \n' ...
    'quality? Which filter gave the clearest result?\n\n'])
fprintf(['As we can see from figure %g the linear filtered gave us clearer result than the other two filters\n' ...
    'applied to our data. Giving us sharper resuldtt and visualizating the reflected waves by removind the\n' ...
    'direct surface waves.\n\n'],n)
fprintf('Question 5: What can be observed about the energy?\n\n')
fprintf('Here we can see that the energy is filtered from the data as we can see in figures %g.\n\n',s1)
fprintf(['Question 6: Repeat Task 5, but now by taking every 2nd, 4th, and 8th trace. Visualize the ' ...
    'three images and the\nimage from Task 5 in one figure using subplot (using the'...
'same parameters as in Task 2). Comparing the\nfour images, explain what happens to the main energy.\n\n'])
fprintf(['The energy is now coming from other areas of the plot this is caused by alisasing, with other words\n' ...
    'different signals become indistinguishable when sampled as can be seen in figure %g.\n\n'],s2)
fprintf(['PD: So I didnt do it for every image created before so now I only used the sinusoidal filter to show the\n' ...
    'result.\n\n'])
%-------------------------------------------------END------------------------------------------------------



