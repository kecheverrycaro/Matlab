%% ASSIGNMENT2 - Programming Assignment 2 for AESM1511, WS22, TU Delft
% Program to visualize seismic measurement data in different domains,
% convert data between domains and filter data with different filters.
%
% Other m-files required:
%   function    apply_filter.m
%   function    get_plot.m
% MAT-files required: none
%
% Author: Malte Leander Schade
% Student ID: 5850282
% Mail: contact@malteschade.com
% Copyright: Copyright 2022, Assignment2-AESM1511-Malte-Schade'
% Version: 1.2.0
% Date created: 07.10.22
% Date modified: 07.10.22
% Status: COMPLETED
%
% Parameter Dictionary:
%   fpath       [-]     filepath to seismic measurement data
%   n_traces    [-]     number of measurement traces
%   n_samples   [-]     number of sampling times
%   t_s         [s]     sampling time
%   s_s         [m]     sampling space
%   f_bp        [Hz]    bandpass filter frequency intervall
%   f_soft      [Hz]    soft filter frequency
%   f_hard      [Hz]    hard filter frequency
%   s_lim       [m]     filter distance limit


%% SETTINGS
clear all;  % clear all previous variables
close all;  % close all figures
clc;        % clear commandline


%% PARAMETERS
fpath = 'refl_3layers_fp50_dx0p5_500_rvz.bin';
n_traces = 401;
n_samples = 1001;
t_s = 0.001;
s_s = 2.5;
f_bp = [1,200];
f_soft = 30;
f_hard = 60;
s_lim = 450;


%% IMPORTS
% import dataset from binary file and restore dimensionality
fid=fopen(fpath,'r');
temp=fread(fid,n_traces*n_samples,'float32',0,'ieee-le');
fclose(fid);
data_refl(:,:)=reshape(temp,n_samples,n_traces);


%% CALCULATIONS
% calculate time frequency sampling
tf_s = 1/(n_samples*t_s);

% calculate space frequency sampling
sf_s = 1/(n_traces*s_s);

% convert filter frequencies to data position
fp_bp = ceil(tf_s*f_bp);
fp_soft = ceil(tf_s*f_soft);
fp_hard = ceil(tf_s*f_hard);
sp_lim = ceil(s_lim/s_s);

% transform to frequency-space domain and multiply with sampling time
data_freq = fft(data_refl,[], 1)*t_s;

% filter frequency band with bandpass filter
data_freq = data_freq(fp_bp(1):fp_bp(2),:);

% filter data by soft and hard filter frequency and space limit with udf
data_filter_1 = apply_filter(data_freq, 'bri', fp_soft, fp_hard, sp_lim);
data_filter_2 = apply_filter(data_freq, 'lin', fp_soft, fp_hard, sp_lim);
data_filter_3 = apply_filter(data_freq, 'sin', fp_soft, fp_hard, sp_lim);

% transform to time-space domain and multiply with number of frequencies,
% time frequency sampling and 2
data_time_1 = real(ifft(data_filter_1,[],1))*(1+f_bp(2)-f_bp(1))*2*tf_s;
data_time_2 = real(ifft(data_filter_2,[],1))*(1+f_bp(2)-f_bp(1))*2*tf_s;
data_time_3 = real(ifft(data_filter_3,[],1))*(1+f_bp(2)-f_bp(1))*2*tf_s;

% transform to frequency-wavenumber domain, shift and multiply with
% sampling space
data_wave_1 = fftshift(fft(data_freq(:,1:1:end),[],2),2)*s_s;
data_wave_2 = fftshift(fft(data_freq(:,1:2:end),[],2),2)*s_s;
data_wave_4 = fftshift(fft(data_freq(:,1:4:end),[],2),2)*s_s;
data_wave_8 = fftshift(fft(data_freq(:,1:8:end),[],2),2)*s_s;


%% VISUALIZATIONS
% Plot 1: Response dataset in frequency-space-domain
subplot(2,2,1);
get_plot(abs(data_freq), [0, s_s*n_traces], [f_bp(1), f_bp(2)],...
    'default', 'Distance [m]', 'Frequency [Hz]',...
    'Unfiltered dataset in frequency-space-domain')

% Plot 2: Filter 2 dataset in frequency-space-domain
subplot(2,2,2);
get_plot(abs(data_filter_1), [0, s_s*n_traces], [f_bp(1), f_bp(2)],...
    'default', 'Distance [m]', 'Frequency [Hz]',...
    'Hard filter dataset in frequency-space-domain')

% Plot 3: Filter 3 dataset in frequency-space-domain
subplot(2,2,3);
get_plot(abs(data_filter_2), [0, s_s*n_traces], [f_bp(1), f_bp(2)],...
    'default', 'Distance [m]', 'Frequency [Hz]',...
    'Linear filter dataset in frequency-space-domain')

% Plot 4: Filter 4 dataset in frequency-space-domain
subplot(2,2,4);
get_plot(abs(data_filter_3), [0, s_s*n_traces], [f_bp(1), f_bp(2)],...
    'default', 'Distance [m]', 'Frequency [Hz]',...
    'Sinusoidal filter dataset in frequency-space-domain')

figure;

% Plot 5: Response dataset in time-space-domain
subplot(2,2,1);
get_plot(real(data_refl), [0, s_s*n_traces], [0, t_s*n_samples],...
    'gray', 'Distance [m]', 'Time [s]',...
    'Unfiltered dataset in time-space-domain')

% Plot 6: Filter 1 dataset in time-space-domain
subplot(2,2,2);
get_plot(real(data_time_1), [0, s_s*n_traces], [0, t_s*n_samples],...
    'gray', 'Distance [m]', 'Time [s]',...
    'Hard filter dataset in time-space-domain')

% Plot 7: Filter 2 dataset in time-space-domain
subplot(2,2,3);
get_plot(real(data_time_2), [0, s_s*n_traces], [0, t_s*n_samples],...
    'gray', 'Distance [m]', 'Time [s]',...
    'Linear filter dataset in time-space-domain')

% Plot 8: Filter 3 dataset in time-space-domain
subplot(2,2,4);
get_plot(real(data_time_3), [0, s_s*n_traces], [0, t_s*n_samples],...
    'gray', 'Distance [m]', 'Time [s]',...
    'Sinusoidal filter dataset in time-space-domain')

figure;

% Plot 9: Every Response dataset in frequency-wavenumber-domain
subplot(2,2,1);
get_plot(abs(data_wave_1), [-0.5*sf_s*n_traces, 0.5*sf_s*n_traces], ...
    [f_bp(1), f_bp(2)], ...
    'default', 'Wavenumber [1/m]', 'Frequency [Hz]',...
    'Every trace in frequency-wavenumber-domain')

% Plot 10: Every 2nd Response dataset in frequency-wavenumber-domain
subplot(2,2,2);
get_plot(abs(data_wave_2), [-0.5*sf_s*n_traces, 0.5*sf_s*n_traces], ...
    [f_bp(1), f_bp(2)], ...
    'default', 'Wavenumber [1/m]', 'Frequency [Hz]',...
    'Every 2nd trace in frequency-wavenumber-domain')

% Plot 11: Every 4th Response dataset in frequency-wavenumber-domain
subplot(2,2,3);
get_plot(abs(data_wave_4), [-0.5*sf_s*n_traces, 0.5*sf_s*n_traces], ...
    [f_bp(1), f_bp(2)], ...
    'default', 'Wavenumber [1/m]', 'Frequency [Hz]',...
    'Every 4th trace in frequency-wavenumber-domain')

% Plot 12: Response dataset in frequency-wavenumber-domain
subplot(2,2,4);
get_plot(abs(data_wave_8), [-0.5*sf_s*n_traces, 0.5*sf_s*n_traces], ...
    [f_bp(1), f_bp(2)], ...
    'default', 'Wavenumber [1/m]', 'Frequency [Hz]',...
    'Every 8th trace in frequency-wavenumber-domain')


%% DISPLAYS
disp('Question 1: What happens when attempting to plot ')
disp('and why that happens?');
disp('It is impossible to plot the result as it is a complex number.');
disp('It is necessary to get the absolute value.')
disp(' ');
disp('Question 2: Should this be done for all horizontal distances?');
disp('No a limit to the direct wave x measurement range is sufficient');
disp('and saves relevant data.');
disp(' ')
disp('Question 3: To what events in');
disp('the time-space domain does that energy correspond?');
disp('The offset is placed due to the end of the experiment time and that')
disp('the direct wave has not reached further receivers at that time.');
disp(' ');
disp('Question 4: What events are preserved and what');
disp('is their preservation quality?');
disp('Low freq. direct waves at the experiment beginning are reduced.');
disp('Brickwall filtering at 60hz reduces quality substantially.');
disp('The Sinusoidal filter gives the best results.');
disp(' ');
disp('Question 5: What happens to the main energy?');
disp('The main energy drags around the plot and comes in from');
disp('the other side. This is an effect of alisasing.');