% Example of Multiscale Adaptive Gabor Expansion (MAGE)
%
% 1. Create simulated target signal with known ground truth, adds small noise
% 2. Extracts Natom gabors from the target signal via MAGE
% 3. Filter with Gaussian-envelope Gabor atom (FIR)

% Author: Ryan T. Canolty (2019) rcanolty@gmail.com http://accl.psy.vanderbilt.edu/
% For additional information, visit https://www.biorxiv.org to read:
% TITLE: Multiscale Adaptive Gabor Expansion (MAGE):
% Improved Detection of Transient Oscillatory Burst Amplitude and Phase
% AUTHORS: Ryan T. Canolty and Thilo Womelsdorf  
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


% clear
% clc
% close all
format short g

sd=20190211;
rng(sd);

% 1. Create simulated target signal with known ground truth, adds small noise
% target signal sampling rate srate and time support t
srate=1000;
t=(1/srate)*(0:2^13-1);
p=get_signal_parameters(...
    'sampling_rate',srate,...
    'number_points_time_domain',length(t));

% Establish support for gabor parameters
tbounds=t([1 end]); % sec
vbounds=[1 50]; % Hz
sbounds=[-6 0]; % au
Abounds=[1 3]; % au
rpbounds=[-pi pi]; % radians

% Select parameters for Ntarget target gabors gt_act (actual)
% random selection from uniform distribution
% interval=maximum-minimum;
% samples=minimum+rand(number_of_samples,1)*interval;

Ntarget=100;
% center time tt_act
tt_act=tbounds(1)+rand(Ntarget,1)*(tbounds(2)-tbounds(1));
% center frequency vt_act
vt_act=vbounds(1)+rand(Ntarget,1)*(vbounds(2)-vbounds(1));
% duration parameter st_act
st_act=sbounds(1)+rand(Ntarget,1)*(sbounds(2)-sbounds(1));
% amplitude At_act
At_act=Abounds(1)+rand(Ntarget,1)*(Abounds(2)-Abounds(1));
% real-value phase rpt_act
rpt_act=rpbounds(1)+rand(Ntarget,1)*(rpbounds(2)-rpbounds(1));

% actual (ground-truth) target gabor parameters
% [center_time center_frequency duration_parameter amplitude(>0) phase[-pi,pi)]
gparams_act=[tt_act vt_act st_act At_act rpt_act];
% hand adjust act duration parameter for gabor # 3
gparams_act(3,3)=-1+gparams_act(3,3);

disp(gparams_act);

% synth signal (bursts only)
gsignal=mage_params2signal(gparams_act,p);

noise=0.1*randn(size(gsignal));

% signal to be given to MAGE
target_signal=real(gsignal)+noise;
raw_signal=target_signal;

% 2. extracts Natom gabors from the target signal via MAGE

Natoms=100; % number of time-freq-scale atoms to extract from signal

% input mage parameters as ndgrids (parameter support)
mage_parameters=make_mage_parameter_ndgrids(...
    2,200,64,0.5,... % min_freq,max_freq,num_freqs,min_freq_step,...
    0.1,0.6,8,... % min_fbw,max_fbw,num_fbw,...
    p.number_points_time_domain,p.sampling_rate);

% mage decomposition
mage_decomp_output=mage_decomp(raw_signal,Natoms,mage_parameters);

% mage filter 40-50 Hz
min_f=40; % Hz, low cut-off frequency
max_f=50; % Hz  high cut-off frequency
mage_filtered_output=mage_filter(mage_decomp_output,min_f,max_f);

fig=figure;
set(fig,'Color','w');
t=mage_filtered_output.time_support;
h=subplot(4,1,1);
y1=real(mage_filtered_output.raw_signal);
ybounds=max(abs(y1))*[-1 1];
plot(t,y1,'k');
ylim(ybounds);
title('Original Signal');
xlim(p.t([1 end]));
h=subplot(4,1,2);
y2=real(mage_filtered_output.filtered_signal);
plot(t,y2,'r');
ylim(ybounds);
title('Filtered Signal 40-50 Hz');
xlim(p.t([1 end]));
h=subplot(4,1,3);
y3=real(mage_filtered_output.excluded_signal);
plot(t,y3,'r');
ylim(ybounds);
title('Excluded Signal <40, >50 Hz');
xlim(p.t([1 end]));
h=subplot(4,1,4);
y4=real(mage_filtered_output.residual);
plot(t,y4,'b');
ylim(ybounds);
title('Residual');
xlabel('Time (sec)');
xlim(p.t([1 end]));

% 3. Filter with Gaussian-envelope Gabor atom (FIR)
vk=40; % Hz, center frequency of FIR filter
sk=-5; % unitless, duration parameter

clear g0
g0.center_time=0;
g0.center_frequency=vk;
g0.duration_parameter=sk;
g0.chirp_rate=0;
g0=complete_chirplet_parameters(g0);

fbw=g0.fractional_bandwidth; % unitless, bandwidth/center_frequency
f_sk=gabor_filter3_sk(target_signal',p.sampling_rate,vk,sk);
f_fbw=gabor_filter2_fbw(target_signal',p.sampling_rate,vk,fbw);

fig=figure;
set(fig,'Color','w');
h=subplot(2,1,1);
plot(p.t,target_signal,'k');
title('Target Signal = Bursts + Noise');
xlim(p.t([1 end]));
ylim([-1 1]);
h=subplot(2,1,2);
plot(p.t,real(f_sk),'r');
hold on;
plot(p.t,zeros(size(p.t)),':k');
hold off;
title('FIR Gaussian-envelope Gabor atom filtering example');
xlabel('Time (seconds)');
xlim(p.t([1 end]));


