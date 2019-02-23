function ...
    [filtered_signal,residual,gdict,gparams]=linear_filter_burst_finder(...
    raw,srate,center_frequency,num_bursts)

%   [filtered_signal,residual,gdict,gparams]=linear_filter_burst_finder(...
%     raw,srate,center_frequency,num_bursts)

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


p=get_signal_parameters('sampling_rate',srate,... % Hz
    'number_points_time_domain',length(raw));

% filter with initial probe gabor
clear g
g.center_frequency=center_frequency;
g.fractional_bandwidth=0.2;
g.chirp_rate=0;
g=make_chirplet(...
    'chirplet_structure',g,...
    'signal_parameters',p);
fs=filter_with_chirplet(...
    'raw_signal',raw,...
    'signal_parameters',p,...
    'chirplet',g);

% only look at amp local maxima for linear filter
[maxvals,maxinds]=find_maxima(abs(fs.time_domain));
% 
% filt_signal=eegfilt(raw,srate,min_freq,max_freq);
% s=make_signal_structure(...
%       'raw_signal',filt_signal,...
%       'output_type','analytic',...
%       'signal_parameters',p);
% f2=s.time_domain;  
% [maxvals,maxinds]=find_maxima(abs(fs2));

[maxvals,maxinds]=find_maxima(abs(fs.time_domain));
tt_est=p.time_support(maxinds);
vt_est=repmat(g.center_frequency,1,length(maxinds));
st_est=repmat(g.duration_parameter,1,length(maxinds));
gparams=[tt_est; vt_est; st_est;...
    abs(fs.time_domain(maxinds));...
    angle(fs.time_domain(maxinds))];

gparams=gparams(:,1:num_bursts);

gdict=mage_params2dict(gparams,p);

filtered_signal=mage_dict2signal(gdict,p);

residual=raw-filtered_signal;


%     [filtered_signal,residual,gdict,gparams]=linear_filter_burst_finder(...
%     raw,srate,min_freq,max_freq,num_bursts)




