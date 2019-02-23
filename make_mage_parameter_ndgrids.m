function ...
    mage_parameters=make_mage_parameter_ndgrids(...
    min_freq,max_freq,num_freqs,min_freq_step,...
    min_fbw,max_fbw,num_fbw,...
    num_tsamples,srate)

%     mage_parameters=make_mage_parameter_ndgrids(...
%     min_freq,max_freq,num_freqs,min_freq_step,...
%     min_fbw,max_fbw,num_fbw,...
%     num_tsamples,srate)

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


% make decomp parameter space ndgrids tt_array, vt_array, st_array
cf_list=make_center_frequencies(min_freq,max_freq,num_freqs,min_freq_step);
fbw_list=linspace(min_fbw,max_fbw,num_fbw);
p=get_signal_parameters(...
      'sampling_rate',srate,...
      'number_points_time_domain',num_tsamples);
 
[tt_array,vt_array,fbw_array]=ndgrid(p.time_support,cf_list,fbw_list);
st_array=vpfbw2sp(vt_array,fbw_array);

clear mage_parameters
mage_parameters.center_time_support=p.time_support;
mage_parameters.center_frequency_support=cf_list;
mage_parameters.fbw_support=fbw_list;
mage_parameters.tt_array=single(tt_array);
mage_parameters.vt_array=single(vt_array);
mage_parameters.st_array=single(st_array);
mage_parameters.srate=srate;
