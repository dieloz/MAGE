function ...
    filter_output=mage_filter(mage_decomp_output,min_f,max_f)

% filter_output=mage_filter(mage_decomp_output,min_f,max_f)

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


 p=get_signal_parameters(...
      'sampling_rate',mage_decomp_output.srate,...
      'number_points_time_domain',length(mage_decomp_output.time_support));

clear filter_output
filter_output.gparams_included=[];
filter_output.gparams_excluded=[];
filter_output.raw_signal=mage_decomp_output.raw_signal;
filter_output.filtered_signal=[]; % from gparams_included
filter_output.excluded_signal=[]; % from gparams_excluded
filter_output.residual=mage_decomp_output.residual;
filter_output.time_support=mage_decomp_output.time_support;
filter_output.srate=mage_decomp_output.srate;

gps=mage_decomp_output.gparams;
drop_inds1=find(gps(:,2)<min_f);
drop_inds2=find(gps(:,2)>max_f);
drop_inds=unique([drop_inds1(:);drop_inds2(:)]);
drop_inds=drop_inds(:);
keep_inds=setdiff(1:size(gps,1),drop_inds);

filter_output.gparams_included=gps(keep_inds,:);
filter_output.gparams_excluded=gps(drop_inds,:);

filter_output.filtered_signal=mage_params2signal(filter_output.gparams_included,p);
filter_output.excluded_signal=mage_params2signal(filter_output.gparams_excluded,p);



