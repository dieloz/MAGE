function ...
    g=complete_chirplet_parameters(g)

% g=complete_chirplet_parameters(g)
%
% preferred usage:
% g.center_time=0;
% g.center_frequency=40;
% g.duration_parameter=-4;
% g.chirp_rate=0;
% g=complete_chirplet_parameters(g);

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


% must have these:
if isfield(g,'center_time'); t0=g.center_time; else t0=0; end
if isfield(g,'center_frequency'); v0=g.center_frequency; end
if isfield(g,'chirp_rate'); c0=g.chirp_rate; end

% assign or compute duration parameter:
if isfield(g,'duration_parameter')
    s0=g.duration_parameter;
elseif isfield(g,'time_domain_standard_deviation')
    tstd=g.time_domain_standard_deviation;
    s0=log(4*pi*tstd^2);
elseif isfield(g,'frequency_domain_standard_deviation')&&(c0==0)
    fstd=g.frequency_domain_standard_deviation;
    s0=-log(4*pi*fstd^2);    
elseif isfield(g,'fractional_bandwidth')&&(c0==0)
    fbw=g.fractional_bandwidth;
    s0=log((2*log(2))/(fbw^2*pi*v0^2));
else
    error('If g.chirp_rate~=0, then must specify either g.duration_parameter or g.time_domain_standard_deviation; g.frequency_domain_standard_deviation and g.fractional_bandwidth not permitted');
end

% start with new version:
clear g
g.center_time=[];
g.center_frequency=[];
g.duration_parameter=[];
g.chirp_rate=[];
g.fractional_bandwidth=[];
g.time_domain_standard_deviation=[];
g.frequency_domain_standard_deviation=[];
g.time_frequency_covariance=[];
g.time_frequency_correlation_coefficient=[];
g.std_multiple_for_support=[];

% reenter given parameters:
g.center_time=t0;
g.center_frequency=v0;
g.duration_parameter=s0;
g.chirp_rate=c0;
% fixed parameter:
g.std_multiple_for_support=6;
% calculate other properties:
g.time_domain_standard_deviation=...
    sqrt(exp(s0)/(4*pi));
g.frequency_domain_standard_deviation=...
    sqrt((exp(-s0)+c0^2*exp(s0))/(4*pi));
g.time_frequency_covariance=(c0*exp(s0))/(4*pi);
g.time_frequency_correlation_coefficient=...
    g.time_frequency_covariance*(...
    g.time_domain_standard_deviation*...
    g.frequency_domain_standard_deviation)^-1;
g.fractional_bandwidth=2*sqrt(2*log(2))*...
    g.frequency_domain_standard_deviation/v0;% fwhm/center_frequency
