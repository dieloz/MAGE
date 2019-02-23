function ...
    fs=filter_with_chirplet(varargin)

%     fs=filter_with_chirplet(...
%     'signal_structure',s,...
%     'signal_parameters',p,...
%     'chirplet',g);
% OR
%     fs=filter_with_chirplet(...
%     'raw_signal',raw_signal,...
%     'signal_parameters',p,...
%     'chirplet',g);

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


% test to see if the cell varargin was passed directly from
% another function; if so, it needs to be 'unwrapped' one layer
if length(varargin)==1 % should have at least 2 elements
    varargin=varargin{1};
end

for n=1:2:length(varargin)-1
    switch lower(varargin{n})
        case 'signal_structure'
            s=varargin{n+1};
        case 'signal_parameters'
            p=varargin{n+1};
        case 'chirplet'
            g=varargin{n+1};
    end
end

% in case raw, time-domain signal is input rather than signal structure
for n=1:2:length(varargin)-1
    switch lower(varargin{n})
        case 'raw_signal'
            raw_signal=varargin{n+1};
            if size(raw_signal,1)>size(raw_signal,2)
                raw_signal=raw_signal';
            end
    end
end
if exist('raw_signal','var')
    s=make_signal_structure(...
    'raw_signal',raw_signal,...
    'output_type','analytic',...
    'signal_parameters',p);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
fs.time_domain=zeros(size(s.time_domain));
fs.frequency_domain=zeros(size(s.frequency_domain));

fs.frequency_domain(g.signal_frequency_support_indices)=...
    s.frequency_domain(g.signal_frequency_support_indices).*...
    g.filter;

fs.time_domain=ifft(fs.frequency_domain,p.number_points_frequency_domain);
fs.time_domain=fs.time_domain(1:p.number_points_time_domain);
