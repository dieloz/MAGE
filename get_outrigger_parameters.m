function ...
    [gabor_parameters,difference_in_parameters]=get_outrigger_parameters(varargin)

%     [gabor_parameters,difference_in_parameters]=get_outrigger_parameters(...
%                      'focus_gabor',focus_gabor,...
%                      'ip_magnitude',ip_magnitude);
%
% INPUTS --
%   focus_gabor: struct from make_chirplet.m 
%   parameters, [center_time; center_frequency; duration_parameter]
%   ip_magnitude:   real-valued number in range [0,1)

% OUTPUTS --
%   gabor_parameters: 3x8 matrix of gabor parameters, each column one gabor
%         eg, gabor_parameters(1,2)== center time of outrigger gabor_2
%             gabor_parameters(2,2)== center frequency of outrigger gabor_2
%   difference_in_parameters: 3x8 matrix, gabor_parameters -
%                   repmat(focus_gabor_parameters,1,8)

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
%         case 'focus_gabor_parameters'
%             focus_gabor_parameters=varargin{n+1};
        case 'focus_gabor'
            focus_gabor=varargin{n+1};
        case 'ip_magnitude'
            ip_magnitude=varargin{n+1};
    end
end

% t0=focus_gabor_parameters(1);
% v0=focus_gabor_parameters(2);
% s0=focus_gabor_parameters(3);

t0=focus_gabor.center_time;
v0=focus_gabor.center_frequency;
s0=focus_gabor.duration_parameter;

focus_gabor_parameters=[t0 v0 s0]';

dmat=zeros(3,8);
temp1=exp(s0/2);
temp2=sqrt(2/pi);
temp3=sqrt(log(1/ip_magnitude));
dmat(1,1)=temp1*temp2*temp3;
dmat(1,2)=-dmat(1,1);
dmat(1,3)=dmat(1,1)/2;
dmat(1,4)=dmat(1,2)/2;
temp4=exp(-s0/2)*sqrt(3/(2*pi));
dmat(2,3)=temp3*temp4;
dmat(2,4)=dmat(2,3);
dmat(1,5)=dmat(1,3);
dmat(1,6)=dmat(1,4);
dmat(2,5)=-dmat(2,3);
dmat(2,6)=-dmat(2,4);
dmat(3,7)=2*asech(ip_magnitude^2);
dmat(3,8)=-2*asech(ip_magnitude^2);

gabor_parameters=dmat+repmat(focus_gabor_parameters,1,8);
difference_in_parameters=dmat;






























