function ...
    fs=gabor_filter2_fbw(rs,srate,vk,fbw)

% function ...
%     fs=gabor_filter2_fbw(rs,srate,vk,fbw);
%
% input --
%     rs: 1xN real-valued time series of raw signal
%     srate: number of samples per second
%     vk: mean frequency value of gaussian envelope of Gabor time-frequency atom
%     fbw: bandwidth parameter, frequency-domain standard deviation proportional to vk*fbw
%
% output --
%     fs: 1xN complex-valued time series of filtered signal

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


sk=log((2*log(2))/(fbw^2*pi*vk^2));

Nt=length(rs);
Nv=2^ceil(log2(Nt));

v=single((srate/Nv)*(0:Nv-1));
inds=(v>(srate/2));
v(inds)=v(inds)-srate;

rs_fd=fft(rs,Nv);
rs_fd(v<0)=0;
rs_fd(v>0)=2*rs_fd(v>0);

Gk=exp(-pi*exp(sk)*(v-vk).^2); % un-normalized
Gk=Gk*(sqrt(length(v))/norm(Gk)); % because of discrete sampling and different time/freq sample numbers

fs_fd=rs_fd.*Gk;
fs=ifft(fs_fd,Nv);
fs=fs(1:Nt);

end







