function ...
    synth_signal=mage_dict2csignal(gdict,p)

   % synth_signal=mage_dict2csignal(gdict,p)

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


synth_signal=zeros(1,p.number_points_time_domain);
Ntarget=length(gdict);
for n=1:Ntarget
    if mod(n,10)==0; disp(n); end
    
     clear g
     g=gdict{n};
    
    synth_signal(g.signal_time_support_indices)=...
        synth_signal(g.signal_time_support_indices)...
        +g.wt*g.time_domain;
        
end

