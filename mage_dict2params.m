function ...
    gparams = mage_dict2params(gdict,p)
 
% gparams = mage_dict2params(gdict,p)

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


Natoms=length(gdict);
gparams=zeros(Natoms,5);
for n=1:Natoms
    clear g
    g=gdict{n};
    gparams(1,n)=g.center_time;
    gparams(2,n)=g.center_frequency;
    gparams(3,n)=g.duration_parameter;
    gparams(4,n)=abs(g.wt);
    gparams(5,n)=angle(g.wt);
end
