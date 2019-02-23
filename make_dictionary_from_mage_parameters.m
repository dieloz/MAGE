function ...
    gdict = make_dictionary_from_mage_parameters(gabor_params,p)

%   gdict = make_dictionary_from_mage_parameters(gabor_params,p)

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


if size(gabor_params,1)~=5
    gabor_params=gabor_params';
end

% numparams=size(gabor_params,1);
Natoms=size(gabor_params,2);
for n=1:Natoms
    clear g
    g.center_time=gabor_params(1,n);
    g.center_frequency=gabor_params(2,n);
    g.duration_parameter=gabor_params(3,n);
    g.chirp_rate=0;
    g=make_chirplet(...
        'chirplet_structure',g,...
        'signal_parameters',p);
    
    wt=gabor_params(4,n)*exp(1i*gabor_params(5,n));
    g.wt=wt;
    gdict{n}=g;
end





