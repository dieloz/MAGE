function ...
    gsignal=mage_params2signal(gparams,p)
  
% gsignal=mage_params2signal(gparams,p)

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


gsignal=zeros(1,p.number_points_time_domain);
t=p.time_support;
Ntarget=size(gparams,1);
for k=1:Ntarget
    tk=gparams(k,1);
    vk=gparams(k,2);
    sk=gparams(k,3);
    Ak=gparams(k,4);
    rpk=gparams(k,5);
    
    tstd=sqrt(exp(sk)/(4*pi));
    ind_start=find(p.time_support >= tk-6*tstd,1,'first');
    ind_stop=find(p.time_support <= tk+6*tstd,1,'last');
    inds=ind_start:ind_stop;
    t=p.time_support(inds);
    
    gk=exp(-exp(-sk)*pi*(t-tk).^2).*... % amplitude envelope
        exp(2*pi*1i*vk*(t-tk)); % frequency modulation
    gk=gk/norm(gk); % need to normalize due to discrete sampling
    wgk=Ak*exp(1i*rpk)*gk;
    
    gsignal(inds)=gsignal(inds)+wgk;
end

gsignal=gsignal(:);
