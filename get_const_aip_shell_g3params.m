function ...
    [t_outrigger,v_outrigger,s_outrigger,params,dparams]=...
    get_const_aip_shell_g3params(t_focus,v_focus,s_focus,aip,num_points)

%     [t_outrigger,v_outrigger,s_outrigger,params,dparams]=...
%     get_const_aip_shell_g3params(t_focus,v_focus,s_focus,aip,num_points)

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

%     [t_outrigger,v_outrigger,s_outrigger,params,dparams]=...
%     get_const_aip_shell_g3params(t_focus,v_focus,s_focus,aip,num_points)


% get const-aip shell gabor 3params in terms of focus 3params and spherical
% coordinates (r,theta,phi) MATH CONVENTION, NOT PHYSICS (SEE WIKIPEDIA)
pp=fibonacci_sphere_sph(num_points);
theta=pp(:,2);
phi=pp(:,3);

% value of r is function of aip
r=sqrt(2)*sqrt(log(1./(aip)));

% change of coordinates from spherical coords (egocentric relative to
% focus), to standardized (essentially z-score normalized) coords
% (allocentric frame of reference)
[zdt,zdv,zds]=sph2z(r,theta,phi);

ds=zds2ds(zds);
% [ip_tstd,ip_vstd]=spds2ipstd(ds,s_focus);
[dt,dv]=z2dtdv(zdt,zdv,zds,s_focus);

s_outrigger=s_focus+ds;
t_outrigger=t_focus+dt;
v_outrigger=v_focus+dv;

params=[t_outrigger v_outrigger s_outrigger];

dparams=[dt dv ds];
