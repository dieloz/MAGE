function ...
    gp_params=const_aip_shell_probe_gabors(gt_params,aip,N)

%  gp_params=const_aip_shell_probe_gabors(gt_params,aip,N);

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

% get const-aip shell gabor 3params in terms of focus 3params and spherical
% coordinates (r,theta,phi) MATH CONVENTION, NOT PHYSICS (SEE WIKIPEDIA)


offset=2/N;
increment=pi*(3-sqrt(5));
for n=1:N
    y = -1 + (n-(1/2))*offset;
    r = sqrt(1-y^2);
    phi = increment * mod(n,N);
    x = r*cos(phi);
    z = r*sin(phi);    
end
r=sqrt(x.^2+y.^2+z.^2);
theta=atan2(y,x);
phi=acos(z./r);
% reassign value of r so it is a function of aip
r=sqrt(2)*sqrt(log(1./(aip)));

% change of coordinates
% from spherical coords (egocentric relative to center target (tt,vt,st)
% to standardized (essentially z-score normalized) coords (allocentric frame of reference)
zds=r.*cos(phi);
zdt=r.*cos(theta).*sin(phi);
zdv=r.*sin(theta).*sin(phi);
temp=zds./sqrt(zds.^2);
ds=2*temp.*acosh(exp(zds.^2));
st=gt_params(3);
ip_tstd=sqrt((1/(2*pi))*(exp(st).*(1+exp(ds))));
ip_vstd=sqrt((1/(2*pi))*(exp(-st-ds).*(1+exp(ds))));
dt=zdt.*ip_tstd;
dv=zdv.*ip_vstd;

gp_params=zeros(N,3);
gp_params(:,3)=gt_params(3)+ds;
gp_params(:,2)=gt_params(2)+dv;
gp_params(:,1)=gt_params(1)+dt;

