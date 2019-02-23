function ...
    boolval=gns_in_gk_neighborhood(tk,vk,sk,tns,vns,sns,aipthresh)

% boolval=gns_in_gk_neighborhood(tk,vk,sk,tns,vns,sns,aipthresh)

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


dt=tns-tk;
dv=vns-vk;
ds=sns-sk;

sub2=(exp(-sk)*pi)./(1+exp(ds));
sub3=(exp(ds+sk)*pi)./(1+exp(ds));
p1=sqrt(sech(ds/2));
p2=exp((-dt.^2).*sub2);
p3=exp((-dv.^2).*sub3);
aipval=p1.*p2.*p3;

boolval=(aipval>=aipthresh);

