function ...
    ip=ip2g(tt,vt,st,tp,vp,sp)

%     ip=ip2g(tt,vt,st,tp,vp,sp);
% 
% computes inner product between two Gabor time-frequency atoms (Gaussian envelope)
% 
% inputs --
% tt: 1xN vector, center time of target Gabor atom
% vt: 1xN vector, center frequency of target Gabor atom
% st: 1xN vector, duration parameter of target Gabor atom
% tp: 1xN vector, center time of probe Gabor atom
% vp: 1xN vector, center frequency of probe Gabor atom
% sp: 1xN vector, duration parameter of probe Gabor atom
% 
% output --
% ip: inner product between target and probe Gabor atoms (complex)

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

% % correct orientation:
% if size(tp,1)>size(tp,2)
%     tp=tp';
% end
% if size(vp,1)>size(vp,2)
%     vp=vp';
% end
% if size(sp,1)>size(sp,2)
%     sp=sp';
% end
% if size(tt,1)>size(tt,2)
%     tt=tt';
% end
% if size(vt,1)>size(vt,2)
%     vt=vt';
% end
% if size(st,1)>size(st,2)
%     st=st';
% end
% 
% match array sizes:
if isequal(size(tt,1),[1 1])
    tt=repmat(tt(1),size(tp));
end
if isequal(size(vt,1),[1 1])
    vt=repmat(vt(1),size(tp));
end
if isequal(size(st,1),[1 1])
    st=repmat(st(1),size(tp));
end

if isequal(size(vp,1),[1 1])
    vp=repmat(vp(1),size(tp));
end
if isequal(size(sp,1),[1 1])
    sp=repmat(sp(1),size(tp));
end

% % inner product of 2 Gabor atoms:
% set chirp rates to zero, compute inner product for Gaussian chirplets:
ct=0;
cp=0;
dt=tp-tt;
ps=sp+st;
p1=4*1i*pi*dt.^2;
p2=exp(sp).*(1i*ps-4*pi*dt.*(cp.*dt-2*vp));
p3=exp(st).*(1i*ps+4*pi*dt.*(ct.*dt+2*vt));
p4=1i*exp(ps).*...
    (-1i*ct.*ps-8*pi*ct.*dt.*vp+4*pi*(vp-vt).^2+...
    cp.*(1i*ps+4*pi*dt.*(ct.*dt+2*vt)));
p5=-4*1i*exp(sp)+4*exp(st).*(-1i+(cp-ct).*exp(sp));
ip=sqrt(2)*exp((p1+p2+p3+p4)./p5)./...
    sqrt(1i*(cp-ct)+exp(-sp)+exp(-st));



