function ...
    tfs_array=filter_with_gabor_filterbank1(rs,srate,cf_list,fbw_list)

%     tfs_array=filter_with_gabor_filterbank1(rs,srate,cf_list,fbw_list)

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

Nt=length(rs);
Nv=2^ceil(log2(Nt));

v=single((srate/Nv)*(0:Nv-1));
inds=(v>(srate/2));
v(inds)=v(inds)-srate;

rs_fd=fft(rs,Nv);
rs_fd(v<0)=0;
rs_fd(v>0)=2*rs_fd(v>0);

[vk_array,fbw_array]=ndgrid(cf_list,fbw_list);
g2params_list=[vk_array(:) fbw_array(:)];
Nfilters=size(g2params_list,1);
num_freqs=length(cf_list);
num_fbws=length(fbw_list);

gbank=zeros(Nv,num_freqs,num_fbws,'single');
for n=1:Nfilters
    vk=g2params_list(n,1);
    vk_ind=find(vk==cf_list);
    fbw=g2params_list(n,2);
    fbw_ind=find(fbw==fbw_list);
    
    sk=log((2*log(2))/(fbw^2*pi*vk^2));
    
    Gk=exp(-pi*exp(sk)*(v-vk).^2); % un-normalized
    Gk=Gk*(sqrt(length(v))/norm(Gk)); % because of discrete sampling and different time/freq sample numbers
    
    gbank(:,vk_ind,fbw_ind)=Gk;
end

signal_array=repmat(rs_fd(:),1,num_freqs,num_fbws); % should now be same size as gbank

fs_fd=signal_array.*gbank;
tfs_array=ifft(fs_fd,Nv,1);
tfs_array=tfs_array(1:Nt,:,:);


