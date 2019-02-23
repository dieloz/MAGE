function ...
    mage_output=mage_decomp(raw_signal,Natoms,mage_parameters)

%     mage_output=mage_decomp(raw_signal,Natoms,mage_parameters)

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


raw_signal=raw_signal(:);

% initialize output structure
clear mage_output
mage_output.gparams=zeros(Natoms,5,'single');
mage_output.raw_signal=single(raw_signal);
mage_output.synth_signal=zeros(size(raw_signal),'single');
mage_output.residual=single(raw_signal);
mage_output.decay=zeros(Natoms,1,'single');
mage_output.SNR=zeros(Natoms,1,'single');
mage_output.Natoms=Natoms;

 p=get_signal_parameters(...
      'sampling_rate',mage_parameters.srate,...
      'number_points_time_domain',length(mage_parameters.center_time_support));
mage_output.time_support=p.time_support;
mage_output.srate=p.sampling_rate;
  
% perform MAGE decomp
for iA=1:Natoms
    disp(Natoms-iA+1);
    tfs_array=filter_with_gabor_filterbank1(...
        mage_output.residual,...
        mage_parameters.srate,...
        mage_parameters.center_frequency_support,...
        mage_parameters.fbw_support);
    [best_val,best_ind]=max(abs(tfs_array(:)));
    best_est_params=[...
        mage_parameters.tt_array(best_ind) ...
        mage_parameters.vt_array(best_ind) ...
        mage_parameters.st_array(best_ind) ...
        best_val ...
        angle(tfs_array(best_ind))];

    % given best parameter values from MP, now reassign via MAGE
    temp=mage_reassign_gparams_td(...
        mage_output.residual,...
        mage_parameters.srate,...
        best_est_params);
    % only reassign if amplitude of reassigned point is larger
    if temp(4)>best_est_params(4)
        best_est_params=temp;
    end
    
    mage_output.gparams(iA,:)=best_est_params;
    
    rgksignal=mage_params2rsignal(best_est_params,p);
    rgksignal=rgksignal(:);
    
    mage_output.synth_signal=mage_output.synth_signal+rgksignal;
    mage_output.residual=mage_output.residual-rgksignal;
    mage_output.decay(iA)=mage_output.gparams(end,4);
    mage_output.SNR(iA)=norm(mage_output.synth_signal)/norm(mage_output.residual);
end

