function ...
    mp_output=mp_decomp(raw_signal,Natoms,mp_parameters)

%     mp_output=mp_decomp(raw_signal,Natoms,mp_parameters)

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
clear mp_output
mp_output.gparams=zeros(Natoms,5,'single');
mp_output.raw_signal=single(raw_signal);
mp_output.synth_signal=zeros(size(raw_signal),'single');
mp_output.residual=single(raw_signal);
mp_output.decay=zeros(Natoms,1,'single');
mp_output.SNR=zeros(Natoms,1,'single');
mp_output.Natoms=Natoms;

 p=get_signal_parameters(...
      'sampling_rate',mp_parameters.srate,...
      'number_points_time_domain',length(mp_parameters.center_time_support));
mp_output.time_support=p.time_support;
mp_output.srate=p.sampling_rate;
  
% perform MP decomp
for iA=1:Natoms
    disp(Natoms-iA+1);
    tfs_array=filter_with_gabor_filterbank1(...
        mp_output.residual,...
        mp_parameters.srate,...
        mp_parameters.center_frequency_support,...
        mp_parameters.fbw_support);
    [best_val,best_ind]=max(abs(tfs_array(:)));
    best_est_params=[...
        mp_parameters.tt_array(best_ind) ...
        mp_parameters.vt_array(best_ind) ...
        mp_parameters.st_array(best_ind) ...
        best_val ...
        angle(tfs_array(best_ind))];
    
% % %     % given best parameter values from MP, now reassign via MAGE
% % %     temp=mp_reassign_gparams_td(...
% % %         mp_output.residual,...
% % %         mp_parameters.srate,...
% % %         best_est_params);
% % %     % only reassign if amplitude of reassigned point is larger
% % %     if temp(4)>best_est_params(4)
% % %         best_est_params=temp;
% % %     end
    
    mp_output.gparams(iA,:)=best_est_params;
    
    rgksignal=mage_params2rsignal(best_est_params,p);
    rgksignal=rgksignal(:);
    
    mp_output.synth_signal=mp_output.synth_signal+rgksignal;
    mp_output.residual=mp_output.residual-rgksignal;
    mp_output.decay(iA)=mp_output.gparams(end,4);
    mp_output.SNR(iA)=norm(mp_output.synth_signal)/norm(mp_output.residual);
end

