function [PredRotations, volConstructed] = Using_Sychronization_V0_14( inputFileName, outputMRCFileName, apix, refqFileName)
% The Input file should be the file produced by Preprocess_CompCls_V0_14.m
% The output MRC file contains the constructed volume
% apix is pixA
% refqFileName is optional, if provided, it is only used to align the predicted rotations
%
% The 3D reconstruction using Aspire 0.14 by Yoel Shkolnisky, Jan. 2019.
%
% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2019/11

load(inputFileName, 'noisy_projs','clstack','corrstack', 'shift_equations', 'n_theta');

[n, ~, K] = size(noisy_projs);

est_shifts=[shift_equations(:,1:end-1);sparse(1:3,1:3,ones(1,3),3,2*K);eye(2*K)]...
    \[shift_equations(:,end);0;0;0;zeros(2*K,1)];
est_shifts = full(reshape(est_shifts, 2, K)');

%% Estimate orientations using sychronization.

log_message('Estimating projection angles using Sychronization');

% S=cryo_syncmatrix_vote(clstack,n_theta);
% PredRotations=cryo_syncrotations(S); % old version 2013

% Using the new Sychronization method of 2016
PredRotations=cryo_sync3n_estimate_rotations(clstack,n_theta,0);

if (nargin>3)
    load(refqFileName, 'refq');
    K=size(refq,2);
    refRots = qs_to_rots(refq, 1);
    newRefRots = refRots;
    J2=[0 1 0; 1 0 0; 0 0 1]; % matrix to exchange x and y axis
    for i = 1:K
        newRefRots(:,:,i) = refRots(:,:,i)*J2;
    end
    PredRotations = Align2Rots(PredRotations, newRefRots); 
end

%%  3D reconstruction
maxIter = 30;

log_message('Reconstructing 3D density from projections');

vf = recon3d_firm( noisy_projs,...
PredRotations, -est_shifts, 1e-6, maxIter, zeros(n,n,n));

volConstructed=real(vf);

WriteMRC(volConstructed, apix, outputMRCFileName); % Output density map

end 

