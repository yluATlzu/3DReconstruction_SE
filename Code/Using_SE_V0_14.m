function [PredRotations, volConstructed, votedAngle] = Using_SE_V0_14(inputFileName, outputMRCFileName, apix, refqFileName)
% The Input file should be the file produced by Preprocess_CompCls_V0_14.m
% The output MRC file contains the constructed volume
% apix is pixA
% refqFileName is optional, if provided, it is only used to align the predicted rotations
%
% The 3D reconstruction using Spherical Embeddings.
%
% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2019/11

load(inputFileName, 'noisy_projs','clstack','corrstack', 'shift_equations');

[n, ~, K] = size(noisy_projs);

est_shifts=[shift_equations(:,1:end-1);sparse(1:3,1:3,ones(1,3),3,2*K);eye(2*K)]...
    \[shift_equations(:,end);0;0;0;zeros(2*K,1)];
est_shifts = full(reshape(est_shifts, 2, K)');

%% Compute Projection Angles from images (!! Independent from Aspire !!)

log_message('Estimating projection angles using Spherical Embedding');

[PredRotations,  votedAngle] = ComputeAnglesFromProjs_SE_Ver7(clstack, corrstack);

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

