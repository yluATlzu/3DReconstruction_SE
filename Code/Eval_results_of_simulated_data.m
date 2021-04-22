function [MSE_angle, refDist, predDist] = Eval_results_of_simulated_data(refqFileName, PredRotations)
% Evaluation of the projection angle results of the simulated data
% refqFileName - file stores reference rotations by quaternions
% PredRotations - file stores predicted rotation matrices

% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2020/11

% loading reference data
load(refqFileName, 'refq');
K=size(refq,2);
refRots = qs_to_rots(refq, 1);
newRefRots = refRots;
J2=[0 1 0; 1 0 0; 0 0 1]; % matrix to exchange x and y axis
for i = 1:K
    newRefRots(:,:,i) = refRots(:,:,i)*J2;
end

[ref_angles, refDist, refDistXaxes] = computeDistFromRots(newRefRots);

% Computing predicted data
[~, predDist, predDistXaxes]= computeDistFromRots(PredRotations);

%% Evaluation 
% Evaluating alpha beta 
% figure; plot(refDist(:), predDist(:),'.'); 
% title(['Evaluating alpha beta for ', PredictionResultFileName], 'Interpreter','none');

% Evaluating gamma 
% figure; plot(refDistXaxes(:), predDistXaxes(:), '.'); 
% title(['Evaluating gamma for ', PredictionResultFileName], 'Interpreter','none');

% Evaluating all angles 
flagInv =1;
alignedRotations = Align2Rots(PredRotations, newRefRots); 
pred_angles_aligned = rots_to_EulerAngles (alignedRotations, flagInv); 
tmp=mod(abs(pred_angles_aligned-ref_angles), 2*pi);
tmp(tmp>pi)=2*pi-tmp(tmp>pi); 
MSE_angle = mean(tmp.^2,1);

disp(['MSE of Euler Angles is ',  num2str(MSE_angle)]);

end


