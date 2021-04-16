function testSimulatedData(numPtcls, snr)
% Run experiments with the simulated data sets 
% the following two parameters is used to select the data set
%  numPtcls - number of images
%  snr - SNR level

% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2021/04

%% Produce the 3D reconstruction results

refqFileName = ['refq', num2str(numPtcls), '.mat'];

inputfileName = ['EMD3508_', num2str(numPtcls), 'particles@SNR=', num2str(snr), '.mrcs'];
CLfileName = ['CL_', num2str(numPtcls), 'p@SNR=', num2str(snr), '.mat'];

mrcLUDfileName = ['MrcLUD-', num2str(numPtcls), 'p@SNR=', num2str(snr), '.mrc'];
mrcSychfileName = ['MrcSychronization-', num2str(numPtcls), 'p@SNR=', num2str(snr), '.mrc'];
mrcSEfileName = ['MrcSE-', num2str(numPtcls), 'p@SNR=', num2str(snr), '.mrc'];

% Compute common lines
if ~isfile(CLfileName)
    Preprocess_CompCls_V0_14( inputfileName, CLfileName, 0);
    pause(5);
end

% Compute 3D reconstructions 
PredRotationsLUD = Using_LUD_V0_14( CLfileName, mrcLUDfileName, 1.084, refqFileName);
pause(5);
PredRotationsSync = Using_Sychronization_V0_14( CLfileName, mrcSychfileName, 1.084, refqFileName);
pause(5);
PredRotationsSE = Using_SE_V0_14( CLfileName, mrcSEfileName, 1.084, refqFileName);
pause(5);

%% Evaluation of the predicted angles 

[mseLUD, refDist, predDistLUD] = Eval_results_of_simulated_data(refqFileName, PredRotationsLUD);
[mseSych, ~, predDistSych] = Eval_results_of_simulated_data(refqFileName, PredRotationsSync);
[mseSE, ~, predDistSE] = Eval_results_of_simulated_data(refqFileName, PredRotationsSE);


%% Evaluation with FSC

[fscLUD, ~] = Eval_w_FSC('ref_map_emd_3508.mrc', mrcLUDfileName, 1.084, true);
[fscSych, ~] = Eval_w_FSC('ref_map_emd_3508.mrc', mrcSychfileName, 1.084, true);
[fscSE, freqs] = Eval_w_FSC('ref_map_emd_3508.mrc', mrcSEfileName, 1.084, true);

outputfileName = ['Results_All_', num2str(numPtcls), 'p@SNR=', num2str(snr), '.mat'];

save(outputfileName, 'refDist', 'predDistLUD', 'predDistSych', 'predDistSE', ...
    'mseLUD', 'mseSych', 'mseSE', 'freqs', 'fscLUD', 'fscSych', 'fscSE');


%% plot the FSC curves
figure; plot(freqs, fscSE, 'k');
hold;
plot(freqs, fscLUD, 'b');
plot(freqs, fscSych, 'g');
ns = size(freqs, 2);
plot(freqs, ones(ns,1)*0.143, 'r--');
plot(freqs, ones(ns,1)*0.5, 'r--');
title(inputfileName);
legend('SE_OurVer7','LUD', 'Syc');

end