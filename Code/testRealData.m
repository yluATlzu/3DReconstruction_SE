% Run experiments with the real data set 

%% Produce the reconstruction results

Preprocess_CompCls_V0_14( '531p.mrcs',  'RealProjectionCLstack.mat', 5);

pause(5);
Using_LUD_V0_14( 'RealProjectionCLstack.mat', 'Mrc_LUD_real.mrc', 2.68);
pause(5);
Using_Sychronization_V0_14( 'RealProjectionCLstack.mat', 'Mrc_Sychronization_real.mrc', 2.68);
pause(5);
Using_SE_V0_14( 'RealProjectionCLstack.mat', 'Mrc_SE_real.mrc', 2.68);
pause(5);

%% Evaluation with FSC

[fscLUD, ~] = Eval_w_FSC('Ref_map_emd_2660.mrc', 'Mrc_LUD_real.mrc', 2.68);
[fscSych, ~] = Eval_w_FSC('Ref_map_emd_2660.mrc', 'Mrc_Sychronization_real.mrc', 2.68);
[fscSE, freqs] = Eval_w_FSC('Ref_map_emd_2660.mrc', 'Mrc_SE_real.mrc', 2.68);

save('FscCurves_realdata.mat', 'freqs', 'fscLUD', 'fscSE', 'fscSych');

figure; plot(freqs, fscSE, 'k');
hold;
plot(freqs, fscLUD, 'b');
plot(freqs, fscSych, 'g');
ns = size(freqs, 2);
plot(freqs, ones(ns,1)*0.143, 'r--');
plot(freqs, ones(ns,1)*0.5, 'r--');
title('FSCs for the Real datasets');
legend('SE-OurVer7','LUD', 'Syc');
