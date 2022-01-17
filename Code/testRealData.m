function testRealData(realdata_name)
% Run experiments with the simulated data sets 
% the following parameter is used to select the data set
%  realdata_name - can be 'EMPIAR-10028' or 'EMPIAR-10328'

% Author: Yonggang Lu (ylu@lzu.edu.cn)
% Updated 2022/01

%% Produce the reconstruction results

if strcmp(realdata_name, 'EMPIAR-10028')
    apix = 2.68;
    refmapFnm = 'Ref_map_emd_2660.mrc';
    Preprocess_CompCls_V0_14( '10028_531p.mrcs',  'RealProjectionCLstack.mat', 5, false);
elseif strcmp(realdata_name, 'EMPIAR-10328')
    apix = 1.059;
    refmapFnm = 'emd_22689.map';
    Preprocess_CompCls_V0_14( '10328_390p.mrcs',  'RealProjectionCLstack.mat', 5, false);
end

pause(5);
Using_LUD_V0_14( 'RealProjectionCLstack.mat', 'Mrc_LUD_real.mrc', apix);
pause(5);
Using_Sychronization_V0_14( 'RealProjectionCLstack.mat', 'Mrc_Sychronization_real.mrc', apix);
pause(5);
Using_SE_V0_14( 'RealProjectionCLstack.mat', 'Mrc_SE_real.mrc', apix);
pause(5);

%% Evaluation with FSC

[fscLUD, ~] = Eval_w_FSC(refmapFnm, 'Mrc_LUD_real.mrc', apix);
[fscSych, ~] = Eval_w_FSC(refmapFnm, 'Mrc_Sychronization_real.mrc', apix);
[fscSE, freqs] = Eval_w_FSC(refmapFnm, 'Mrc_SE_real.mrc', apix);

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

end
