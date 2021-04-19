function produceSimulatedProjections(numPtcls, snr)
% To produce the simulated projections for emd_3508.map
% numPtcls is the number of images
% snr is the signal to noise ratio

% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2021/04

apix = 1.084;
n=260; % size of the volume: nXnXn

tmpmap = ReadMRC('emd_3508.map');
volref = zeros(n, n, n, 'single');
volref(6:255,:,14:246)=tmpmap;
if ~exist('ref_map_emd_3508.mrc', 'file')
    WriteMRC(volref, apix, 'ref_map_emd_3508.mrc'); 
end

refqFileName = ['refq', num2str(numPtcls), '.mat'];
if ~exist(refqFileName, 'file')
    refq=qrand(numPtcls);  % Generate random uniform quaternions.
    save(refqFileName, 'refq');
else
    load(refqFileName, 'refq');
end

precision='single';
cleanprojs=cryo_project(volref, qs_to_rots(refq, 1), n, precision); 
noisy_projs=cryo_addnoise(cleanprojs, snr, 'gaussian');

outputDatafilename = ['EMD3508_', num2str(numPtcls), 'particles@SNR=', num2str(snr), '.mrcs'];
WriteMRC(noisy_projs, apix, outputDatafilename);

end