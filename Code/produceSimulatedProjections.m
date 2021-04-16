% To produce the simulated projections for emd_3508.map

clear;

% define the parameters
K = 100;   % number of images
SNR = 0.2;  % SNR level
apix = 1.084;
n=260; % size of the volume: nXnXn

tmp = ReadMRC('emd_3508.map');
volref = zeros(n, n, n, 'single');
volref(6:255,:,14:246)=tmp;

precision='single';
refq=qrand(K);  % Generate random uniform quaternions.

cleanprojs=cryo_project(volref, qs_to_rots(refq, 1), n, precision); 

noisy_projs=cryo_addnoise(cleanprojs, SNR, 'gaussian');

outputDatafilename = ['EMD3508_', num2str(K), 'particles@SNR=', num2str(SNR), '.mrcs'];

WriteMRC(noisy_projs, apix, outputDatafilename);
