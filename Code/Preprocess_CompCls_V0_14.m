function Preprocess_CompCls_V0_14( inputFileName,  outputFileName, max_shift)
% The input file is the MRC file containing the projection images
% The output file contains the estimated common lines and the the projection images
% max_shift=5 for real data and max_shift=0 for Simulated data
%
% To compute common lines between images using cross correlation
% Need ASPIRE V0.14
% 
% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2019/11

noisy_projs = ReadMRC(inputFileName);

[n, ~, K] = size(noisy_projs); % size of the volume: nXnXn

% normalization - YLU, added on 20/8/2018
meanSum = sum(noisy_projs(:))/K;
for i=1:K
    tmp=noisy_projs(:,:,i);
    rateSum=meanSum/sum(tmp(:));
    noisy_projs(:,:,i) = noisy_projs(:,:,i)*rateSum;
end
masked_projs=mask_fuzzy(noisy_projs,ceil(n/2)); % Applly circular mask

%% Compute polar Fourier transform, using radial resolution n_r and angular
% resolution n_theta. n_theta is the same as above.
n_theta=360;
n_r= size(noisy_projs,1);

[npf,sampling_freqs]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections   

% Find common lines from projections
shift_step=1;
[clstack,corrstack,shift_equations]= cryo_clmatrix(npf,-1,1,max_shift,shift_step);

save(outputFileName, 'noisy_projs','clstack','corrstack','n_theta','n_r', 'shift_equations');

end
