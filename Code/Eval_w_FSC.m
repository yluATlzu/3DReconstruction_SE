function [fscCurve, freqs] = Eval_w_FSC(RefVolFileNM, ConstructedVolFileNM, apix, flagAligned)
% !! NEED set to actual APIX
%  apix = 1.084 for DataEMD3508
%  apix = 2.68 for the real dataset
%  flagAligned is optional, default is false; if true, no density alignment is performed 

% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2020/11

if (nargin<4)
    flagAligned = false;
end

volref = ReadMRC(RefVolFileNM);
volConstructed = ReadMRC(ConstructedVolFileNM);

if (not(flagAligned))
    % Rotate and align the images 
    [Rest,estdx,volConstructedAlign,reflect]=cryo_align_densities(volref,volConstructed);
    WriteMRC(volConstructedAlign, apix, ['Aligned-', ConstructedVolFileNM]);
else
    volConstructedAlign = volConstructed;
end

% FSC curves
fscCurve = FSCorr(volConstructedAlign, volref);

n=size(volConstructed, 3);

ns = floor(n/2);
freqs = [1:ns]/(n*apix);

% figure; plot(freqs, fscCurve, 'b');
% hold;
% plot(freqs, ones(ns,1)*0.143, 'r--');
% plot(freqs, ones(ns,1)*0.5, 'r--');
% title(['FSC curve for ', ConstructedVolFileNM]);

end
