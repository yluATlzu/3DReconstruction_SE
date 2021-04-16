function [angles, distAB, distXaxes] = computeDistFromRots(rotations)
% Help to compute useful distances from rotations
% rotations - input rotations
%
% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2020/11

    flagInv =1;
    angles = rots_to_EulerAngles(rotations, flagInv);
    distAB = sphricalDist(angles(:,1:2),angles(:,1:2));
    
    K= size(rotations, 3);
    x_axes = zeros(3,K);
    for i=1:K
        x_axes(:,i) = rotations(:,:,i)*[1,0,0]';
    end
    alphas = cart2pol(x_axes(1,:), x_axes(2,:));
    alphas(alphas<0)=alphas(alphas<0)+2*pi;
    x_axes_angles = [real(acos(x_axes(3,:)))' alphas'];
    distXaxes = sphricalDist(x_axes_angles, x_axes_angles);
end