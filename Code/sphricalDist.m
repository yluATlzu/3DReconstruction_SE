function dist = sphricalDist(ydata1, ydata2)
% compute distance on spherical surface
%   ydata(:,1) are polar angles, ydata(:,2) are azimuthal angles 
%   ydata1 needs to be (mX2), ydata2 can be (nX2)  
%
% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2017/11

n = size(ydata1, 1);
m = size(ydata2, 1);
cosdist = cos(ydata1(:,1))*cos(ydata2(:,1)');
cosdist = cosdist + (sin(ydata1(:,1))*sin(ydata2(:,1)').*cos(repmat(ydata1(:,2), 1, m)-repmat(ydata2(:,2)',n,1)));
dist = real(acos(cosdist));

end

