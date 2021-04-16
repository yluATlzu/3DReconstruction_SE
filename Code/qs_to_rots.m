function rots = qs_to_rots(qs, flagInv)
% From quaternions to rotation matrices.
%
% flagInv (optional)
% flagInv=1: do the transpose/inverse
% flagInv=0 (default): do not do the transpose/inverse
%
% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2018/11

if ~exist('flagInv','var')
    flagInv=0;
end

N= size(qs, 2);
rots = zeros(3, 3, N);
for i=1:N
    rots(:,:,i) = q_to_rot(qs(:,i));
    if (flagInv)
        rots(:,:,i)=rots(:,:,i)';
    end
end

end