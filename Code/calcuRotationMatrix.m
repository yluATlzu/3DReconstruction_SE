function [ Rotation ] = calcuRotationMatrix( angle )
%CALCUROTATIONMATRIX Summary of this function goes here
%   Compute rotation from 2D to 3D
%
% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2017/11

N = size(angle,1);
Rotation = zeros(3,3,N);
for i =1:N
    
    a = angle(i,2);
    b = angle(i,1);
    r = angle(i,3);
    
    Ra = [cos(a),sin(a),0;-sin(a),cos(a),0;0,0,1];
    Rb = [cos(b),0,-sin(b);0,1,0;sin(b),0,cos(b)];
    Rr = [cos(r),sin(r),0;-sin(r),cos(r),0;0,0,1];
    
    Rotation(:,:,i) = Rr*Rb*Ra;

end

% !!! Compute the transpose of the rotation- Now rotates from 2D to 3D !!
for i=1:N
    Rotation(:,:,i) = Rotation(:,:,i)';
end

end

