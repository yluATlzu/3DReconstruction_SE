function target_angles_Aligned = MatchPerpedicularAngles(target_angles, ref_angles)
% Rotate target_angles so that they are mostly perpendicular to ref_angles
%
% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2018/11

n = size(target_angles, 1);
c1=target_angles(:, 1);  %theta
p1=target_angles(:, 2);  %phai
v1 = [sin(c1).*cos(p1), sin(c1).*sin(p1), cos(c1)]; %nX3

c2=ref_angles(:, 1);  %theta
p2=ref_angles(:, 2);  %phai
v2 = [sin(c2).*cos(p2), sin(c2).*sin(p2), cos(c2)]; %nX3

A=[  v1(:,1).*v2(:,1)  v1(:,2).*v2(:,1)  v1(:,3).*v2(:,1)];
A=[A v1(:,1).*v2(:,2)  v1(:,2).*v2(:,2)  v1(:,3).*v2(:,2)];
A=[A v1(:,1).*v2(:,3)  v1(:,2).*v2(:,3)  v1(:,3).*v2(:,3)];
B=zeros(n, 1);

[V, ~] =eig(A'*A); %
R=V(:,1); % rotation found using least square method
R=reshape(R, 3,3)';

perpendV1 = R*v1';

disp(['[MatchPerpedicularAngles] Average angle is: ' num2str(acos(sum(abs(dot(v2',perpendV1)))/n)*180/pi)]);

lengthV = (perpendV1(1,:).^2+perpendV1(2,:).^2+perpendV1(3,:).^2).^0.5;
perpendV1= perpendV1./repmat(lengthV, 3, 1);

alphas = cart2pol(perpendV1(1,:), perpendV1(2,:));
alphas(alphas<0)=alphas(alphas<0)+2*pi;
target_angles_Aligned= [real(acos(perpendV1(3,:)))' alphas'];

end
