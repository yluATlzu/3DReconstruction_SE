function allEulerAngles = FindGammaAngles(OutputAngles)
% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2018/11

K = size(OutputAngles, 1)/2;

% find alpha betta
projAngles = OutputAngles(1:K,:);

% compute gamma
pred_Xaxes_new = OutputAngles(K+1:2*K,:);

c1=pred_Xaxes_new(:, 1);  %theta
p1=pred_Xaxes_new(:, 2);  %phai
v1 = [sin(c1).*cos(p1), sin(c1).*sin(p1), cos(c1)]; %nX3

R = calcuRotationMatrix([projAngles(:,1), projAngles(:,2), zeros(K,1)]);

gammas = zeros(K,1);
for i=1:K
    rotatedXaxes = R(:,:,i)'*v1(i,:)';
    gammas(i) = cart2pol(rotatedXaxes(1), rotatedXaxes(2));    
end

gammas(gammas<0)=gammas(gammas<0)+2*pi;

allEulerAngles = [projAngles gammas];

end
    



