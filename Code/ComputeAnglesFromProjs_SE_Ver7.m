%
% Computing Projection Angles from Input Images
%
% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2019/11

function [PredRotations, votedAngle]= ComputeAnglesFromProjs_SE_Ver7( clstack, corrstack)

%  compute maxh and votedAngle
K=size(clstack,2);
[maxh, ang] = computeMaxhandSts_ang(clstack,corrstack,K);
votedAngle = ang.*pi./180;

imgScore = sum(maxh);
[~, idx1] = max(imgScore);
[~, idx2] = max(maxh(idx1,:));
if (idx1>idx2) 
    tmp=idx1;
    idx1=idx2;
    idx2=tmp;
end

for i=1:K
  maxh(i,i)=1;
end

%% First spherical embedding 

p=10;

Pfilter1 = maxh>prctile(maxh(:),p);
Pfilter = Pfilter1.* maxh;

[pred_angles, tmpCost] = elliptic_embed_unitShpere(votedAngle, Pfilter, 2);

disp(['Final error  is: ', num2str(tmpCost)]);
predDist = sphricalDist(pred_angles, pred_angles);  


%% Compute gamma using spherical embedding with: "clstack", "predDist" (without "rotation")

clstackRad = clstack*pi/180;
deltaXaxes = zeros(K, K);
for i=1:K
    for j=1:K
        cosAngle = cos(clstackRad(i,j))*cos(clstackRad(j,i))+sin(clstackRad(i,j))*sin(clstackRad(j,i))*cos(predDist(i,j));
        deltaXaxes(i, j) = abs(acos(cosAngle));
    end
end

% Second spherical embedding of Xaxes
[pred_Xaxes, tmpCost] = elliptic_embed_unitShpere(deltaXaxes, Pfilter, 2);
disp(['Final error  is: ', num2str(tmpCost)]);

pred_Xaxes_Aligned = MatchPerpedicularAngles(pred_Xaxes, pred_angles);  


%% Compute rotation angles

% computing gamma before the third SNE
allEulerAngles = FindGammaAngles([pred_angles; pred_Xaxes_Aligned]);
allEulerAngles(:,3)= allEulerAngles(:,3)+1*pi/180;

PredRotations = calcuRotationMatrix( allEulerAngles );


%% selecting from two mirror rotations - Added by Ylu @ Dec. 06, 2019

J2=[-1 0 0; 0 1 0; 0 0 1];
for k=1:K
    PredRotations2(:,:,k)=PredRotations(:,:,k)*J2';
end

theta12=clstackRad(idx1, idx2);
cmLine12=PredRotations(:,:,idx1)*[cos(theta12), sin(theta12), 0]';
theta21=clstackRad(idx2, idx1);
cmLine21=PredRotations(:,:,idx2)*[cos(theta21), sin(theta21), 0]';

diffAngle1=norm(cmLine21-cmLine12);

theta12=clstackRad(idx1, idx2);
cmLine212=PredRotations2(:,:,idx1)*[cos(theta12), sin(theta12), 0]';
theta21=clstackRad(idx2, idx1);
cmLine221=PredRotations2(:,:,idx2)*[cos(theta21), sin(theta21), 0]';

diffAngle2=norm(cmLine221-cmLine212);

if abs(diffAngle2)<abs(diffAngle1)
    PredRotations=PredRotations2;   
end

end
