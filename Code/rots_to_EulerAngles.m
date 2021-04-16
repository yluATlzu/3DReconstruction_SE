function proj_angles = rots_to_EulerAngles (rots, flagInv)
% This only works for forward roation matrix as shown in (A.4) of the 1987 paper
% (The forward rotation is to rotate a plane's normal vector in 3D to the Z axis, 
% which corresponds to rotate a projection from 3D location to the XY plane)
% {The forward rotation can be used to produce sythetic projections at a
% specific orientation.}
% 
% For the inverse roation matrix , set flagInv =1 !
% The inverse roation matrix rotates a plane's normal vector from the Z axis 
% to the actual 3D orientation {used for 3D reconstruction} 
%
% If flagInv =1, tranpose rots first
%
% Output: proj_angles is the array of [beta, alpha, gamma]
%
% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2017/11

if ~exist('flagInv','var')
    flagInv=0;
end

N=size(rots,3);
proj_angles = zeros(N,3);   % 1 column is B; 2 column is A; 3 column is R
for i = 1:N
    
    rot=rots(:,:,i);
    
    if (flagInv)
        rot=rot';
    end
    
    B = acos(rot(3,3));
    
    sina = rot(3,2)/sin(B);
    cosa = rot(3,1)/sin(B);
    
    if(sina*cosa>0)
        A = asin(abs(sina));
        if(sina<0)
            A = A+pi;
        end
    else if(sina*cosa<0)
           A = acos(cosa);
           if(cosa>0)
               A = 2*pi -A;
           end
        else 
           if(sina~=0)
              A = asin(sina);
           else
              A = acos(cosa);
           end
         end
    end
    
    sinr = rot(2,3)/sin(B);
    cosr = -rot(1,3)/sin(B);
    
    if(sinr*cosr>0)
        R = asin(abs(sinr));
        if(sinr<0)
            R = R+pi;
        end
    else if(sinr*cosr<0)
           R = acos(cosr);
           if(cosr>0)
               R = 2*pi -R;
           end
        else 
           if(sinr~=0)
              R = asin(sinr);
           else
              R = acos(cosr);
           end
         end
    end
    
    proj_angles(i,:) =[B,A,R];
end
