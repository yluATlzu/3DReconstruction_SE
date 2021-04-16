function [angle12, flag] = ComputeAngleCorr(c12,c21,corr12,c13,c31,corr13,c23,c32,corr23,corrThresh, verbose)
% Compute angle between projection 1 and projection 2 using projection 3
% cij is the common line (angle) between projection i and projection j on
% projection i
%
% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2017/11

if (nargin == 9)
    verbose = false;
    corrThresh = -1;
end
if (nargin == 10)
    verbose = false;
end

% compute Euler angle between 1 and 2
if (corrThresh>0 && (corr12<=corrThresh || corr13<=corrThresh || corr23<=corrThresh)) 
    flag = false;
    angle12 = -1;    
    if (verbose) disp('Cannot find the angle (1)!'); end;
    return;
else    
    
    a=cos((c32-c31)/180*pi);  % c3
    b=cos((c23-c21)/180*pi);  % c2
    c=cos((c13-c12)/180*pi);  % c1

    if (1+2*a*b*c <= a*a+b*b+c*c )  
        flag = false;
        angle12 = -1;    
        if (verbose) disp('Cannot find the angle (2)!'); end;
        return;
    end
        
    cosAngle = (a-b*c)/((1-b*b)*(1-c*c))^0.5;
    
    % Modified by Zhang, Bianlan
    % cosAngle = abs(cosAngle); % removed
    if(c13>c12)
        if(c13-c12>180)
            cosAngle = -cosAngle;
        end
    else 
        if(c12-c13<180)
            cosAngle = -cosAngle;
        end
    end
    
    if(c23>c21)
       if(c23-c21>180)
            cosAngle = -cosAngle;
       end
    else
       if(c21-c23<180)
           cosAngle = -cosAngle;
       end
    end


    if (cosAngle<=1)
        flag = true;
        angle12 = real(acosd(cosAngle));
        if (verbose) disp(['The angle between 1 and 2 is: ', num2str(angle12)]); end;
    else
        flag = false;
        angle12 = -1;    
        if (verbose) disp('Cannot find the angle (3)!'); end;
        return;
    end
end

end
