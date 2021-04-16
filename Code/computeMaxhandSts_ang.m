function [maxh_matrix, AngleALL]=computeMaxhandSts_ang(c,corr,N)
% Compute angle between two projections by voting
% Reference: 
% Singer, A., Coifman, R.R., Sigworth, F.J. et al., "Detecting Consistent Common Lines
% in Cryo-EM by Voting". Journal of Structural Biology. 2010, 169(3): 312-322.
%
% By Zhang, Bianlan
%
% Modified by Yonggang Lu (ylu@lzu.edu.cn)
%

maxh_matrix = zeros(N,N);
AngleALL = zeros(N,N);

% compute similarity matrix
T=60;
tho=180/T;

for  k1 = 1:N-1
    for k2 = k1+1:N   
         h =zeros(1,T);
         angles = -1*ones(1,N);
         for k3 = 1:N
             if(k3==k1)||(k3==k2)
                 continue;
             end
		     [angle12, flag] = ComputeAngleCorr(c(k1,k2),c(k2,k1),corr(k1,k2),c(k1,k3),c(k3,k1),corr(k1,k3),c(k2,k3),c(k3,k2),corr(k2,k3));
             angles(k3) = angle12;
             
             if flag == true
                ang = (1:T)*180.0/T; 
                h = h + exp(-(ang-angle12).^2./(2*tho*tho));  % Modified by Lu    
             end
         end
         
         h = h./N; % Modified by Lu

         [maxh,idx] = max(h);
         ang = idx*tho - tho/2;
         
         % only use the top triangle
         AngleALL(k1, k2) = ang;
         maxh_matrix(k1, k2) = maxh;
        
    end
end

% use top triangle to creat the whole symmetric matrix -ylu
AngleALL =  AngleALL+AngleALL';
maxh_matrix = maxh_matrix + maxh_matrix';

end
