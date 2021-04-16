function [alignedTargetRot, O, mse, flag] = Align2Rots(targetRot, refRot) 
% Align  targetRot to refRot (both are inverse rotations)
% Rotation matrices R1 and R2 must have the same size: 3X3XK 
% targetRot[3,3,N], refRot[3,3,N]
% 
% O: rotation matrix for the target
% flag = 2 means needs to do mirro tranlation before rotation
%
% Author: Yonggang Lu (ylu@lzu.edu.cn)
% 2017/11

TOL=1.0e-14;
J=[1 0 0; 0 1 0; 0 0 -1]; % Reflection matrix
J2=[0 1 0; 1 0 0; 0 0 1]; % matrix to exchange x and y axis
J3=[-1 0 0; 0 1 0; 0 0 1]; % Reflection matrix 2
J4=[-1 0 0; 0 1 0; 0 0 -1]; % Reflection matrix 3

K=size(refRot, 3);

% Register estimated rotations to the true one, and compute the difference
    % between the two.
    % modified rot2 -- Ylu, 10/26/07
    % Added rot3 and rot4-- Ylu, 10/30/07
    
    rot=zeros(3*K,3);  % The K estimated rotation matrices stacked as a matrix with 3K rows.
    rot1=zeros(3*K,3); % True true K rotation matrices stacked as a matrix with 3K rows.
    rot2=zeros(3*K,3); % Reflected matrices of rot, which are also an estimated rotation matrices.
    rot3=zeros(3*K,3); % matrices of rot exchanged x and y axis, which are also an estimated rotation matrices.
    rot4=zeros(3*K,3); % Reflected matrices of rot exchanged x and y axis, which are also an estimated rotation matrices.
    rot5=zeros(3*K,3); 
    rot6=zeros(3*K,3); 
    
    rot7=zeros(3*K,3); 
     
    for k=1:K
        R=targetRot(:,:,k);
        rot(3*(k-1)+1:3*k,:)=R';
        Rref=refRot(:,:,k);
        rot1(3*(k-1)+1:3*k,:)=Rref';
        rot2(3*(k-1)+1:3*k,:)=(R*J)';
        rot3(3*(k-1)+1:3*k,:)=(R*J2)';
        rot4(3*(k-1)+1:3*k,:)=(R*J*J2)';
        rot5(3*(k-1)+1:3*k,:)=(R*J2*J3)';
        rot6(3*(k-1)+1:3*k,:)=(R*J*J2*J3)';
        rot7(3*(k-1)+1:3*k,:)=(R*J3)';
        
        rot8(3*(k-1)+1:3*k,:)=(R*J4)';
    end
    
    % Compute the possible orthogonal matrices which register the
    % estimated rotations to the true ones.
    O1= rot.'*rot1./K;       
    O2=rot2.'*rot1./K;
    O3=rot3.'*rot1./K;
    O4=rot4.'*rot1./K;
    O5=rot5.'*rot1./K;
    O6=rot6.'*rot1./K; 
    O7=rot7.'*rot1./K;
    
    O8=rot8.'*rot1./K;
        
    err= zeros(1,7);
    err(1)=norm(O1*O1.'-eye(3));
    err(2)=norm(O2*O2.'-eye(3));
    err(3)=norm(O3*O3.'-eye(3));    
    err(4)=norm(O4*O4.'-eye(3));  
    err(5)=norm(O5*O5.'-eye(3)); 
    err(6)=norm(O6*O6.'-eye(3));
    err(7)=norm(O7*O7.'-eye(3));
    
    err(8)=norm(O8*O8.'-eye(3));
    
    % Find the best O.
    [minErr, flag] = min(err);
    
    if (minErr>TOL) 
        fprintf('Registering matrix is not orthogonal, err=%e  tol=%e\n',...
            minErr,TOL);
    end
    
    errd=zeros(1,7);  
    errd(1)=abs(det(O1)-1);
    errd(2)=abs(det(O2)-1);
    errd(3)=abs(det(O3)-1);
    errd(4)=abs(det(O4)-1);
    errd(5)=abs(det(O5)-1);
    errd(6)=abs(det(O6)-1);
    errd(7)=abs(det(O7)-1);
    
    errd(8)=abs(det(O8)-1);
    
    [minErrd, flagTemp] = min(errd);
    
    % Adjust the selected flag --  Ylu, 10/31/07
    if (err(flagTemp)==err(flag))
        flag = flagTemp;
    else
        minErrd = errd(flag);
    end
   
    if (minErrd>TOL) 
        fprintf('Determinant of registering matrix is not 1, err=%e  tol=%e\n',...
            minErrd,TOL);
    end
    
    
    if flag == 1
        [U,~,V]=svd(O1); % Use O1 as the registering matrix
    elseif flag ==2
        [U,~,V]=svd(O2); % Use O2 as the registering matrix
    elseif flag ==3
        [U,~,V]=svd(O3); % Use O3 as the registering matrix 
    elseif flag ==4
        [U,~,V]=svd(O4); % Use O4 as the registering matrix
    elseif flag ==5
        [U,~,V]=svd(O5); % Use O5 as the registering matrix   
    elseif flag ==6
        [U,~,V]=svd(O6); % Use O5 as the registering matrix 
    elseif flag ==7
        [U,~,V]=svd(O7); % Use O5 as the registering matrix   
    else   % flag ==8
        [U,~,V]=svd(O8); % Use O6 as the registering matrix
    end   
    
    O=V*U';
    
    % Compute aligned rotation and estimation errors
    diff=zeros(K,1);
    mse=0;
    alignedTargetRot=zeros(3,3,K);
    for k=1:K
        R=targetRot(:,:,k);
        Rref=refRot(:,:,k);
        alignedRot = O*R;
        
        if (flag ==2 || flag ==4 || flag ==6)
            alignedRot = alignedRot*J;
        end
        
        if (flag ==3 || flag ==4 || flag ==5 || flag ==6) % modified by Ylu
            alignedRot = alignedRot*J2;
        end
        
        if (flag ==5 || flag ==6 || flag ==7)
            alignedRot = alignedRot*J3;
        end
        
        if (flag ==8)
            alignedRot = alignedRot*J4;
        end
        
        diff(k)=norm(O*R-Rref,'fro');
        mse=mse+diff(k).^2;
        
        alignedTargetRot(:,:,k)=alignedRot;  % modified by Ylu
    end
    mse=mse/K;
    
end
