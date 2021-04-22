function [ydata,e]=elliptic_embed_unitShpere(D, weights, x, initSolution)
% The function will attempt to find an embedding of the points, whose distance matrix 
% is given by D, into a unit spherical surface in (x+1)-dimensional space. 
% The surface itself is actually x-dimensional.
%
% D: input distance matrix
% weights: weights[i,j] gives the weight for D[i,j]
% x: dimension of spherical surface
% initSolution: initial solution
%
% The Matlab code is a modified version of the code for the spherical and 
% hyperbolic embedding described in the paper:
% "Spherical and Hyperbolic Embeddings of Data", Richard C. Wilson, 
% Edwin R. Hancock, Elzbieta Pekalska, Robert P.W. Duin,
% IEEE Transactions on Pattern Analysis and Machine Intelligence, 2014
%
% The optimization method used is an efficient optimisation method based on gradient descent
% on the tanget plane.

% Modified by Yonggang Lu (ylu@lzu.edu.cn) @ 2018/4/5

sz=size(D,1);

r=1;
edim = x;

if ~exist('initSolution', 'var')
    
    Z=cos(D/r);
    Z=0.5*(Z+Z');

    % Correction to eigenvalues
    [U,E]=eig(Z);
    if( ~issorted(diag(E)) )
        error('Eigenvalues not in sorted order');
    end

    v=diag(E)-1;
    mv=v(sz-edim-1);
    if mv>=-1
        d=diag(E);
    else
        d=1-(diag(E)-1)/mv;
    end
    d(1:sz-edim-1)=0;
    Z=U*diag(d)*U';

    % vector lengths should be 1 at this point, but need renormalising if we
    % are using reduced dimensions
    S=diag(1./sqrt(diag(Z)));
    Z=S*Z*S';

    % result should be symmetric, this just ensures it is precisely so
    Z=0.5*(Z+Z');
    
    r=1;
    [U,E]=eig(Z);
    X=r*U(:,sz-edim:sz)*sqrt(max(0,E(sz-edim:sz,sz-edim:sz)));
else
    X=initSolution;
end
    

%% Optimisation - this can be time consuming, could be improved?

disp('      Iteration       Error');

Z=X;
errOld = Inf;

MaxIter = 1000;
for iter=1 : MaxIter
    ZZ=cell(1,sz);
    for i=1:sz
        x=X(i,:);
        for j=1:sz
            y=X(j,:);
            ZZ{i}(j,:)=LogMapE(y,x);
        end
    end
        
    
    for i=1:sz
        x=X(i,:);

        Z=ZZ{i};
        idx=1:sz;
        idx(i)=[];
        [dx,resnorm]=doptimise(Z(idx,:),D(idx,i),r*edim*0.0001, weights(i,idx));
        
        y=ExpMapE(dx,x);
        
        X(i,:)=y;
    end
    
    if ~rem(iter, 10) 
        Dc=r*acos(max(-1,min(1,X*X'/r^2)));
        err=norm(weights.*(D-Dc),'fro');
        disp(['Iteration ' num2str(iter) ': error is ' num2str(err)]);
        
        if (abs(err/errOld-1) < 1e-6 ) 
            disp(['Stop early at Iteration ' num2str(iter) ': error is ' num2str(err)]);
            break;
        end
        errOld = err;
    end
end

Z = X*X'/r^2;

%% 

% value clipping, Z should always be in range [-1,1]
Z=min(1,max(-1,0.5*(Z+Z')));

% some data about the quality of the embedding
Dn=acos(Z);
e=NormRMSError(Dn,D);

%  Calculate and plot embedding coordinates
[U,E]=eig(0.5*(Z+Z'));
eX=U*sqrt(max(0,E));
V1 = eX(:,sz-2:sz); % x,y,z

alphas = cart2pol(V1(:,1), V1(:,2));
alphas(alphas<0)=alphas(alphas<0)+2*pi;
ydata= [real(acos(V1(:,3))) alphas];

 
%% Optimization functions
function [nx,resnorm]=doptimise(Y,d,TolX, weights)
[samples,dim]=size(Y);
nx=zeros(1,dim);
dsq=zeros(samples,1);

resnorm=0;
for iter=1:1
    for p=1:samples
        xd=nx-Y(p,:);
        dsq(p)=xd*xd';
    end
    
    g=4*(weights.*(dsq-d.^2)')*(repmat(nx,samples,1)-Y);
    dx=g;
    
    c=zeros(1,4);
    dE2=dx*dx';
    c(1)=4*samples*dE2.^2;

    for p=1:samples
        
        xd=(nx-Y(p,:));
       
        fp=(xd*xd'-d(p)^2)*weights(p);
        gp=dx*xd';
                
        c(2)=c(2)+12*gp*dE2;
        c(3)=c(3)+8*gp^2+4*fp*dE2;
        c(4)=c(4)+4*fp*gp;
        
    end
    
    c(1) = max(realmin, c(1));     
    roots=RealCubicRoot(c);
    step=max(roots(imag(roots)==0 & roots<0));
    if (~isempty(step))
        nx=nx+step(1)*dx;
    end

    if norm(step*dx)<TolX
        return;
    end
end


function x=RealCubicRoot(c)
A=[ -c(2)/c(1) -c(3)/c(1) -c(4)/c(1); 1 0 0; 0 1 0 ];
x=eig(A);

%% Log and exp maps for sphere

function x=LogMapE(v,m)
theta=acos(max(-1,min(1,dot(v,m)/dot(m,m))));
if abs(theta)<1e-8
    x=zeros(size(v));
    return;
end
x=theta*(v-m*cos(theta))/sin(theta);

function v=ExpMapE(x,m)
r=sqrt(dot(m,m));
theta=norm(x)/r;
if abs(theta)<1e-8
    v=m;
    return;
end
v=x*sin(theta)/theta+m*cos(theta);
