function err=NormRMSError(Dn,D)
% err=NormRMSError(Dn,D) computes the RMS error of the distance matrix Dn
% with respect to D, after the distances are normalised so that each has
% unit average distance

sf=mean(mean(D));
ND=D/sf;

NDn=Dn/sf;

err=norm(ND-NDn,'fro')/size(D,1);
