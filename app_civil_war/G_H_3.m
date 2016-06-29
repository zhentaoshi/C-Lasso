function [logl, grad, Hess, G1,G2,G3]=G_H_3(mle,fe,YL,YR,x2, x3)
% this script is developed upon Dhaene and Jochmans (2015)
[T,N]=size(YL);

R=ones(T,1)*fe'+mle(1)*YR+mle(2)*x2 + mle(3) * x3;
R(R<-4) = -4;
R(R>4) = 4;

F=normcdf(R); 
A=1-F; 
logF=log(F); 
logA=log(A);
logFA=logF+logA; 

logf=-0.5*(log(2*pi)+R.*R); 


B=exp(logf-logFA); 
C=-R.*B; 
D=(R.*R-1).*B;
E=YL-F; 
EB=E.*B; 
EC=E.*C; 
ED=E.*D; 
H=EC-EB.*EB; 
J=ED-3*EB.*EC+2*EB.*EB.*EB; % 3rd order 

% 2nd derivative. D-FE-1, D-FE-2, D-FE-3
DFE1    =-ones(T,1)*(sum(YR.*H)./sum(H));
DFE2    =-ones(T,1)*(sum( x2.*H ) ./sum(H));
DFE3    =-ones(T,1)*(sum( x3.*H ) ./sum(H));


DFE11  = ones(T,1)*((sum(J.*(YR+DFE1)).*sum(YR.*H)-sum(YR.*J.*(YR+DFE1)).*sum(H))./((sum(H)).^2));
DFE22  = ones(T,1)*((sum(J.*(x2+DFE2)).*sum(x2.*H)-sum(x2.*J.*(x2+DFE2)).*sum(H))./((sum(H)).^2));
DFE33  = ones(T,1)*((sum(J.*(x3+DFE3)).*sum(x3.*H)-sum(x3.*J.*(x3+DFE3)).*sum(H))./((sum(H)).^2));

DFE12  = ones(T,1)*((sum(J.*(x2+DFE2)).*sum(YR.*H)-sum(YR.*J.*(x2+DFE2)).*sum(H))./((sum(H)).^2));
DFE13  = ones(T,1)*((sum(J.*(x3+DFE3)).*sum(YR.*H)-sum(YR.*J.*(x3+DFE3)).*sum(H))./((sum(H)).^2));
DFE23  = ones(T,1)*((sum(J.*(x3+DFE3)).*sum(x2.*H)-sum(x2.*J.*(x3+DFE3)).*sum(H))./((sum(H)).^2));



logl=mean(mean(YL.*logF+(1-YL).*logA));

grad1 = mean(mean(EB.*(YR+DFE1)));
grad2 = mean(mean(EB.*(x2+DFE2)));
grad3 = mean(mean(EB.*(x3+DFE3)));
grad = [grad1; grad2; grad3];

H11 = mean(mean(H.*(YR+DFE1).^2+EB.*DFE11));
H22 = mean(mean(H.*(x2+DFE2).^2+EB.*DFE22));
H33 = mean(mean(H.*(x3+DFE3).^2+EB.*DFE33));
H12 = mean(mean(H.*(YR+DFE1).*(x2+DFE2)+EB.*DFE12));
H13 = mean(mean(H.*(YR+DFE3).*(x3+DFE3)+EB.*DFE13));
H23 = mean(mean(H.*(x2+DFE2).*(x3+DFE3)+EB.*DFE23));
H21 = H12; 
H31 = H13;
H32 = H23;


Hess = [H11, H12, H13; H21, H22, H23; H31, H32, H33];
Hess( isnan(Hess)) = 0;
Hess( isinf(Hess)) = 0;

G1 = EB.*(YR+DFE1);
G2 = EB.*(x2+DFE2);
G3 = EB.*(x3+DFE3);


function [logl grad Hess mle]=FELoglProbitARX1(fe,mle,YL,YR,X)

[T,N]=size(YL); 
R=ones(T,1)*fe'+mle(1)*YR+mle(2)*X; 
F=normcdf(R); 
A=1-F; 
logF=log(F); 
logA=log(A); 
logFA=logF+logA; 
logf=-0.5*(log(2*pi)+R.*R); 
B=exp(logf-logFA); 
C=-R.*B;
E=YL-F; 
EB=E.*B; 
EC=E.*C; 
H=EC-EB.*EB;

logl=sum(sum(YL.*logF+(1-YL).*logA)); 
grad=sum(EB)';
Hess=sum(H)';