function [logl grad Hess Grho Gbeta]=G_H(mle,fe,YL,YR,X)




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
D=(R.*R-1).*B;
E=YL-F;
EB=E.*B;
EC=E.*C;
ED=E.*D;
H=EC-EB.*EB;



J=ED-3*EB.*EC+2*EB.*EB.*EB;




DFERHO     =-ones(T,1)*(sum(YR.*H)./sum(H));
DFEBETA    =-ones(T,1)*(sum(X.*H) ./sum(H));






DFERHORHO  = ones(T,1)*((sum(J.*(YR+DFERHO)).*sum(YR.*H)-sum(YR.*J.*(YR+DFERHO)).*sum(H))./((sum(H)).^2));
DFEBETABETA= ones(T,1)*((sum(J.*(X+DFEBETA)).*sum(X.*H) -sum(X .*J.*(X+DFEBETA)).*sum(H))./((sum(H)).^2));
DFERHOBETA = ones(T,1)*((sum(J.*(X+DFEBETA)).*sum(YR.*H)-sum(YR.*J.*(X+DFEBETA)).*sum(H))./((sum(H)).^2));





logl=mean(mean(YL.*logF+(1-YL).*logA));
grad=[mean(mean(EB.*(YR+DFERHO)));mean(mean(EB.*(X+DFEBETA)))];
Hess=[mean(mean(H.*(YR+DFERHO).^2+EB.*DFERHORHO))           , mean(mean(H.*(YR+DFERHO).*(X+DFEBETA)+EB.*DFERHOBETA));
      mean(mean(H.*(YR+DFERHO).*(X+DFEBETA)+EB.*DFERHOBETA)), mean(mean(H.*(X+DFEBETA).^2+EB.*DFEBETABETA))        ];





Grho = EB.*(YR+DFERHO);
Gbeta= EB.*(X+DFEBETA);


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
