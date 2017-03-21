clear all 
load fig41 
  
V=VV(NV);
U=UU(:,NV); 
f1=n0*log(1+exp((mu-E)./kT));
f2=n0*log(1+exp((mu-V-E)./kT)); 
sig1=zeros(Np);
sig2=zeros(Np);
rhoE=zeros(Np); 
for k=1:NE 
    ck=1-((E(k)+zplus-U(1))/(2*t));
    ka=acos(ck); 
    sig1(1,1)=-t*exp(i*ka);
    gam1=i*(sig1-sig1'); 
    ck=1-((E(k)+zplus-U(Np))/(2*t));
    ka=acos(ck); 
    sig2(Np,Np)=-t*exp(i*ka);
    gam2=i*(sig2-sig2'); 

    G=inv(((E(k)+zplus)*eye(Np))-T-diag(U)-sig1-sig2); 
    A1=G'*gam1*G;
    A2=G'*gam2*G; 
    rhoE=((f1(k)*A1)+(f2(k)*A2))/(2*pi); 
    JJ=(-.5*q)*diag((rhoE*Jop)+(Jop*rhoE)); 
    J1(k)=sum(JJ); 
    trans(k)=real(trace(gam1*A2)); 
    J2(k)=(q*q/(2*pi*hbar))*trans(k)*(f1(k)-f2(k)); 
    trans(k)=real(trace(gam2*A1)); 
    J3(k)=(q*q/(2*pi*hbar))*trans(k)*(f1(k)-f2(k)); 
end 
  
save fig41a 
figure;
hold on 
plot(J1,E) 
plot(J2,E,'x') 
plot(J3,E,'+') 