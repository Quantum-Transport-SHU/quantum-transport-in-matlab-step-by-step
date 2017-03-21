clear all 
load fig53 
Jop=(q*t/hbar)*i*((diag(ones(Np-1,1),-1))-(diag(ones(Np-1,1),1))); 
  
V=VV(NV);U=UU(:,NV); 
f1=n0*log(1+exp((mu-E)./kT));f2=n0*log(1+exp((mu-V-E)./kT)); 
sig1=zeros(Np);sig2=zeros(Np);rhoE=zeros(Np); 
for k=1:NE 
            fs=n0*log(1+exp((Fn-E(k))./kT)); 
            sigin=fs.*gams;sigin=diag(sigin); 
            ck=1-((E(k)+zplus-U(1))/(2*t));ka=acos(ck); 
            sig1(1,1)=-t*exp(i*ka);gam1=i*(sig1-sig1'); 
            ck=1-((E(k)+zplus-U(Np))/(2*t));ka=acos(ck); 
            sig2(Np,Np)=-t*exp(i*ka);gam2=i*(sig2-sig2'); 
  
            G=inv(((E(k)+zplus)*eye(Np))-T-diag(U)-sig1-sig2-sigs); 
            A1=G'*gam1*G;A2=G'*gam2*G; 
            rhoE=((f1(k)*A1)+(f2(k)*A2)+(G'*sigin*G))/(2*pi); 
            JJ(:,k)=(-.5*q)*diag((rhoE*Jop)+(Jop*rhoE));k 
end 
 save fig53a 
hold on 
plot(JJ(2,:),E) 
plot(JJ(Np-1,:),E,'--') 
 
 