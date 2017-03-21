function fig52 
clear all 
  
%Constants (all MKS, except energy which is in eV) 
hbar=1.06e-34;q=1.6e-19;epsil=10*8.85E-12;kT=.025; 
m=.25*9.1e-31;n0=2*m*kT*q/(2*pi*(hbar^2)); 
  
%inputs 
a=3e-10;t=(hbar^2)/(2*m*(a^2)*q);beta=q*a*a/epsil; 
Ns=15;Nc=70;Np=Ns+Nc+Ns;XX=a*1e9*[1:1:Np]; 
mu=.318;Nd=2*((n0/2)^1.5)*Fhalf(mu/kT) 
            Nd=Nd*[ones(Ns,1);.5*ones(Nc,1);ones(Ns,1)]; 
  
%d2/dx2 matrix for Poisson solution 
D2=-(2*diag(ones(1,Np)))+(diag(ones(1,Np-1),1))+...
(diag(ones(1,Np-1),-1)); 
D2(1,1)=-1;D2(Np,Np)=-1;%zero field condition 
  
%Hamiltonian matrix 
T=(2*t*diag(ones(1,Np)))-(t*diag(ones(1,Np-1),1))-...
(t*diag(ones(1,Np-1),-1)); 
Jop=(q*t/((Np-1)*hbar))*i*((diag(ones(Np-1,1),-1))-...
(diag(ones(Np-1,1),1))); 
  
%energy grid 
NE=301;E=linspace(-.25,.5,NE);dE=E(2)-E(1),zplus=i*1e-12; 
f0=n0*log(1+exp((mu-E)./kT)); 
  
%initial guess for U 
U=[zeros(Ns,1);.2*ones(Nc,1);zeros(Ns,1)]; 
  
%voltage bias steps 
NV=5;VV=linspace(0,.25,NV);dV=VV(2)-VV(1);Fn=mu*ones(Np,1); 
for kV=1:NV 
V=VV(kV) 
f1=n0*log(1+exp((mu-E)./kT));f2=n0*log(1+exp((mu-V-E)./kT)); 
  
            in=10; 
            while in>.01 
                        sig1=zeros(Np);sig2=zeros(Np);rho=zeros(Np); 
  
sigs=-i*.0125*[zeros(Ns,1);ones(Nc,1);zeros(Ns,1)];sigs=diag(sigs); 
                        gams=i*(sigs-sigs');gams=diag(gams); 
                        for k=1:NE 
                                    fs=n0*log(1+exp((Fn-E(k))./kT)); 
                                    sigin=fs.*gams;sigin=diag(sigin); 
                                    ck=1-((E(k)+zplus-U(1))/(2*t));ka=acos(ck); 
                                    sig1(1,1)=-t*exp(i*ka);gam1=i*(sig1-sig1'); 
                                    ck=1-((E(k)+zplus-U(Np))/(2*t));ka=acos(ck); 
                                    sig2(Np,Np)=-t*exp(i*ka);gam2=i*(sig2-sig2'); 
  
G=inv(((E(k)+zplus)*eye(Np))-T-diag(U)-sig1-sig2-sigs); 
                                    A1=G'*gam1*G;A2=G'*gam2*G; 
  
            rho=rho+(dE*((f1(k)*A1)+(f2(k)*A2)+(G'*sigin*G))/(2*pi)); 
                        end 
                                                n=(1/a)*real(diag(rho)); 
  
                        %correction dU from Poisson 
                        D=zeros(Np,1); 
                        for k=1:Np 
                                    z=(Fn(k)-U(k))/kT; 
D(k,1)=2*((n0/2)^1.5)*((Fhalf(z+.1)-Fhalf(z))/.1)/kT; 
                        end 
                                    dN=n-Nd+((1/beta)*D2*U); 
dU=(-beta)*(inv(D2-(beta*diag(D))))*dN;U=U+dU; 
  
                        %check for convergence 
                        in=(max(max(abs(dN))))/(max(max(Nd))) 
            end 
            UU(:,kV)=U; 
            J(:,kV)=(-.5*q)*diag((rho*Jop)+(Jop*rho)); 
            Fn=Fn+[zeros(1,Ns) linspace(0,-dV,Nc) -dV*ones(1,Ns)]'; 
end 
II=sum(J); 
Fn=Fn-[zeros(1,Ns) linspace(0,-dV,Nc) -dV*ones(1,Ns)]'; 
  
save fig52 
hold on 
plot(XX,J(:,NV)) 
%plot(XX,Fn) 
%plot(XX,UU(:,NV)) 
  
