%Fig.3.1 
function fig31 
clear all 
  
%Constants (all MKS, except energy which is in eV) 
hbar=1.06e-34;
q=1.6e-19;
epsil=10*8.85E-12;
kT=.025; 
m=.25*9.1e-31;
n0=2*m*kT*q/(2*pi*(hbar^2)); 
  
%inputs 
a=3e-10;
t=(hbar^2)/(2*m*(a^2)*q);
beta=q*a*a/epsil; 
Ns=15;
Nc=70;
Np=Ns+Nc+Ns;
XX=a*1e9*[1:1:Np]; 
mu=.318;
Fn=mu*ones(Np,1); 
Nd=2*((n0/2)^1.5)*Fhalf(mu/kT) 
Nd=Nd*[ones(Ns,1);.5*ones(Nc,1);ones(Ns,1)]; 
  
%d2/dx2 matrix for Poisson solution 
D2=-(2*diag(ones(1,Np)))+(diag(ones(1,Np-1),1))+(diag(ones(1,Np-1),-1)); 
D2(1,1)=-1;
D2(Np,Np)=-1; %zero field condition 
  
%Hamiltonian matrix 
T=(2*t*diag(ones(1,Np)))-(t*diag(ones(1,Np-1),1))-(t*diag(ones(1,Np-1),-1)); 
  
%energy grid 
NE=301;
E=linspace(-.25,.5,NE);
dE=E(2)-E(1);
zplus=i*1e-12; 
f0=n0*log(1+exp((mu-E)./kT)); 
  
%self-consistent calculation 
U=[zeros(Ns,1);.2*ones(Nc,1);zeros(Ns,1)];%guess for U 
dU=0;
ind=10;
while ind>0.01 
    %from U to n 
    sig1=zeros(Np);sig2=zeros(Np);n=zeros(Np,1); 
    for k=1:NE 
        ck=1-((E(k)+zplus-U(1))/(2*t));
        ka=acos(ck); 
        sig1(1,1)=-t*exp(i*ka);
        gam1=i*(sig1-sig1'); 
        ck=1-((E(k)+zplus-U(Np))/(2*t));
        ka=acos(ck); 
        sig2(Np,Np)=-t*exp(i*ka);gam2=i*(sig2-sig2'); 

        G=inv(((E(k)+zplus)*eye(Np))-T-diag(U)-sig1-sig2); 
        A=i*(G-G');
        rhoE=f0(k)*diag(A)./(2*pi); 
        n=n+((dE/a)*real(rhoE)); 
    end 
  
    %correction dU from Poisson 
    D=zeros(Np,1); 
    for k=1:Np 
        z=(Fn(k)-U(k))/kT; 
        D(k,1)=2*((n0/2)^1.5)*((Fhalf(z+.1)-Fhalf(z))/.1)/kT; 
    end 
    dN=n-Nd+((1/beta)*D2*U); 
    dU=(-beta)*(inv(D2-(beta*diag(D))))*dN;U=U+dU; 

    %check for convergence 
    ind=(max(max(abs(dN))))/(max(max(Nd))) 
end 

plot(XX,n) 
%plot(XX,U) 
%plot(XX,F,':') 
  
save fig31 
hold on 
  
%For periodic boundary conditions, add 
 %           T(1,Np)=-t;T(Np,1)=-t; 
%and replace section entitled "from U to n" with 
            %from U to n 
            %[P,D]=eig(T+diag(U));D=diag(D); 
            %rho=log(1+exp((Fn-D)./kT));rho=P*diag(rho)*P'; 
            %n=n0*diag(rho);n=n./a; 
  


