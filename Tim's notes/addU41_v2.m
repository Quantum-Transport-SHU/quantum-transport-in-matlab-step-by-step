% differences from v_1
% directly get U from Poisson equation
% not adding U0

function addU41_v2
clear all; 

% Constants (all MKS, except energy which is in eV) [[Page 261
hbar=1.06e-34;
q=1.6e-19;
epsil=10*8.85E-12;
kT=.025;
m=.25*9.1e-31;
n0=2*m*kT*q/(2*pi*(hbar^2));
a=3e-10;
t=(hbar^2)/(2*m*(a^2)*q);
beta=q*a*a/epsil; %??? q^2/epsil
Ns=15;  % sample points # of source
Nc=70;
Np=Ns+Nc+Ns;
XX=a*1e9*[1:1:Np];
mu=.318;
Nd=2*((n0/2)^1.5)*Fhalf(mu/kT); % ionized donor density of n-yped semiconductor (above 1.1)
Nd=Nd * [ones(Ns,1); .5*ones(Nc,1); ones(Ns,1)];

% d2/dx2 matrix for Poisson solution 
D2=-(2*diag(ones(1,Np)))+(diag(ones(1,Np-1),1))+(diag(ones(1,Np-1),-1));
D2(1,1)=-1; D2(Np,Np)=-1; % zero field condition 

% Hamiltonian matrix 
T=(2*t*diag(ones(1,Np)))-(t*diag(ones(1,Np-1),1))-(t*diag(ones(1,Np-1),-1)); % (2.1)
Jop=(q*t/((Np-1)*hbar))*i*((diag(ones(Np-1,1),-1))-(diag(ones(Np-1,1),1)));

% energy grid 
NE=301;
E=linspace(-.25,.5,NE);
dE=E(2)-E(1);
zplus=i*1e-12;
f0=n0*log(1+exp((mu-E)./kT));  % fermi function

% initial guess for U
U=[zeros(Ns,1);.2*ones(Nc,1);zeros(Ns,1)];
% U0=[0; 0; 0; .1; .1; .1; .1; .1; zeros(2*Ns + Nc - 8, 1)]; %%%%%%%%%%%%
% % U0 = zeros(2*Ns + Nc, 1);

% voltage bias steps 
NV=5;
VV=linspace(0,.25,NV);

for kV=1:NV 
    V=VV(kV) 
    Fn=[mu*ones(Ns,1); (mu-(.5*V))*ones(Nc,1); (mu-V)*ones(Ns,1)];
    f1=n0*log(1+exp((mu-E)./kT));
    f2=n0*log(1+exp((mu-V-E)./kT));

    in=10; 
    while in>.01 
        sig1=zeros(Np);
        sig2=zeros(Np);
        rho=zeros(Np);
        for k=1:NE 
            ck=1-((E(k)+zplus-U(1))/(2*t));
            ka=acos(ck); 
            sig1(1,1)=-t*exp(i*ka); % (3.15a)
            gam1=i*(sig1-sig1'); % (4.4)
            ck=1-((E(k)+zplus-U(Np))/(2*t));
            ka=acos(ck);
            sig2(Np,Np)=-t*exp(i*ka); % (3.15a)
            gam2=i*(sig2-sig2'); % (4.4)
            G=inv(((E(k)+zplus)*eye(Np))-T-diag(U)-sig1-sig2); % (4.3) %%%%%%%%%%%%
            A1=G'*gam1*G;
            A2=G'*gam2*G;
            rho=rho+(dE*((f1(k)*A1)+(f2(k)*A2))/(2*pi)); % (4.1)
        end

        n=(1/a)*real(diag(rho)); % (1.6)

        % correction dU from Poisson 
        % D=zeros(Np,1); 
        % for k=1:Np 
        %     z=(Fn(k)-U(k) - U0(k))/kT;
        %     D(k,1)=2*((n0/2)^1.5)*((Fhalf(z+.1)-Fhalf(z))/.1)/kT; 
        % end
        % dN=n - (Nd - (1/beta)*D2*(U)); %%%%%%%%%%%%
        % dU=(-beta)*(inv(D2-(beta*diag(D))))*dN;
        % U=U+dU; % (1.1)

        dN=n - (Nd - (1/beta)*D2*(U));
        U = beta * inv(D2) * (Nd - n);

        % check for convergence 
        in=(max(max(abs(dN))))/(max(max(Nd)))
    end 
    UU(:,kV)=U; 
    J(:,kV)=(-.5*q)*diag((rho*Jop)+(Jop*rho)); 
end 

II=sum(J); 
  
% save fig41 
hold on 
%%
% plot(VV,II) 

plot(U)