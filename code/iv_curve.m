function currents=iv_curve(voltages, mol_H, extends)
%% This function output an i-v curve of a transport model including
%% electrodes and local states, the iv curve is obtained from the
%% integration of transmission coefficients without electrostatic field
%% selfconsistency 

%% input requirement: 
%% voltages : vector of bias voltage
%% mol_H: a model hamiltonian of the central system
%% extends: the enviroment of the 'mol_H', including electrodes and local states whose number is not limited.
%%          each cell element of 'extends' is a local state or electrode with its members defined by the following.
%%          part_type:  1 for electrode, 0 for local state
%%          epsilon:   onsite energy
%%          coupling:  vector of the coupling hamiltonian with the molecule
%%          coupling_index : index of the coupling, meaning to which orbital or grid it couples with.
%%          sigma:  selfenergy, pure imaginary matrix with the same dimension as 'epsilon'
%%          fermi: fermi level of the electrode, not useful for local state


n0 = length(mol_H);
N = length(extends);
H = zeros(n0+N, n0+N);

nlead = 0;
for i=1:N
    if extends{i}.part_type == 1
        nlead = nlead+1;
    end
end

Sigma = zeros([nlead, n0+N, n0+N]);
Gamma = zeros([nlead, n0+N, n0+N]);

H(1:n0,1:n0) = mol_H;
ilead = 0;
for i = 1:N
    H(n0+i, n0+i) = extends{i}.epsilon;
    ind = extends{i}.coupling_index;
    H(ind, n0+i) = extends{i}.coupling';
    H(n0+i, ind) = extends{i}.coupling;
    if extends{i}.part_type == 1
        H(n0+i, n0+i) = H(n0+i, n0+i) + extends{i}.sigma;
        ilead = ilead + 1;
        Sigma(ilead,n0+i,n0+i) = Sigma(ilead,n0+i,n0+i) + extends{i}.sigma;         
    end   
end

for i =1:nlead
    Gamma(i,:,:) = 1.j*(Sigma(i,:,:) - Sigma(i,:,:)');
end
    
maxv = max(abs(voltages));
fermi = extends{1}.fermi;
NE = 201;
energies = linspace(fermi-maxv/2, fermi+maxv/2, NE);
eta = 1e-6;
S = eye(n0+N, n0+N);
T = zeros(nE);
for i=1:NE
    G = inv((energies(i)+1j*eta)*S-H);
    T(i) = real(trace(squeeze(Gamma(1,:,:))*G*squeeze(Gamma(2,:,:))*G'));
end

nv = len(voltages);
currents = zeros(nv);
dE = energies(2) - energies(1);
for i = 1:nv
    factors = fermidistribution(energies, fermi+voltages(i)/2, 0) - fermidistribution(energies, fermi-voltages(i)/2, 0);
    currents(i) = sum(T*factors)*dE;
end

end
























