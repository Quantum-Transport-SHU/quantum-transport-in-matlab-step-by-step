function H = get_selfconsistent_hamiltonian(S, H0, Sigma, voltage, fermi, kt)

N = length(S);
dU = zeros(N) + 1;
tol = 1e-4;

D0 = get_density_matrix(S, H0, Sigma, fermi, kt, 0);

H = H0;

D2=-(2*diag(ones(1,N)))+(diag(ones(1,N-1),1))+(diag(ones(1,N-1),-1)); 

while max(abs(dU)) > tol
    D = get_density_matrix(S, H, Sigma, fermi, kt, voltage);
    drho = -4*pi*diag(D-D0);
    addition = zeros(len(drho));
    addition(1) = voltage/2;
    addition(5) = -voltage/2;
    dU = D2\(drho + addition);
    H = H + diag(dU);
end

end



