function H = get_selfconsistent_hamiltonian(S, H0, Sigma, D2, voltage, fermi, kt)

N = length(S);
dU = zeros(N,1) + 1;
tol = 1e-4;

D0 = get_density_matrix(S, H0, Sigma, fermi, kt, 0);

H = H0;

n=0;
U0 = zeros(N,1) + 1;

mx = mixer(0.01, 6, false);

while max(abs(dU)) > tol
    D = get_density_matrix(S, H, Sigma, fermi, kt, voltage);
    drho = -4*pi*diag(D-D0);
    addition = zeros(length(drho),1);
    addition(1) = voltage/2;
    addition(5) = -voltage/2;
    U = inv(D2)*(drho + addition);

    dU = U - U0;
    
    U = mx.mix(U);
    
    U0 = U;
    H = H0 + diag(U);   
    diff = max(abs(dU));
    n = n+1;
%     if diff > 1
%         fprintf('difference %f of step %d at voltage %f \n', diff,  n, voltage);
%     end
end

end



