function [H,U,rho] = get_selfconsistent_hamiltonian(S, H0, Sigma, D2, Laplacian, voltage, fermi, kt)

N = length(S);
dU = zeros(N,1) + 1;
tol = 1e-6;

D0 = get_density_matrix(S, H0, Sigma, fermi, kt, 0);

H = H0;
H0_length = length(H0);

n=0;
U0 = zeros(N,1) + 1;

mx = mixer(0.1, 6, false);

addition = zeros(2,1);
addition(1) = voltage/2;
addition(2) = -voltage/2;
rho_addition = -Laplacian(:,[H0_length+1,H0_length+2])*addition;
rho_addition = rho_addition(1:H0_length);

while max(abs(dU)) > tol
    D = get_density_matrix(S, H, Sigma, fermi, kt, voltage);
    drho = -4*pi*diag(D-D0);

    U = inv(D2)*(drho + rho_addition);

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
rho = -4*pi*diag(D);
end



