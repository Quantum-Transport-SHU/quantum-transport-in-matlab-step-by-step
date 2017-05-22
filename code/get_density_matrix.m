function D = get_density_matrix(S, H, Sigma, fermi, kt, voltage)

E_min = -20+fermi;
E_max = 20+fermi;
tol = 1e-8;

D = my_quad(@(x)get_negf(x, S, H, Sigma, fermi, kt, voltage), E_min, E_max, tol);

D = 1.j*squeeze(D)/2/pi;

D = real(D);
end


