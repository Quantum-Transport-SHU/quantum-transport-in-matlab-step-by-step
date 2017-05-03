function D = get_density_matrix(S, H, Sigma, fermi, kt, voltage)

E_min = -5;
E_max = 5;
tol = 1e-4;

D = my_quad(@(x)get_negf(x, S, H, Sigma, fermi, kt, voltage), E_min, E_max, tol)/(2*pi);

D = squeeze(D);



