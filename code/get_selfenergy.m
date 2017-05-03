function sigma=get_selfenergy(energy, S_00, S_01, H_00, H_01)
tau_ij = energy*S_01-H_01;
tau_ji = energy*S_01'-H_01';
g=get_surface_green_function(energy, S_00, S_01,H_00,H_01);
sigma = tau_ji*g*tau_ij;
end
