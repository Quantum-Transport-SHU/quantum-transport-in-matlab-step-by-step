clear;
h0 = 0;
t0 = -1.0;
H0 = [h0,t0,0,0,0; t0,h0,t0,0,0; 0,t0,h0,t0,0; 0,0,t0,h0,t0; 0,0,0,t0,h0];
% 电极在后�?第二�?coupling _index �?

coupling_sample = 15;
v_sample = 71;
currents = zeros(coupling_sample, v_sample);
couplings = linspace(0.0001, .3, coupling_sample);
voltages = linspace(-2.0,2.0, v_sample)/2;
for i =  1:numel(couplings)
    extends = {
        struct('part_type', 0, 'epsilon', 0, 'coupling', couplings(i), 'coupling_index', 2, 'sigma', 0, 'fermi', 0),
        struct('part_type', 1, 'epsilon', 0, 'coupling', 0.1, 'coupling_index', 1, 'sigma', .1*1.j, 'fermi', 0), 
        struct('part_type', 1, 'epsilon', 0, 'coupling', 0.1, 'coupling_index', 5, 'sigma', .1*1.j, 'fermi', 0)
    };


    kt = 0.1;
    i
    currents(i, :) = sc_iv_curve(voltages, H0, extends, kt) * .0066236 % To A
end

surf(voltages, couplings,  currents)