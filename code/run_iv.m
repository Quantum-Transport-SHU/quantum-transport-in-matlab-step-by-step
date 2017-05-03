h0 = 0;
t0 = -1.0;
H0 = [h0,t0,0,0,0; t0,h0,t0,0,0; 0,t0,h0,t0,0; 0,0,t0,h0,t0; 0,0,0,t0,h0];
% 电极在后， 第二个 coupling _index 为
extends = {
    struct('part_type', 0, 'epsilon', 0, 'coupling', 0.5, 'coupling_index', 5, 'sigma', 0, 'fermi', 2.3),
    struct('part_type', 1, 'epsilon', 0, 'coupling', 0.1, 'coupling_index', 1, 'sigma', .1*1.j, 'fermi', 2.3), 
    struct('part_type', 1, 'epsilon', 0, 'coupling', 0.1, 'coupling_index', 5, 'sigma', .1*1.j, 'fermi', 2.3)
};

N = 201;
voltages = linspace(-2.0,2.0, N);
%currents = iv_curve(voltages, H0, extends);

fermi = 1.0;
kt = 0.01;

currents = sc_iv_curve(voltages, H0, extends, fermi, kt);

plot(voltages, currents);
hold;
plot(voltages, -fliplr(currents), 'r');
plot(voltages, currents + fliplr(currents), 'g');