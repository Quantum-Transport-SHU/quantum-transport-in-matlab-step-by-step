clear;
h0 = 0;
t0 = -1.0;
H0 = [h0,t0,0,0,0; t0,h0,t0,0,0; 0,t0,h0,t0,0; 0,0,t0,h0,t0; 0,0,0,t0,h0];
% 电极在后， 第二个 coupling _index 为
extends = {
    struct('part_type', 0, 'epsilon', 0, 'coupling', 0.1, 'coupling_index', 5, 'sigma', 0, 'fermi', 0),
    struct('part_type', 1, 'epsilon', 0, 'coupling', 0.1, 'coupling_index', 1, 'sigma', .1*1.j, 'fermi', 0), 
    struct('part_type', 1, 'epsilon', 0, 'coupling', 0.1, 'coupling_index', 5, 'sigma', .1*1.j, 'fermi', 0)
};

N = 41;
voltages = linspace(-2.0,2.0, N)/2;
%voltages = [-0.09];
%currents = iv_curve(voltages, H0, extends);

kt = 0.1;

currents = sc_iv_curve(voltages, H0, extends, kt);
%currents = iv_curve(voltages, H0, extends);

plot(voltages, currents);
hold;
plot(voltages, -fliplr(currents), 'r');
plot(voltages, currents + fliplr(currents), 'g');
