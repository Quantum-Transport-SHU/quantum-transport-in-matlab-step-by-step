clear;
h0 = 0.5;
t0 = -0.1;

% H0_dimension = 10;
% 
% H0 = zeros(H0_dimension, H0_dimension);
% 


H0 = [h0,t0,0,0,0; t0,h0,t0,0,0; 0,t0,h0,t0,0; 0,0,t0,h0,t0; 0,0,0,t0,h0];
% ????????????????????????coupling _index ???
extends = {
    struct('part_type', 0, 'epsilon', 0.5, 'coupling', -0.05, 'coupling_index', 3, 'sigma', 0, 'fermi', 0),
    struct('part_type', 1, 'epsilon', 0.5, 'coupling', -0.1, 'coupling_index', 1, 'sigma', 0.1*1.j, 'fermi', 0), 
    struct('part_type', 1, 'epsilon', 0.5, 'coupling', -0.1, 'coupling_index', 5, 'sigma', 0.1*1.j, 'fermi', 0)
};

N = 201;
voltages = linspace(-2.0,2.0, N)/2;
%voltages = [0.5];
%currents = iv_curve(voltages, H0, extends);
%N = 11;
%voltages = linspace(0,0.2,N);

kt = 0.1;

[currents, TT, DOS, energies, U, rho, fermiFactors] = sc_iv_curve(voltages, H0, extends, kt);
%currents = iv_curve(voltages, H0, extends);

currents = currents * .0066236; % To A

figure
plot(voltages, currents, 'r');
xlabel('Voltage (V)');
ylabel('Current (A)');
hold;
hline(0,'black');
vline(0,'black');

%print('epsilon_0_coupling_0d1_cp_index_2 _N201','-dpng');

[X, Y] = meshgrid(voltages,energies);

figure
s1 = surf(X, Y, TT,'EdgeColor','none');
xlabel('Voltage (V)');
ylabel('Energies (eV)');
zlabel('Transmission coefficient');

figure
s2 = surf(X, Y, DOS, 'EdgeColor','none');
xlabel('Voltage (V)');
ylabel('Energies (eV)');
zlabel('DOS');

figure
surf(X, Y, fermiFactors,'EdgeColor','none');
xlabel('Voltage (V)');
ylabel('Energies (eV)');
zlabel('Fermi distribuction');

figure
surf(X, Y, TT .* fermiFactors,'EdgeColor','none');
xlabel('Voltage (V)');
ylabel('Energies (eV)');
zlabel('Transmission coefficient * Fermi factor');