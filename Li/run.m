Rs = [[0,0];[0,1];[0,-1];[1,0];[-1,0]];
H = zeros(8, 8);
for i = 1:5
    R = Rs(i,:);
    h = get_realspace_hamiltonian(R);
    H = H + h;
end
eta = 1e-6;
G = inv(1.j*eta*eye(8)-H);
