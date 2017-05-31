function G_lesser = get_negf(energy, S, H, Sigma, fermi, kt, voltage)

fl = fermidistribution(energy, fermi+voltage/2, kt);
fr = fermidistribution(energy, fermi-voltage/2, kt);
fc = fermidistribution(energy, fermi, kt);

sigmal = squeeze(Sigma(1,:,:));
sigmar = squeeze(Sigma(2,:,:));
sigma = sigmal + sigmar;

N = length(S);
NE = length(energy);

eta = 1e-5;

G_lesser = zeros(NE, N, N);
unit_mat = eye(N);
for i =1:NE
    G0 = inv((energy(i) + 1.j*eta) * S -H);
    G = inv((energy(i) + 1.j*eta) * S - H - sigmal - sigmar);
    sigmal_lesser = fl(i)*(sigmal' - sigmal);
    sigmar_lesser = fr(i)*(sigmar' - sigmar);
    GL0 = fc(i) * (G0 - G0');
    delta = (unit_mat +  G*sigma)* GL0 *(unit_mat + sigma'*G');
    G_lesser(i,:,:) = G*(sigmal_lesser + sigmar_lesser)*G' + delta;
end

end