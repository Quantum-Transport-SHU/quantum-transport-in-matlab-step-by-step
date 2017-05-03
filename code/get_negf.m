function G_lesser = get_negf(energy, S, H, Sigma, fermi, kt, voltage)

fl = fermidistribution(energy, fermi+voltage/2, kt);
fr = fermidistribution(energy, fermi-voltage/2, kt);

sigmal = squeeze(Sigma(1,:,:));
sigmar = squeeze(Sigma(2,:,:));

N = length(S);
NE = length(energy);

G_lesser = zeros(NE, N, N);

for i =1:NE
    G = inv(energy(i) * S - H - sigmal - sigmar);
    sigmal_lesser = 1.j*fl(i)*(sigmal' - sigmal);
    sigmar_lesser = 1.j*fr(i)*(sigmar' - sigmar);
    G_lesser(i,:,:) = G*(sigmal_lesser + sigmar_lesser)*G';
end