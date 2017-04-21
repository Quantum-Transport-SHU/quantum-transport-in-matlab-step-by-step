%This script is just to show how a hamiltonian matrix writes
%if there is serveral electrodes and several local states

h0 = 0;
t0 = -1.0;
H0 = [h0,t0,0,0,0; t0,h0,t0,0,0; 0,t0,h0,t0,0; 0,0,t0,h0,t0; 0,0,0,t0,h0];
%H0 describe a chain system, a molecule can be more complicate, but still can be described by a hamiltonian matrix.

h1 = 1;
t1 = -0.2;
T1 = [0,0,0,0,t1];
H1 = [H0, T1'; T1, h1];
%H1 describe a chain system coupled with a local state, where the onsite energy of the local state is h1, 
%and the local state has a relatively week coupling to the system (t1=-0.2) compared to the inner part of the chain.
%further more, the local state coupled to the fifth index of the H0 matrix, if the basis sequence in H0 is from
%left to right, the local state sits on the right.

h2 = 0.5;
t2 = -0.4;
T2 = [t2,0,0,0,0];
Sigma2 = 0.1*1.i;
H01 = H0;
H01(1,1) = H01(1,1)+Sigma2;
H2 = [h1,T2,0;T2',H01,T1';0,T2,h2];
%H2 describe a chain system coupled with a local state and an electrode whose selfenergy is Sigma2. the electrode has
%energy level (h2) and coupling (t2) to the chain, in this setting, the electrode has no interaction with previous local state.

h3 = 0.4;
t3 = -0.3;
t31 = -0.01;
T3 = [0,0,0,0,0,t3,t31];
Sigma3 = 0.05*1.j;
H02 = H2;
H02(6,6) = H02(6,6)+Sigma3;
H3 = [H2,T3';T3,h3]
%H3 describe a chain system coupled with a local state and two electrodes. The second electrode has very weak coupling to
%the local state.

%the case of more electrodes and more local states can be deduced from this example.


