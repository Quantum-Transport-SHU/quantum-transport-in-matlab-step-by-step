function H=get_realspace_hamiltonian(R)
if all(R==[0,0])
    H = [0,1,1,0,0,0,0,0;
         1,0,0,1,0,0,0,0;
         1,0,0,0,1,0,0,0;
         0,1,0,0,0,1,0,0;
         0,0,1,0,0,1,1,0;
         0,0,0,1,1,0,0,1;
         0,0,0,0,1,0,0,0;
         0,0,0,0,0,1,0,0];
elseif all(R==[0,1])
    H = zeros(8,8);
    H(7,1) = 1;
    H(8,2) = 1;
elseif all(R==[0,-1])
    H = zeros(8,8);
    H(1,7) = 1; 
    H(2,8) = 1;
elseif all(R==[1,0])
    H =zeros(8,8);
    H(3,4) = 1;
    H(7,8) = 1;
elseif all(R==[-1,0])
    H = zeros(8,8);
    H(4,3) = 1;
    H(8,7) = 1;
else
    disp('unkown displacement')
end


    
     