%a simple file using Boryden mixer
classdef mixer < handle
    properties
        beta;
        nmaxold;
        step;
        R_iM;
        D_iM;
        eta_M;
        c_M;
        v_M;
        u_M;
        dN;
        verbose;
    end
    methods
        function self=mixer(beta,nmaxold,verbose)
            if nargin<3
                verbose=false;
            end
            if nargin<2
                nmaxold=6;
            end
            if nargin<1
                beta=0.1;
            end
            self.verbose = verbose;
            self.beta = beta;
            self.nmaxold = nmaxold;
            self.step=0;
        end
    
        function M = mix(self, M)
            if self.step>2
                self.D_iM(1) = [];
            end
            if self.step>0
                self.D_iM{end+1} = M - self.R_iM{end};
                fmin=sum(reshape(self.D_iM{end}.^2,1,[]));
                if self.verbose
                    fprintf('Mixer: broyden: fmin_G = %f\n',fmin);
                end
            end
            if self.step == 0
                self.eta_M = zeros(size(M));
            else
                if self.step >=2
                    self.c_M(:) = [];
                    if length(self.v_M) >= self.nmaxold
                        self.u_M(1) = [];
                        self.v_M(1) = [];
                    end    
                    temp_M = self.D_iM{2}-self.D_iM{1};
                    self.v_M{end+1}=temp_M/sum(reshape(temp_M.^2,1,[]));
                    if length(self.v_M) < self.nmaxold
                        nstep = self.step -1;
                    else
                        nstep = self.nmaxold;
                    end
                    for i=1:nstep
                        self.c_M{end+1} = sum(reshape(self.v_M{i}.*self.D_iM{2},1,[]));
                    end    
                    self.u_M{end+1} = self.beta*temp_M + self.R_iM{2}-self.R_iM{1};
                    usize = length(self.u_M);
                    for i=1:usize-1
                        a_M = sum(reshape(self.v_M{i}.*temp_M,1,[]));
                        self.u_M{usize} = self.u_M{usize}-a_M*self.u_M{i};
                    end
                end
                self.eta_M = self.beta*self.D_iM{end};
                usize = length(self.u_M);
                for i=1:usize
                    self.eta_M = self.eta_M - self.c_M{i}*self.u_M{i};
                end     
                M = M - self.D_iM{end} + self.eta_M;
                if self.step>=2
                    self.R_iM(1)=[];
                end
           end     
           self.R_iM{end+1} = M;
           self.step = self.step+1;
        end        
    end
end    

            
            


                                                    
                       
