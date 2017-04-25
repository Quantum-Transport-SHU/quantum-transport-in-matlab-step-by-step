function y=fermidistribution(E, mu, kt)
if nargin < 3
    kt = 0;
end
x = (E-mu)/kt;
y = 1./(1+exp(x));
y(find(x>200))=0;
y(find(x<-200))=1;
end
