function y=fermidistribution(E, mu, kt)
if nargin < 3
   kt = 0;
end
x = (E-mu)/kt;
y = 1./(1+exp(x));
y(find(x>100))=0;
y(find(x<-100))=1;
y(find((E-mu)==0))=0.5;
end
