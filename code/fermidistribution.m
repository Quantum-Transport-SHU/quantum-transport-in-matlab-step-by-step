function y=fermidistribution(E, mu, kt)
x = (E-mu)/kt;
x = clip(x, -200,200);
y = 1./(1+exp(x));
end