function y=fermidistribution(E, mu, kt)
if kt == 0
    if real(E)>mu
       y=0;
    else
       y=1;
    end
else
    y = 1./(1+exp((E-mu)/kt));
end
end