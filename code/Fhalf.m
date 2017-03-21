function y=Fhalf(x) 
    xx=linspace(0,abs(x)+10,251);dx=xx(2)-xx(1); 
    fx=(2*dx/sqrt(pi))*sqrt(xx)./(1+exp(xx-x)); 
    y=sum(fx); 
end 