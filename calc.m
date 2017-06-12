# not a function file
g2=1
T=1

	div1 = (1-exp(-T))*(1-exp(-g2*T));
c1 =  sqrt(g2/div1);
c2 = -2/(1+g2);
c3= exp(-T*(1+g2)/2)-1;

tot = (c1*c2*c3)^2
