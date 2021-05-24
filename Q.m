function [S] = Q(x,z)
options = odeset('reltol',1e-9,'abstol',1e-9) ;
[T,Y] = ode45(@florenz,[0 5],x,options);
[D,peri] = P(T,Y,x,z);
S = D-x;
end