function [S] = Q(x,z)

%Funció Q(x) = P(x)-x
% Donat un punt x i una alçada del pla {Z = z}, la funció retorna la
% diferència entre el punt x i la seva imatge per la funció de Poincaré
% P(x).
% S: La diferència entre P(x) i x.
%

options = odeset('reltol',1e-9,'abstol',1e-9) ;
[T,Y] = ode45(@florenz,[0 5],x,options);
[D,peri] = P(T,Y,x,z)
S = D-x
end