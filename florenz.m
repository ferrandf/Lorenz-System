function v = florenz(t,x)
global  Pr b r
L = [-Pr Pr 0 ; r -1 0 ; 0 0 -b] ;
v = L*x + [0 ; -x(1)*x(3) ; x(1)*x(2)] ;
