function [S,t2] = P(T,Y,x0,z)
%Funció de Poincaré -------------------------------------------
%
% Aquesta funció troba, donada la integració per ode45 des d'un punt x0
% guardada en [T,Y], el punt que talla el pla {Z = z} després de fer una
% "volta" sencera. 
% S és el punt en qüestió i t2 és el temps des de x0 fins que ha arribat a
% S (és el període de la "volta").
%

global Pr r b

K = 300;

while(Y(K,3)<z || norm(x0-Y(K,:)')>15) %Trobem quin punt està prop del que tallaria el pla
    K = K+1;
end

e0 = [Y(K-3,1),Y(K-3,2),Y(K-3,3),T(K-3)];
e1 = [Y(K-2,1),Y(K-2,2),Y(K-2,3),T(K-2)];
e2 = [Y(K-1,1),Y(K-1,2),Y(K-1,3),T(K-1)];
e3 = [Y(K,1),Y(K,2),Y(K,3),T(K)];
e4 = [Y(K+1,1),Y(K+1,2),Y(K+1,3),T(K+1)];
e5 = [Y(K+2,1),Y(K+2,2),Y(K+2,3),T(K+2)];


t = [e0(4) e1(4) e2(4) e3(4) e4(4) e5(4)];
e_var_2 = [e0(2) e1(2) e2(2) e3(2) e4(2) e5(2)];
e_var_1 = [e0(1) e1(1) e2(1) e3(1) e4(1) e5(1)];
e_var_3 = [e0(3)-z e1(3)-z e2(3)-z e3(3)-z e4(3)-z e5(3)-z];

p = polyfit(t,e_var_3,5);
t2_amb_complexes = roots(p);
figure(10) %Representa els punts amb els quals hi apliquem el polinomi interpolador i també en verd el valor del període de la "volta".
plot(t,e_var_3,'o')
t2 = 0;
for k =1:1:length(t2_amb_complexes) %Trobem el temps que ens interessa de totes les arrels de roots(p).
   if imag(t2_amb_complexes(k)) == 0 
       if real(t2_amb_complexes(k)) > t(3) && real(t2_amb_complexes(k)) < t(4)
           t2 = t2_amb_complexes(k);
           break
       end
   end
end
hold on
plot(t2,0,'o','Color','green')
xlabel('$t$')
ylabel('$z-r+1$')
% roots(p);



q1 = polyfit(t,e_var_1,5);
q2 = polyfit(t,e_var_2,5);

S1 = polyval(q1,t2);
S2 = polyval(q2,t2);
S3 = z;

S = [S1, S2, S3]


end