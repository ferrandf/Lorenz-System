clear all
close all
format long e;
format compact;
global Pr b r Zref Tor 
Pr = 10 ; b = 8/3 ; 
options = odeset('reltol',1e-9,'abstol',1e-9) ;

%% 
%Q1: Per quin valor de r l'origen presenta una bifurcació? Quin tipus de
%bifurcació és? Està clar que x0 = (0,0,0) és un punt d'equilibri, ja que
%florenz(x0) = 0;


x0 = [0;0;0]';

%Estudiant l'origen tenim que és estable per r<1 i inestable per r>1.
%Els altres dos punts d'equilibri són:

qplus = [sqrt(b*(r-1)),sqrt(b*(r-1)),r-1];
qminus = [-sqrt(b*(r-1)),-sqrt(b*(r-1)),r-1];

%Amb una cerca dicotòmica podem establir a quin valor de r els valors
%propis passen de ser reals a complexos.

down = 1.2+1e-8;
up = 1.4;
T = 1;
it = 0;
while(T ==1)
    mid = (up+down)/2
    pol = [1  (Pr+b+1)  b*(Pr+mid)  2*Pr*b*(mid-1)];
    RR = roots(pol)
    if(imag(RR(2)) == 0)
        mid = mid+1e-5;
        pol = [1  (Pr+b+1)  b*(Pr+mid)  2*Pr*b*(mid-1)];
        RRR = roots(pol);
        if imag(RRR(2)) == 0
            down = mid;
            it = it+1
        else
            T=0;
            R = mid;
        end
    else
        up = mid;
    end
end

%Per tant, per Rc = R els valors propis passen de ser reals a tenir part
%imaginària diferent de zero.

%Ara, podem representar la bifurcació de Pitchfork de r = 1. El que fem és
%integrar a partir d'una certa condició inicial i veure a quin punt
%d'equilibri "convergeix".
X1 = []; X2 = []; X3 = []; X4 = []; X5 = []; X6 = [];
for k = 1e-1:5e-2 + 1e-6:2 %Fem un for per valors propers a r = 1;
    r = k;
    if r>1 %Si r és més gran que 1 canviem les condicions inicials.
        x0 = [sqrt(b*(r-1)),sqrt(b*(r-1)),r+5];
        x1 = [-sqrt(b*(r-1)),-sqrt(b*(r-1)),r+5];
    else
        x0 = [1,1,1];
        x1 = [1,1,1];
    end
    [T,Y] = ode45(@florenz,[0 100],x0,options);
    [U,Z] = ode45(@florenz,[0 100],x1,options);
    X1 = [X1 Y(end,1)]; X2 = [X2 Y(end,2)]; X3 = [X3 Y(end,3)];
    X4 = [X4 Z(end,1)]; X5 = [X5 Z(end,2)]; X6 = [X6 Z(end,3)];
end

R = 1e-1:5e-2 + 1e-6:2;

figure(11) %Es veu la bifurcació PITCHFORK a r = 1, en vermell q_+ i en blau q_-;
subplot(1,3,1)
plot(R,X1,'-k','linewidth',2,'Color','red');
hold on
plot(R,X4,'-k','linewidth',2,'Color','blue');
plot(linspace(0,1-1e-5),zeros(length(linspace(0,1-1e-5))),'-k','linewidth',2,'Color','black');
plot(linspace(1,2),zeros(length(linspace(1,2))),'--','linewidth',2,'Color','black');
xlabel('$r$'); ylabel('$x$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

subplot(1,3,2)
plot(R,X2,'-b','linewidth',2,'Color','red') ;
hold on
plot(R,X5,'-k','linewidth',2,'Color','blue');
plot(linspace(0,1-1e-5),zeros(length(linspace(0,1-1e-5))),'-k','linewidth',2,'Color','black');
plot(linspace(1,2),zeros(length(linspace(1,2))),'--','linewidth',2,'Color','black');

xlabel('$r$'); ylabel('$y$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

subplot(1,3,3)
plot(R,X3,'-r','linewidth',2,'Color','red') ;
hold on
plot(R,X6,'-k','linewidth',2,'Color','blue');
xlabel('$r$'); ylabel('$z$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

figure(12)
plot3(R,X1,X2,'-k')
grid on

%Definició de l'r crític:

rc = (Pr*(Pr+b+3))/(Pr-1-b);


%% Estudi de r crític







%% 

%NEWTON (Pts. fixos)-------------------------------------------------

set(0,'defaulttextinterpreter','latex')



% r = 19.5 ; xe = [-3.422720,0,19.5]'; %x d'estudi per r = 19.5
 r = 19.5 ; xe = [-4.7,0,19.2]'; %x d'estudi per r = 19.5
% r = 22 ;x0 = [3.7373397847413;2.1028363724307;20.00000892744];

% pla = florenz(0,xe)'
% w = null(pla); % Find two orthonormal vectors which are orthogonal to v
%    [P,Q] = meshgrid(-50:50); % Provide a gridwork (you choose the size)
%    X = xe(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
%    Y = xe(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
%    Z = xe(3)+w(3,1)*P+w(3,2)*Q;
%    surf(X,Y,Z,'EdgeColor', 'none', 'FaceColor', 'texturemap', 'FaceAlpha', 0.3)
 
% [X,FVAL] = fsolve(@(t,x)florenz(t,x),xe)










%% 

%DIBUIXOS -----------------------------------------------------------

set(0,'defaulttextinterpreter','latex')

global Pr b r Zref Tor 
Pr = 10 ; b = 8/3 ; 

options = odeset('reltol',1e-9,'abstol',1e-9) ;

% r = 24 ; x0 = [3;2;5]';
% r = 22 ;x0 = [3.7373397847413;2.1028363724307;20.00000892744];
% r = .5; x0 = [3;2;5]';
r = 19.5; x0 = [3;2;5]';
% r = 25; x0 = [3;2;5]';

%     r = 19; x0 = [4*sqrt(3),4*sqrt(3),24.4]; %Veiem que els punts C+,C-
%     són estables

x1 = [-3.422720,0,19.5];

[T,Y] = ode45(@florenz,[0 2.2],xe,options);

A = Y

figure(1)
subplot(3,1,1)
plot(T(:,1),Y(:,1),'-k','linewidth',2) ;
xlabel('$t$'); ylabel('$x$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;
%daspect([1 1 1])

subplot(3,1,2)
plot(T(:,1),Y(:,2),'-b','linewidth',2) ;
xlabel('$t$'); ylabel('$y$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

subplot(3,1,3)
plot(T(:,1),Y(:,3),'-r','linewidth',2) ;
xlabel('$t$'); ylabel('$z$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

figure(2)
plot3(Y(:,1),Y(:,2),Y(:,3),'-k','linewidth',2)


hold on

plot3(sqrt(b*(r-1)), sqrt(b*(r-1)), r-1, 'o', 'Color', 'red') %Punts C+,C-
plot3(-sqrt(b*(r-1)), -sqrt(b*(r-1)), r-1, 'o', 'Color', 'red')
plot3(x0(1),x0(2),x0(3),'o','Color','green')
plot3(xe(1),xe(2),xe(3),'o','Color','cyan')
hold on
pla = florenz(0,xe)'
w = null(pla); % Find two orthonormal vectors which are orthogonal to v
   [P,Q] = meshgrid(-50:50); % Provide a gridwork (you choose the size)
   X = xe(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
   Y = xe(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
   Z = xe(3)+w(3,1)*P+w(3,2)*Q;
   surf(X,Y,Z,'EdgeColor', 'none', 'FaceColor', 'texturemap', 'FaceAlpha', 0.3)

xlabel('$x$'); ylabel('$y$') ;  zlabel('$z$'), grid on  
a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;
hold off


%% 

%Prova per 19.5 ------------------------------:
r

zp = 6.93375;
xp = [-4.3;-4.3;zp]

pla = florenz(0,xp)';

[T,B] = ode45(@florenz,[0 100],xp,options);

B

figure(3)
subplot(3,1,1)
plot(T(:,1),B(:,1),'-k','linewidth',2) ;
xlabel('$t$'); ylabel('$x$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;
%daspect([1 1 1])

subplot(3,1,2)
plot(T(:,1),B(:,2),'-b','linewidth',2) ;
xlabel('$t$'); ylabel('$y$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

subplot(3,1,3)
plot(T(:,1),B(:,3),'-r','linewidth',2) ;
xlabel('$t$'); ylabel('$z$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

figure(4)
plot3(B(:,1),B(:,2),B(:,3),'-k','linewidth',2)


hold on

plot3(sqrt(b*(r-1)), sqrt(b*(r-1)), r-1, 'o', 'Color', 'red') %Punts C+,C-
plot3(-sqrt(b*(r-1)), -sqrt(b*(r-1)), r-1, 'o', 'Color', 'red')
plot3(xp(1),xp(2),xp(3),'o','Color','cyan')
hold on
pla = florenz(0,xp)'
w = null(pla); % Find two orthonormal vectors which are orthogonal to v
   [P,Q] = meshgrid(-50:50); % Provide a gridwork (you choose the size)
   X = xp(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
   Y = xp(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
   Z = xp(3)+w(3,1)*P+w(3,2)*Q;
   surf(X,Y,Z,'EdgeColor', 'none', 'FaceColor', 'texturemap', 'FaceAlpha', 0.3)

xlabel('$x$'); ylabel('$y$') ;  zlabel('$z$'), grid on  
a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;
hold off

%Pts interés 19.5:

e0 = [-2.536965104432087e+00    -4.160924541554467e+00     6.755942099858802e+00  9.518633916894746e-01];
e1 = [-2.566697491412065e+00    -4.212484871544123e+00     6.742639737411485e+00  9.536820512400535e-01];
e2 = [-2.596829229629287e+00    -4.264706130262444e+00     6.729872142397670e+00  9.555007107906323e-01];
e3 = [-2.627365095756077e+00    -4.317595300523674e+00     6.717648639956355e+00  9.573193703412111e-01];
e4 = [-2.658309907203616e+00    -4.371159356482967e+00     6.705978806266057e+00  9.591380298917899e-01];
e5 = [-2.689953881491605e+00    -4.425898734240063e+00     6.694774679608803e+00  9.609731287792871e-01];


t = [e0(4) e1(4) e2(4) e3(4) e4(4) e5(4)];
e_var_2 = [e0(2)+4.3 e1(2)+4.3 e2(2)+4.3 e3(2)+4.3 e4(2)+4.3 e5(2)+4.3];
e_var_1 = [e0(1) e1(1) e2(1) e3(1) e4(1) e5(1)];
e_var_3 = [e0(3) e1(3) e2(3) e3(3) e4(3) e5(3)];

p = polyfit(t,e_var_2,5)
t2_amb_complexes = roots(p)
figure(5)
plot(t,e_var_2,'o')
d1 = linspace(9.518633916894746e-01,9.609731287792871e-01,30);
d2 = polyval(p,d1);
hold on
plot(d1,d2)


t2 = 9.567168929072364e-01;

q1 = polyfit(t,e_var_1,5);
q3 = polyfit(t,e_var_3,5);

S1 = polyval(q1,t2)
S2 = -4.3
S3 = polyval(q3,t2)

S =[S1  S2  S3]'; 
figure(6)
plot3(B(:,1),B(:,2),B(:,3),'-k','linewidth',2)
hold on
plot3(sqrt(b*(r-1)), sqrt(b*(r-1)), r-1, 'o', 'Color', 'red') %Punts C+,C-
plot3(-sqrt(b*(r-1)), -sqrt(b*(r-1)), r-1, 'o', 'Color', 'red')
plot3(xp(1),xp(2),xp(3),'o','Color','cyan')
hold on
pla2 = florenz(0,S)'
w = null(pla2); % Find two orthonormal vectors which are orthogonal to v
[P,Q] = meshgrid(-50:50); % Provide a gridwork (you choose the size)
X = S(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
Y = S(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
Z = S(3)+w(3,1)*P+w(3,2)*Q;
surf(X,Y,Z,'EdgeColor', 'none', 'FaceColor', 'texturemap', 'FaceAlpha', 0.3)


[T,C] = ode45(@florenz,[0 1.2],S,options);


min = 100;
for k = 100:1:length(C)
    if min> norm(S-C(k,:)') && dot(C(k,:),pla2)-dot(S,pla2) < 0
        min = norm(S-C(k,:)');
        K = k;
    end
end

e0 = [C(K-2,1),C(K-2,2),C(K-2,3),T(K-2)];
e1 = [C(K-1,1),C(K-1,2),C(K-1,3),T(K-1)];
e2 = [C(K,1),C(K,2),C(K,3),T(K)];
e3 = [C(K+1,1),C(K+1,2),C(K+1,3),T(K+1)];
e4 = [C(K+2,1),C(K+2,2),C(K+2,3),T(K+2)];
e5 = [C(K+3,1),C(K+3,2),C(K+3,3),T(K+3)];

t = [e0(4) e1(4) e2(4) e3(4) e4(4) e5(4)];
e_var_2 = [e0(2) e1(2) e2(2) e3(2) e4(2) e5(2)];
e_var_1 = [e0(1) e1(1) e2(1) e3(1) e4(1) e5(1)];
e_var_3 = [e0(3) e1(3) e2(3) e3(3) e4(3) e5(3)];

p = polyfit(t,e_var_2,5)
t2_amb_complexes = roots(p)
figure(7)
plot(t,e_var_2,'o')
d1 = linspace(1.076,1.085,30);
d2 = polyval(p,d1);
hold on
plot(d1,d2)

% roots(p);

t2 = 6.503183494249812e-01;

q1 = polyfit(t,e_var_1,5);
q3 = polyfit(t,e_var_3,5);

SS1 = polyval(q1,t2)
SS2 = -2.914358560216681e+01
SS3 = polyval(q3,t2)

figure(8)

plot3(C(:,1),C(:,2),C(:,3),'-k','linewidth',2)
hold on
plot3(sqrt(b*(r-1)), sqrt(b*(r-1)), r-1, 'o', 'Color', 'red') %Punts C+,C-
plot3(-sqrt(b*(r-1)), -sqrt(b*(r-1)), r-1, 'o', 'Color', 'red')
plot3(S(1),S(2),S(3),'o','Color','cyan')
hold on
surf(X,Y,Z,'EdgeColor', 'none', 'FaceColor', 'texturemap', 'FaceAlpha', 0.3)
plot3(C(776,1),C(776,2),C(776,3),'o')
plot3(C(777,1),C(777,2),C(777,3),'o')



%% 

x6 = [3.7373397847413  2.1028363724307  20.00000892744];
r = 22;
[T,W] = ode45(@florenz,[0 10],x6,options);


figure(9)
subplot(3,1,1)
plot(T(:,1),W(:,1),'-k','linewidth',2) ;
xlabel('$t$'); ylabel('$x$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;
%daspect([1 1 1])

subplot(3,1,2)
plot(T(:,1),W(:,2),'-b','linewidth',2) ;
xlabel('$t$'); ylabel('$y$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

subplot(3,1,3)
plot(T(:,1),W(:,3),'-r','linewidth',2) ;
xlabel('$t$'); ylabel('$z$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

figure(10)
plot3(W(:,1),W(:,2),W(:,3),'-k','linewidth',2)

xlabel('$x$'); ylabel('$y$') ;  zlabel('$z$'), grid on  
a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;


%Última part: 

Q0= @(s1,s2) [S1-s1;S3-s2];

[X,FVAL] = fsolve(Q0,[xp(1) xp(3)])

X

