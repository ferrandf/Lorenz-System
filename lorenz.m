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
%Els altres dos punts d'equilibri són (per r>1):

qplus = [sqrt(b*(r-1)),sqrt(b*(r-1)),r-1];
qminus = [-sqrt(b*(r-1)),-sqrt(b*(r-1)),r-1];

%Amb una cerca dicotòmica podem establir a quin valor de r els valors
%propis passen de ser reals a complexos.

down = 1.2+1e-8;
up = 1.4;
T = 1;
it = 0;
while(T ==1)
    mid = (up+down)/2;
    pol = [1  (Pr+b+1)  b*(Pr+mid)  2*Pr*b*(mid-1)];
    RR = roots(pol);
    if(imag(RR(2)) == 0)
        mid = mid+1e-5;
        pol = [1  (Pr+b+1)  b*(Pr+mid)  2*Pr*b*(mid-1)];
        RRR = roots(pol);
        if imag(RRR(2)) == 0
            down = mid;
            it = it+1;
        else
            T=0;
            R = mid;
        end
    else
        up = mid;
    end
end


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
title('Bifurcació de Pitchfork')
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

%Definició de l'r crític:

rc = (Pr*(Pr+b+3))/(Pr-1-b);


%% Q4: r = 18
r = 18; x0 = [12;14;r-1]';
x1 = [5.5,10,7]';

h = r-1; %alçada del pla \Sigma = {z = h} que utilitzem per definir la funció de Poincaré.

[T,Y] = ode45(@florenz,[0 10],x0,options);


figure(1)

subplot(3,1,1)
plot(T(:,1),Y(:,1),'-k','linewidth',2) ;
xlabel('$t$'); ylabel('$x$','rotation',0) ; grid on  
title('Valors de òrbita per x0 = [12;14;r-1] amb r = 18')
a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

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
title('$r = 18$, $x0 = [12;14;r-1]$')

hold on

plot3(sqrt(b*(r-1)), sqrt(b*(r-1)), r-1, 'o', 'Color', 'red') %Punts C+,C-
plot3(-sqrt(b*(r-1)), -sqrt(b*(r-1)), r-1, 'o', 'Color', 'red')
plot3(x0(1),x0(2),x0(3),'diamond','Color','green')
hold on
pla = [0 0 1];
w = null(pla); 
   [E1,E2] = meshgrid(-50:50); 
   W1 = 0+w(1,1)*E1+w(1,2)*E2; 
   W2 = 0+w(2,1)*E1+w(2,2)*E2; 
   W3 = h+w(3,1)*E1+w(3,2)*E2;
   surf(W1,W2,W3,'EdgeColor', 'none', 'FaceColor', 'texturemap', 'FaceAlpha', 0.3)

xlabel('$x$'); ylabel('$y$') ;  zlabel('$z$'), grid on  
a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;
hold off


%Ara trobem el cicle límit:
[S,per] = P(T,Y,x0,h); %Ens acostem, d'un punt qualsevol x0 a un punt que té més sentit.

[X,FVAL] = fsolve(@(x)Q(x,h),S);

[U,Z] = ode45(@florenz,[0 10],X,options);
[V,L] = ode45(@florenz,[0,100],x1,options);
%El període de l'òrbita ens el dóna P:
[D,peri18] = P(U,Z,X,h);
%peri18 és el període del cicle límit per r = 18

figure(3) %Projecció al pla xy del cicle límit per r = 18
plot(Z(:,1),Z(:,2),'-k','linewidth',2,'Color',[0.30 0.56 0.645])
txt=('18');
text(1.5,1,txt,'Fontsize',16,'FontName','Times')
title('Òrbites periòdiques')

figure(4) %Representem els valors x(t),y(t),z(t) que s'assoleixen quan comencem en un punt de la òrbita periòdica
subplot(3,1,1)
plot(U(:,1),Z(:,1),'-k','linewidth',1) ;
xlabel('$t$'); ylabel('$x$','rotation',0) ; grid on  
title('Valors de òrbita periòdica per r = 18')

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

subplot(3,1,2)
plot(U(:,1),Z(:,2),'-b','linewidth',2) ;
xlabel('$t$'); ylabel('$y$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

subplot(3,1,3)
plot(U(:,1),Z(:,3),'-r','linewidth',2) ;
xlabel('$t$'); ylabel('$z$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

figure(5)%Representem el cicle límit en verd i una òrbita amb inici a x1 en negre.

plot3(Z(:,1),Z(:,2),Z(:,3),'-k','linewidth',2,'Color','green')
title('$r = 18$')
hold on
plot3(L(:,1),L(:,2),L(:,3),'-k','linewidth',2)
hold on

plot3(sqrt(b*(r-1)), sqrt(b*(r-1)), r-1, 'o', 'Color', 'red') %Punts C+,C-
plot3(-sqrt(b*(r-1)), -sqrt(b*(r-1)), r-1, 'o', 'Color', 'red')
plot3(X(1),X(2),X(3),'diamond','Color','green')
plot3(x1(1),x1(2),x1(3),'diamond','Color','red')

xlabel('$x$'); ylabel('$y$') ;  zlabel('$z$'), grid on  
a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;
hold off

%% Q4: r = 21;
r = 21; x0 = [12,12,r-1]';
x1 = [1,-5,6]';

h = r-1; %Serà l'alçada del pla \Sigma = {z = h}

[T,Y] = ode45(@florenz,[0 10],x0,options);

figure(6)

subplot(3,1,1)
plot(T(:,1),Y(:,1),'-k','linewidth',2) ;
xlabel('$t$'); ylabel('$x$','rotation',0) ; grid on  
title('Valors de òrbita per $x0 = [12;12;r-1]$ amb $r = 21$')
a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

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

figure(7) %S'hi representa una òrbita qualsevol del sistema de Lorenz, començant per x0.


plot3(Y(:,1),Y(:,2),Y(:,3),'-k','linewidth',2)

title('$r = 21$, $x0 = [12;12;r-1]$')
hold on

plot3(sqrt(b*(r-1)), sqrt(b*(r-1)), r-1, 'o', 'Color', 'red') %Punts C+,C-
plot3(-sqrt(b*(r-1)), -sqrt(b*(r-1)), r-1, 'o', 'Color', 'red')
plot3(x0(1),x0(2),x0(3),'diamond','Color','green') %Punt x0
hold on
pla = [0 0 1];
w = null(pla); % Dibuixem el pla
   [E1,E2] = meshgrid(-50:50); 
   W1 = 0+w(1,1)*E1+w(1,2)*E2; 
   W2 = 0+w(2,1)*E1+w(2,2)*E2; 
   W3 = h+w(3,1)*E1+w(3,2)*E2;
   surf(W1,W2,W3,'EdgeColor', 'none', 'FaceColor', 'texturemap', 'FaceAlpha', 0.2)

xlabel('$x$'); ylabel('$y$') ;  zlabel('$z$'), grid on  
a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;
hold off


%Ara trobem el cícle límit per r = 21:
[S,per] = P(T,Y,x0,h); %Ens acostem, d'un punt qualsevol x0 a un punt que té més sentit.

[X,FVAL] = fsolve(@(x)Q(x,h),S); 

[U,A] = ode45(@florenz,[0 10],X,options);
[V,L] = ode45(@florenz,[0,100],x1,options);
%El període de l'òrbita ens el dóna P:
[D,peri21] = P(U,A,X,h);
 %peri21 és el període del cicle límit en r = 21.

figure(3) %Recuperem la figura 3 on s'hi pot veure la projecció al pla xy dels cicles límit per r = 18 i r = 21.
grid on
hold on
plot(A(:,1),A(:,2),'-k','linewidth',2,'Color',[0.32 0.238 0.693])

xlabel('$x$')
ylabel('$y$')
txt=('21');
text(3,2.5,txt,'Fontsize',16,'FontName','Times')

figure(8) %Representem els valors x(t),y(t),z(t) que s'assoleixen quan comencem en un punt de la òrbita periòdica

subplot(3,1,1)
plot(U(:,1),A(:,1),'-k','linewidth',1) ;
xlabel('$t$'); ylabel('$x$','rotation',0) ; grid on  
title('Valors de òrbita periòdica per r = 21')

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

subplot(3,1,2)
plot(U(:,1),A(:,2),'-b','linewidth',2) ;
xlabel('$t$'); ylabel('$y$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

subplot(3,1,3)
plot(U(:,1),A(:,3),'-r','linewidth',2) ;
xlabel('$t$'); ylabel('$z$','rotation',0) ; grid on  

a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;

figure(9) %Representem el cicle límit en verd i una òrbita amb inici a x1 en negre.

plot3(A(:,1),A(:,2),A(:,3),'-k','linewidth',2,'Color','green')
hold on
plot3(L(:,1),L(:,2),L(:,3),'-k','linewidth',2)
title('$r = 21$')
hold on

plot3(sqrt(b*(r-1)), sqrt(b*(r-1)), r-1, 'o', 'Color', 'red') %Punts C+,C-
plot3(-sqrt(b*(r-1)), -sqrt(b*(r-1)), r-1, 'o', 'Color', 'red')
plot3(X(1),X(2),X(3),'diamond','Color','green')
plot3(x1(1),x1(2),x1(3),'diamond','Color','red')

xlabel('$x$'); ylabel('$y$') ;  zlabel('$z$'), grid on  
a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;
hold off


%% Estudi de r crític:

%Igual que en la bifurcació de Pitchfork, podem intentar representar la
%bifurcació de Hopf que hi ha en r = r_c.
%Aprofitarem que tenim les dades dels cicles límit per r = 18, r = 21 i r = 22 de l'enunciat per
%fer una pseudorepresentació de la bifurcació. Són els gràfics 12 i 13.

r = 22;
v = [3.7373397847413;2.1028363724307;20.00000892744];
[F,G] = ode45(@florenz,[0,3],v,options);

zmax = max(Z(:,1)); %18
zmin = min(Z(:,1));
amax = max(A(:,1)); %21
amin = min(A(:,1));
gmax = max(G(:,1)); %22
gmin = min(G(:,1));

u = [zmax,amax,gmax];
d = [zmin,amin,gmin];

X1 = [];
for k = 0.5:5e-1 + 1e-6:25 

    r = k
    
    x0 = [sqrt(b*(r-1)),sqrt(b*(r-1)),r+1];
    x1 = [-sqrt(b*(r-1)),-sqrt(b*(r-1)),r+1];

    [T,Y] = ode45(@florenz,[0 100],x0,options);
    X1 = [X1 Y(end,1)];
    if  k > rc-3e-1
        u = [u Y(end,1)];
        d = [d Y(end,1)];
    end
    
end

R = 0.5:5e-1 + 1e-6:25;

figure(12) %Veiem en r = 1 la bifurcació Pitchfork i a r = rc la de hopf.
plot(R,X1,'-')
hold on
Op = [18 21 22 rc-1e-3];
plot(Op,u,'--');
plot(Op,d,'--');
xlabel('r')
ylabel('x')

zmax = max(Z(:,2)); %18
zmin = min(Z(:,2));
amax = max(A(:,2)); %21
amin = min(A(:,2));
gmax = max(G(:,2)); %22
gmin = min(G(:,2));

u = [zmax,amax,gmax];
d = [zmin,amin,gmin];

X1 = [];
for k = 0.5:5e-1 + 1e-6:25

    r = k

    x0 = [sqrt(b*(r-1)),sqrt(b*(r-1)),r+1];
    x1 = [-sqrt(b*(r-1)),-sqrt(b*(r-1)),r+1];

    [T,Y] = ode45(@florenz,[0 100],x0,options);
    X1 = [X1 Y(end,2)];
    if  k > rc-3e-1
        u = [u Y(end,2)];
        d = [d Y(end,2)];
    end
    
end

R = 0.5:5e-1 + 1e-6:25;

figure(13)
plot(R,X1,'-')
hold on
Op = [18 21 22 rc];
plot(Op,u,'--');
plot(Op,d,'--');
xlabel('r')
ylabel('y')
