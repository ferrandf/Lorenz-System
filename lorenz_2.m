clear all
close all
format long e;
format compact;
global Pr b r Zref Tor 
Pr = 10 ; b = 8/3 ; 



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











%DIBUIXOS -----------------------------------------------------------

set(0,'defaulttextinterpreter','latex')

global Pr b r Zref Tor 
Pr = 10 ; b = 8/3 ; 

options = odeset('reltol',1e-9,'abstol',1e-9) ;

% r = 24 ; x0 = [3;2;5]';
% r = 22 ;x0 = [3.7373397847413;2.1028363724307;20.00000892744];
% r = .5; x0 = [3;2;5]';
r = 18; x0 = [12;14;r-1]';
x1 = [3,-5,6]';
% r = 25; x0 = [3;2;5]';

%     r = 19; x0 = [4*sqrt(3),4*sqrt(3),24.4]; %Veiem que els punts C+,C-
%     són estables

% x0 = [1,1,15];
h = r-1;

[T,Y] = ode45(@florenz,[0 100],x0,options);


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
plot3(Y(729,1),Y(729,2),Y(729,3),'o','Color','green')
plot3(Y(730,1),Y(730,2),Y(730,3),'o','Color','magenta')
plot3(1.043029781533551e+01,1.454550975243910e+01,15,'o','Color','blue')
hold on
pla = [0 0 1];
w = null(pla); % Find two orthonormal vectors which are orthogonal to v
   [E1,E2] = meshgrid(-50:50); % Provide a gridwork (you choose the size)
   X = 0+w(1,1)*E1+w(1,2)*E2; % Compute the corresponding cartesian coordinates
   W = 0+w(2,1)*E1+w(2,2)*E2; %   using the two vectors in w
   Z = h+w(3,1)*E1+w(3,2)*E2;
   surf(X,W,Z,'EdgeColor', 'none', 'FaceColor', 'texturemap', 'FaceAlpha', 0.3)

xlabel('$x$'); ylabel('$y$') ;  zlabel('$z$'), grid on  
a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;
hold off


%Ara trobem el cícle límit:
S = P(T,Y,x0,h);
[X,FVAL] = fsolve(@(x)Q(x,h),x0);

[U,Z] = ode45(@florenz,[0 10],X,options);
[V,L] = ode45(@florenz,[0,100],x1,options);
%El període de l'òrbita ens el dóna P:
[D,peri] = P(U,Z,X,h);
peri

figure(4)
subplot(3,1,1)
plot(U(:,1),Z(:,1),'-k','linewidth',1) ;
xlabel('$t$'); ylabel('$x$','rotation',0) ; grid on  

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

figure(5)
plot3(Z(:,1),Z(:,2),Z(:,3),'-k','linewidth',2,'Color','green')
hold on
plot3(L(:,1),L(:,2),L(:,3),'-k','linewidth',2)

hold on

plot3(sqrt(b*(r-1)), sqrt(b*(r-1)), r-1, 'o', 'Color', 'red') %Punts C+,C-
plot3(-sqrt(b*(r-1)), -sqrt(b*(r-1)), r-1, 'o', 'Color', 'red')
plot3(X(1),X(2),X(3),'o','Color','green')
hold on
pla = [0 0 1];
w = null(pla); % Find two orthonormal vectors which are orthogonal to v
   [E1,E2] = meshgrid(-50:50); % Provide a gridwork (you choose the size)
   W1 = 0+w(1,1)*E1+w(1,2)*E2; % Compute the corresponding cartesian coordinates
   W2 = 0+w(2,1)*E1+w(2,2)*E2; %   using the two vectors in w
   W3 = r-1+w(3,1)*E1+w(3,2)*E2;
   surf(W1,W2,W3,'EdgeColor', 'none', 'FaceColor', 'texturemap', 'FaceAlpha', 0.3)

xlabel('$x$'); ylabel('$y$') ;  zlabel('$z$'), grid on  
a = gca;
a.TickLabelInterpreter = 'latex';
a.FontSize = 12;
hold off



% 
% 
% function Q = Q(x,z)
% [T,Y] = ode45(@florenz,[0 5],x,options);
% Q = P(T,Y,x,z)-x;
% return
% end