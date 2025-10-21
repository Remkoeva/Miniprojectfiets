clear all; close all;
% model van fietser waarin gekeken wordt wat het effect is van het fietsen
% met/ zonder handen. 

% Segmenten: 1 = crank, 2 = voet (van contactpunt op trapper tot enkel),
% 3 = onderbeen, 4 = bovenbeen, 5 = romp (wat doen we met hoofd?)  en 6 =
% arm (kiezen we hier voor onderarm en bovenarm los en samen, dit bepaalt
% namelijk berekening voor traagheidsmomenten en c).
% Voor segmenten 3-6 wordt gebruik gemaakt van tabel uit syllabus
% mechanische analyse. Voor 1 en 2 wordt gebruik gemaakt van:
% Van Soest, Casius (2000), Hieruit wordt ook de hoogte van het zadel
% afgeleid. Zadelhoogte = 0.96 trochanteric height
% (0.96*(londerbeen+lbovenbeen???)

ltot = 1.8;
mtot = 80;
u0 = 0;

%d is gedeifinieerd als afstand van proximale einde tot massamiddelpunt van
%segment

parms.segparms.L = [0.17 0.17 0.25*ltot 0.24*ltot 0.30*ltot (0.17+0.16*ltot)];
parms.segparms.d = [(0.17-0.09) (0.17-0.12) 0.43*(0.25*ltot) 0.43*(0.24*ltot) 0.41*(0.30*ltot) 0.25];
parms.segparms.m = [0.20 1.23 0.05*mtot 0.11*mtot 0.45*mtot (0.03+0.02)*mtot ];
parms.segparms.J = [0.00 0.01 ((0.28*parms.segparms.L(3))^2)*parms.segparms.m(3) ((0.27*parms.segparms.L(4))^2)*parms.segparms.m(4) ((0.25*parms.segparms.L(5))^2)*parms.segparms.m(5) 0.12];
parms.segparms.j = [0.00 0.01 ((0.28*parms.segparms.L(3))^2)*parms.segparms.m(3) ((0.27*parms.segparms.L(4))^2)*parms.segparms.m(4) ((0.25*parms.segparms.L(5))^2)*parms.segparms.m(5) 0.12];
parms.calculate_outputs = 1;

L = parms.segparms.L;

parms.stick.dostick=1; % 1 for yes, 0 for no
parms.stick.axisvector=[-2 2 -2 2]; % display area
parms.stick.timestep=0.05; % animation timestep
parms.stick.realtime_if_possible=1; % if 1 then wait for real time, if 0 then plot when ready
parms.stick.fiindex=[1 2 3 4 5 6]; % vector of indices of segment angles in vector state that are to be plotted
parms.stick.baseindex= [13 14]; % vector containing index of x and y base
                           % position in vector state; supply NaN in case of fixed base

mijnodestick = @(t,state,flag) odestick(t,state,flag,u0,parms);
% parameters van de integrator instellen; zie helpfile bij odeopt


odeopt = odeset('abstol',1e-8,'reltol',1e-8,'outputfcn',mijnodestick);

%Initiële toestand de 15 zeros op het einde zijn de activaties van de
%spieren
%state0 = [pi*-0.5 pi*0.75 pi*0.25 pi*0.75 pi*0.3 -pi*0.5 -4*pi -3.356844886608543 -3.356844886608543 -2.175732796875910 0 0 0 0 0 0 zeros(1,15)];
%goede positie
%state0 = [pi*-0.5 pi*0.85 pi*0.35 pi*0.70 pi*0.3 -pi*0.5 -4*pi -3.878240711416198 -3.878240711416198 -0.806805175964899 0 0 0 0 0 0 zeros(1,15)];
state0 = [pi*-0.5 pi*0.85 pi*0.35 pi*0.70 pi*0.3 -pi*0.5 -4*pi -3.878240711416198 -3.878240711416198 -0.806805175964899 0 0 0 0 0 0 ];


% L(1)*sin(state0(1))*-4*pi =  -L(2)*sin(state0(2))*phip(2) -L(3)*sin(state0(3))*phip(3) -L(4)*sin(state0(4))*phip(4); %xpheup
% - L(1)*cos(state0(1))*-4*pi =  L(1)*cos(state0(1))*phip(2) + L(1)*cos(state0(1))*phip(3) + L(1)*cos(state0(1))*phip(4); %ypheup
% 0 = phip(2) - phip(3);

A = [-L(2)*sin(state0(2)) -L(3)*sin(state0(3)) -L(4)*sin(state0(4)); L(2)*cos(state0(2))  L(3)*cos(state0(3))  L(4)*cos(state0(4)); 1 -1 0];
B = [L(1)*sin(state0(1))*-4*pi; - L(1)*cos(state0(1))*-4*pi; 0];

x = pinv(A)*B;
 



% Definieer tijdsreeks (bijv. 0 tot 0.5 met stap 0.001)
% tspan = 0:0.001:0.5;
% 
% % Simulatie
% [t, state] = ode113(@(t,state) segdynshellminiproject(t, state, u0, parms),tspan, state0, odeopt);



mijnode = @(t,state) segdynshellminiproject(t,state,u0,parms);
[t,state]=ode113(mijnode,[0 0.5],state0,odeopt); 



%In het model is de locatie van de trapas gefixeerd, de hoeksnelheid van de crank constant, de positie van de
%trapas gefixeerd, de hoek van de enkel gefixeerd, de locatie van de heup
%gefixeerd, de hoek van de romp gefixeerd, de locatie van het uiteinde van
%de arm gefixeerd. de locatie van het uiteinde van de arm en de romphoek
%kunnen we mogelijk nog naar wens aanpassen.
parms.calculate_outputs = true;
for i=1:length(t)
 [statedot(:,i),y(:,i)]=segdynshellminiproject(t(i),state(i,:)',u0,parms);
end

statedot=statedot';
y=y';

fitot = state(:,1:6);
fiptot = state(:,7:12);
pos = state(:,13:14);
posp = state(:,15:16);

[epot,ekinx,ekiny,erot,etot]=energy(fitot,fiptot,pos,posp,parms.segparms);
% 
% plot(t,etot)
% figure
% plot(t,fitot)
% figure
% plot(t,y(:,34))
% figure
% 
% vcrankx = -L(1).*sin(state(:,1)).*state(:,7);
% vcranky = L(1).*cos(state(:,1)).*state(:,7);
% vcrank = sqrt(vcrankx.^2+vcranky.^2);
% 
% plot(t,y(:,34).*state(:,7))

% figure
% plot(t,y(:,38))
% figure
% plot(t,state(:,3),t,state(:,4))
% legend('onderbeen', 'bovenbeen')
% figure
% plot(t,(pi-state(:,3)+state(:,4)))
% figure
% plot(t,y(:,48),t,y(:,49))

% figure
% plot((360-(pi-state(:,3)+state(:,4))*((360/(2*pi)))),y(:,51))
% title('Momentsarm van biceps femoris short head als functie van de kniehoek')
% 
% figure
% plot(t,y(:,53))
% title('Het moment wat de biceps femoris short head levert over de tijd')
% figure
% plot(t,y(:,52))
% title('activatie biceps femoris short head')
% 
% figure 
% plot(t,y(:,34))
% plot(t,y(:,48),t,y(:,49),t,y(:,50));
% legend('l sub', 'l mid','l inf')
% figure
% plot(t,y(:,51),t,y(:,52),t,y(:,53));
% legend('v sub', 'v mid', 'v inf')
% figure
% plot(t,state(:,17),t,state(:,18),t,state(:,19),t,state(:,20),t,state(:,21),t,state(:,22),t,state(:,23),t,state(:,24),t,state(:,25),t,state(:,26),t,state(:,27),t,state(:,28),t,state(:,29),t,state(:,30),t,state(:,31))
% figure

 % plot(t,y(:,48),t,y(:,49),t,y(:,50),t,y(:,51));
 % legend('bicfemsh','semimem','semiten','bicfemlh')
 % 
 %   figure
 %   plot(t,y(:,52),t,y(:,53),t,y(:,54));
 %   legend('gmaxsup','gmaxmid','gmaxinf')
 %   figure
 %    plot(t,y(:,55),t,y(:,56),t,y(:,57));
 %   legend('ilia','psoas','recfem')
 %   figure
 %      plot(t,y(:,58),t,y(:,58),t,y(:,60));
 %   legend('vasmed','vasint','vaslat')
 %   figure
 % plot(t,y(:,61),t,y(:,62));
 %   legend('gasmed','gaslat')
 %   figure
 % plot(t,y(:,63),t,y(:,64));
 %   legend('fipheup','fipknie')

%  plot(t,y(:,7),t,y(:,14))
%  legend('F hand op stuur')
%  % plot(t,y(:,50))

% figure
% plot(state(:,1),y(:,48))
% title('crankhoek en vermogen bicfemsh')
figure

Momfr45 = (-cos(0.3*pi)*0.2214)*y(:,5) - (-sin(0.3*pi)*0.2214)*y(:,12);
Momfr56 = cos(0.3*pi)*(0.54-0.2214)*-y(:,6) - (sin(0.3*pi)*(0.54-0.2214))*-y(:,13);
Momtot = Momfr45+Momfr56;
plot(t,Momtot)
title('momtot rond heup')
figure
plot(t,y(:,65))

title('momcontrol')
figure
plot(t,y(:,66))

title('momheuptotspier')

figure
plot(t,y(:,76:79))
legend('P_bicfemlh',' P_semimem','  P_semiten')
figure
plot(t,y(:,79:84))
legend(' P_gmaxsup','P_gmaxmid','P_gmaxinf','P_ilia','P_psoas','P_recfem')
figure
plot(t,y(:,85))
title('v gmax')

