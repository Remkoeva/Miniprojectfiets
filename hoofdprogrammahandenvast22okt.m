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

%d is gedefinieerd als afstand van proximale einde tot massamiddelpunt van
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


state0 = [pi*-0.5 pi*0.85 pi*0.35 pi*0.70 pi*0.4 -pi*0.5 -4*pi -3.878240711416198 -3.878240711416198 -0.806805175964899 0 0 0 0 0 0 ];

%berekenen wat G op romp zou zijn, dit toevoegen als extra massa in arm
%aangezien we Fexty op romp gebruiken om de locatie van de romp op zijn plek te
%houden

Moment_Romp = parms.segparms.d(5)*cos(state0(5))*parms.segparms.m(5)*-9.81;
Moment_Arm = parms.segparms.L(5)*cos(state0(5))*parms.segparms.m(6)*-9.81;
Totale_Moment = Moment_Romp + Moment_Arm;
Totale_Kracht = Totale_Moment / (parms.segparms.L(5)*cos(state0(5)));
Extra_Massa = Totale_Kracht / -9.81;
parms.segparms.m(6) = Extra_Massa;
%Uitrekenen van initial state

A = [-L(2)*sin(state0(2)) -L(3)*sin(state0(3)) -L(4)*sin(state0(4)); L(2)*cos(state0(2))  L(3)*cos(state0(3))  L(4)*cos(state0(4)); 1 -1 0];
B = [L(1)*sin(state0(1))*-4*pi; - L(1)*cos(state0(1))*-4*pi; 0];

x = pinv(A)*B;
 



% Definieer tijdsreeks (bijv. 0 tot 0.5 met stap 0.001)
% tspan = 0:0.001:0.5;
% 
% % Simulatie
% [t, state] = ode113(@(t,state) segdynshellminiproject(t, state, u0, parms),tspan, state0, odeopt);



mijnode = @(t,state) segdynshellminiproject(t,state,u0,parms);
[t,state]=ode113(mijnode,[0 5],state0,odeopt); 



%In het model is de locatie van de trapas gefixeerd, de hoeksnelheid van de crank constant, de positie van de
%trapas gefixeerd, de hoek van de enkel gefixeerd, de locatie van de heup
%gefixeerd, de hoek van de romp gefixeerd, de locatie van het uiteinde van
%de arm gefixeerd. de locatie van het uiteinde van de arm en de romphoek
%kunnen we mogelijk nog naar wens aanpassen.
parms.calculate_outputs = true;
for i=1:length(t)
 [statedot(:,i),y(:,i)]=segdynshellhandenvast(t(i),state(i,:)',u0,parms);
end

statedot=statedot';
y=y';

fitot = state(:,1:6);
fiptot = state(:,7:12);
pos = state(:,13:14);
posp = state(:,15:16);

[epot,ekinx,ekiny,erot,etot]=energy(fitot,fiptot,pos,posp,parms.segparms);


figure
plot(t,y(:,80))
title('P knie')
figure
plot(t,y(:,81))
title('P Heup')
figure
plot(t,y(:,34).*state(7))

%kracht op crank
% figure
% plot(t,sqrt(y(:,2).^2 +y(:,9).^2))
% figure
% plot(t,-y(:,9))
% figure
% plot(t,y(:,97:106))
% 
% vxcrank = -parms.segparms.L(1).*sin(state(:,1)).*state(:,7);
% vycrank = parms.segparms.L(1).*cos(state(:,1)).*state(:,7);
% 
% Crankpower = vxcrank.*-y(:,2) + vycrank.*-y(:,9);
% 
% figure
% plot(t,Crankpower)
% figure
% plot(t,y(:,107:121))
% legend_names = {'bicfemsh', 'semimem', 'semiten', 'bicfemlh', ...
%                 'gmaxsup', 'gmaxmid', 'gmaxinf', 'ilia', 'psoas', ...
%                 'recfem', 'vasmed', 'vasint', 'vaslat', 'gasmed', 'gaslat'};
% 
% legend(legend_names)
% 
% figure
% plot(t,y(:,122))
% title('fipknie')
% figure
% plot(t,y(:,123))
% title('fipheup')
% 
% figure
% plot(t,y(:,124))
% title('v vasmed')
% figure
% plot(t,y(:,125))
% title('f vasmed')

figure
plot(t,y(:,126))
title('P tot')
figure
plot(t,y(:,80)+y(:,81))
title('p knie plus p heup')
figure
figure
plot(t,y(:,127:135))
legend_names = {'bicfemlh', 'semimem', 'semiten', 'gmaxsup', 'gmaxmid', 'gmaxinf', 'ilia', 'psoas', 'recfem'};

legend(legend_names)
figure
plot(t,y(:,136:150))
legend_names = {'bicfemsh', 'semimem', 'semiten', 'bicfemlh', ...
                'gmaxsup', 'gmaxmid', 'gmaxinf', 'ilia', 'psoas', ...
                'recfem', 'vasmed', 'vasint', 'vaslat', 'gasmed', 'gaslat'};

legend(legend_names)
title('verkortingssnelheden')

%HAMSTRINGS LIJKEN AAN HET BEGIN EEN NEGATIEF VERMOGEN ROND DE HEUP TE
%HEBBEN, MISSCHIEN OMDAT ZE AAN STAAN OMDAT ZE ROND DE KNIE VERKORTEN, MAAR
%ROND DE HEUP VERLENGEN
