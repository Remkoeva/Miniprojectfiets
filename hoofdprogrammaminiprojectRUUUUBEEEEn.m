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


state0 = [pi*-0.5 pi*0.85 pi*0.35 pi*0.70 pi*0.40 -pi*0.5 -4*pi -3.878240711416198 -3.878240711416198 -0.806805175964899 0 0 0 0 0 0 ];

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


% 
% figure
% plot(t,y(:,91))
% title('P knie')
% figure
% plot(t,y(:,92))
% title('P Heup')


% figure
% plot(t,y(:,94:99))
% title('Spieractivaties rond heup met handen vast')
% ylim([0 1.2])
% legend('semimembranosus','semitendinosus','biceps femoris long head','gluteus maximus superior','gluteus maximus middle','gluteus maximus inferior')
% xlabel('tijd [s]')
% 
% figure
% plot(t,y(:,100:105))
% title('Potentiële vermogens heupstrekkers bij benodigde moment')
% ylabel('vermogen [W]')
% xlabel('tijd [s]')
% legend('semimembranosus','semitendinosus','biceps femoris long head','gluteus maximus superior','gluteus maximus middle','gluteus maximus inferior')


for i = 1: length(t)
    vcrank = [-L(1)*sin(state(i,1))*state(i,7); L(1)*cos(state(i,1))*state(i,7)];
    evcrank = vcrank/norm(vcrank);
    Fvoetcrank = [-y(i,2),-y(i,9)];
    Ftancrank(i) = dot(evcrank,Fvoetcrank);

    evrad = [-cos(state(1)); -sin(state(1))];
    Fradcrank(i) = dot(evrad,Fvoetcrank);
end

% load('Fradtanvast','Fradcrankvast','Ftancrankvast','tfcrank')
% 
% figure
% plot(t,Ftancrank)
% hold on
% plot(tfcrank,Ftancrankvast)
% legend('handen los','handen vast')
% title('tangentiële kracht op trapper')
% xlabel('tijd [s]');
% ylabel('kracht [N]')
% 
% 
% figure
% plot(t,Fradcrank)
% hold on
% plot(tfcrank,Fradcrankvast)
% legend('handen los','handen vast')
% title('radiale kracht op trapper')
% xlabel('tijd [s]');
% ylabel('kracht [N]')
% 
% load('Ptotvasthand','Ptotvast','tptotvast')
% 
% figure
% plot(t,y(:,105))
% hold on
% plot(tptotvast,Ptotvast)
% legend('handen los','handen vast')
% title('Totale spiervermogens (Fspier*verkortingssnelheid)')
% xlabel('tijd [s]');
% ylabel('spiervermogens [W]')
% load('gewrichtsvermogensvast','Pknievast','Pheupvast','tPgewricht')
% 
% figure
% plot(t,y(:,104))
% hold on
% plot(tPgewricht,Pknievast)
% legend('handen los','handen vast')
% title('Vermogens rond kniegewricht (Mgewricht*hoeksnelheid)')
% xlabel('tijd [s]');
% ylabel('gewrichtsvermogens [W]')
% 
% figure
% plot(t,y(:,105))
% hold on
% plot(tPgewricht,Pheupvast)
% legend('handen los','handen vast')
% title('Vermogens rond heupgewricht (Mgewricht*hoeksnelheid)')
% xlabel('tijd [s]');
% ylabel('gewrichtsvermogens [W]')

% figure
% plot(t,y(:,107:109))
% title('Vermogens van bi-articulaire hamstrings handen los')
% legend('semimembranosus','semitendinosus','biceps femoris long head')
% xlabel('tijd [s]');
% ylabel('spiervermogens [W]')
% 
% figure
% plot(t,y(:,110:112))
% title('Vermogens van gluteus maximus koppen handen los')
% legend('gluteus maximus superior','gluteus maximus middle','gluteus maximus inferior')
% xlabel('tijd [s]');
% ylabel('spiervermogens [W]')

%Vermogen op crank

Pcranklos = state(:,7).*-y(:,34);

load('Pcrankvast','Pcrankvast','tcrankvast')

figure
plot(t,Pcranklos)
hold on
plot(tcrankvast,Pcrankvast)
title('Vermogen op de crank (hoeksnelheid*-Mcrankvoet)')
xlabel('tijd [s]');
ylabel('crankvermogen [W]')
legend('handen los','handen vast')

Gemiddeldcrankvermogenhandenvast = mean(Pcrankvast);
GemiddeldcrankvermogenhandenLOS = mean(Pcranklos);

figure
plot(t,etot,t,epot)
figure
plot(t,y(:,113))
title('fip heup')
figure
plot(t,y(:,114))
title('fip knie')
figure
plot(t,y(:,115))
title('v semimem')
figure
plot(t,y(:,116))
title('v gmax')

