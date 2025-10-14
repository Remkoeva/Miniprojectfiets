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
parms.calculate_outputs = 1;

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

%Initiële toestand
state0 = [pi*0.25 pi*0.75 pi*0.25 pi*0.75 pi*0.25 pi*-0.25 -2*pi 0 0 0 0 0 0 0 0 0];

mijnode = @(t,state) segdynshell(t,state,u0,parms);
[t,state]=ode113(mijnode,[0 1],state0,odeopt); 


 
%In het model is de hoeksnelheid van de crank constant, de positie van de
%trapas gefixeerd, de hoek van de enkel gefixeerd, de locatie van de heup
%gefixeerd, de hoek van de romp gefixeerd, de locatie van het uiteinde van
%de arm gefixeerd. de locatie van het uiteinde van de arm en de romphoek
%kunnen we mogelijk nog naar wens aanpassen.

