function [statedot,y]=segdynshell(t,state,u0,parms)
%function [statedot,y]=segdynshell(t,state,u0,parms)
%
% deze functie kan worden aangeroepen door een standaard integrator;
% gebruikmakend van segdyn (Casius et al., 2004) wordt de versnelling van 
% een keten van segmenten berekend
%
% inputs: 
% t (scalar, huidige tijd)
% state (vector, huidige waarde van de toestand van het systeem; zie paragraaf 5.5 van Casius et al.)
% u0 (vector, steady state waarden voor de ingangen van het systeem)
% parms (struct, bevat alle verder benodigde parameters, enige eis is dat
% de velden binnen deze struct in de verschillende functies consistent met
% elkaar zijn
%
% outputs:
% statedot (KOLOMvector, afgeleide van state en dus zelfde lengte als state)
% y (KOLOMvector, alle uitgangen van het systeem)
%   NB y wordt genegeerd tijdens het uitvoeren van een simulatie met een standaard MATLAB integrator
%   daarom is het ivm rekentijd verstandig y alleen te berekenen als dat
%   nodig is; dit is te regelen vanuit het hoofdprogramma via
%   parms.calculate_outputs; zie voorbeeld hoofdprogramma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: segdynstate (volledige vector coordinaten (zie segdyn)) maken uit state
% als er geen additionele dynamica is dan zullen segdynstate en state identiek zijn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = parms.segparms.L;
segdynstate = state;

fi = state(1:6);
fip = state(7:12);

%Berekenen van locaties van knie en heup

%locatie enkel, knie en heup
rE = [L(1)*cos(fi(1)) + L(2)*cos(fi(2)); L(1)*sin(fi(1)) + L(2)*sin(fi(2)) ];
rK = [L(1)*cos(fi(1)) + L(2)*cos(fi(2)) + L(3)*cos(fi(3)); L(1)*sin(fi(1)) + L(2)*sin(fi(2)) + L(3)*sin(fi(3))];
rH = [L(1)*cos(fi(1)) + L(2)*cos(fi(2)) + L(3)*cos(fi(3)) + L(4)*cos(fi(4)) ; L(1)*sin(fi(1)) + L(2)*sin(fi(2)) + L(3)*sin(fi(3)) + L(4)*sin(fi(4))];
rS = [L(1)*cos(fi(1)) + L(2)*cos(fi(2)) + L(3)*cos(fi(3)) + L(4)*cos(fi(4)) + L(5)*cos(fi(5))  ; L(1)*sin(fi(1)) + L(2)*sin(fi(2)) + L(3)*sin(fi(3)) + L(4)*sin(fi(4))+ L(5)*sin(fi(5))];

vE = [-L(1)*sin(fi(1))*fip(1) + -L(2)*sin(fi(2))*fip(2) ; L(1)*cos(fi(1))*fip(1) + L(2)*cos(fi(2))*fip(2)];
vK = [-L(1)*sin(fi(1))*fip(1) + -L(2)*sin(fi(2))*fip(2) + -L(3)*sin(fi(3))*fip(3); L(1)*cos(fi(1))*fip(1) + L(2)*cos(fi(2))*fip(2) + L(3)*cos(fi(3))*fip(3)];
vH = [-L(1)*sin(fi(1))*fip(1) + -L(2)*sin(fi(2))*fip(2) + -L(3)*sin(fi(3))*fip(3) + -L(4)*sin(fi(4))*fip(4); L(1)*cos(fi(1))*fip(1) + L(2)*cos(fi(2))*fip(2) + L(3)*cos(fi(3))*fip(3) +  L(4)*cos(fi(4))*fip(4)];
vS = [-L(1)*sin(fi(1))*fip(1) + -L(2)*sin(fi(2))*fip(2) + -L(3)*sin(fi(3))*fip(3) + -L(4)*sin(fi(4))*fip(4) + -L(5)*sin(fi(5))*fip(5); L(1)*cos(fi(1))*fip(1) + L(2)*cos(fi(2))*fip(2) + L(3)*cos(fi(3))*fip(3) +  L(4)*cos(fi(4))*fip(4) +  L(5)*cos(fi(5))*fip(5)];
%lokale assen onderbeen
ey_lok_ob = (rK-rE)./norm(rK-rE);
ex_lok_ob = [ey_lok_ob(2);-ey_lok_ob(1)];

%lokale assen bovenbeen
ey_lok_bb = (rH-rK)./norm(rH-rK);
ex_lok_bb = [ey_lok_bb(2);-ey_lok_bb(1)];

%lokale assen romp
ey_lok_romp = (rS-rH)./norm(rS-rH);
ex_lok_romp = [ey_lok_romp(2);-ey_lok_romp(1)];

%Rotatiematrices
R_ob = [ex_lok_ob ey_lok_ob];
R_bb = [ex_lok_bb ey_lok_bb];
R_romp = [ex_lok_romp ey_lok_romp];




%Locaties van origo en insertie in lokale assenstelsels definiëren (data
%van Delp)
% PELVIS:  The pelvic reference frame is fixed at the midpoint of the line
% connecting the two ASIS.
% FEMUR:  The femoral frame is fixed at the center of the femoral head
% TIBIA:   The tibial segment is located at the mid point of the line between
% the medial and lateral femoral epicondyles (see note belowİ).
% PATELLA: The patellar frame is located at the most distal point of the patella.
% CALCANeUS: The calcaneal frame is located at the most distal, inferior point
% on the posterior surface of the calcaneus. HIER MOETEN WE MISSCHIEN EEN
% PUNT TEN OPZICHTE VAN DE ENKEL VOOR VINDEN

% In the anatomical position, the X-axes point anteriorly, the Y-axes point
% superiorly, and the Z-axes point laterally.  Also note that this muscle
% file must be used with a joint file that has the same reference segments.


%Uit Delp nemen we ook standaard aan dat:
vmax = 10; %vezellengtes per seconde
%semimembranosus
semimem_or = [-0.1192; -0.1015]; %pelvis
semimem_in = [-0.0243 -0.0536]; %tibia
Fmax_semimem = 1030;

%semitendinosus, hier worden meerdere inserties gegeven, er is gekozen voor
%de middelste waarden omdat er drie waren
semiten_or = [-0.1237; -0.1043]; %pelvis
semiten_in = [-0.0113; -0.0193]; %tibia
Fmax_semiten = 328;

%biceps femoris long head
bicfemlh_or = [-0.1244; -0.1001]; %pelvis
bicfemlh_in = [-0.0081; -0.0729]; %tibia
Fmax_bicfemlh = 717;

%biceps femoris short head
bicfemsh_or = [0.0050; -0.2111]; %femur
bicfemsh_in = [-0.0101; -0.0725]; %tibia
Fmax_bicfemsh = 402;

%gluteus maximus superior middle en inferior, hier keuze gemaakt uit genoemde
%waarden
gmaxsup_or = [-0.1195; 0.0700]; %pelvis
gmaxsup_in = [-0.0457; -0.0248]; % femur
Fmax_gmaxsup = 382;

gmaxmid_or = [-0.1349; 0.0176]; % pelvis
gmaxmid_in = [-0.0426; -0.0530]; %femur
Fmax_gmaxmid = 546;

gmaxinf_or = [-0.1556; -0.0314];%pelvis
gmaxinf_in = [-0.0299; -0.1041];%femur
Fmax_gmaxinf = 368;

%iliacus en psoas deze hebben geen wrapping point uit file
ilia_or = [-0.0674 ; 0.0365]; %pelvis
ilia_ins = [0.0017 ;-0.0543]; %femur
ilia_wrap = 
Fmax_ilia = 429;

psoas_or = [-0.0647 ; 0.0887]; %pelvis
psoas_in = [0.0016 ;-0.0507]; %femur
psoas_wrap = 
Fmax_psoas = 371;

%quadriceps
%Hoe gaan we van patella aanhechting naar aanhechting op onderbeen???
%ten opzichte van wat is die wrapping gedefinieerd???
recfem_or = [-0.0295 ;-0.0311]; %pelvis
recfem_wrap = [0.0334 ;-0.4030]; %femur range KNEEang (-3.0, -1.46)/*radians*/
recfem_in = [0.0121 ; 0.0437]; %patella
Fmax_recfem = 779;

% vasmed_wrap1 = [0.0370; -0.4048];
% vasmed_wrap2 = [0.0274; -0.4255];
% werken voor kniehoeken (-3.0, -1.21) en (-3.0, -1.78)
vasmed_or = [0.0356 ; -0.2769]; %femur
vasmed_wrap = [0.0274 ; -0.4255]; %femur
vasmed_in = [0.0063 ; 0.0445];%patella
Fmax_vasmed = 1294;

vasint_or = [0.0290 ;-0.1924];%femur
vasint_wrap = [0.0343; -0.4030];%femur kniehoek (-3.0, -1.42)
vasint_in = [0.0058 ; 0.0480];%patella 
Fmax_vasint = 1365;

% vaslat_wrap1 = [0.0361; -0.4030];
% vaslat_wrap2 = [0.0253; -0.4243];
% kniehoek (-3.0, -1.21) en (-3.0, -1.92)
vaslat_or = [0.0048 ;-0.1854];%femur
vaslat_wrap = [0.0361; -0.4030];%femur
vaslat_in = [0.0103 ; 0.0423];%patella
Fmax_vaslat = 1871;

%gastrocnemius medialis
gasmed_or = [-0.0127 ;-0.3929]; %femur
gasmed_in = [0.0044 ; 0.0310]; %calcaneus
gasmed_wrap = [-0.0239; -0.4022]; % voor kniehoek (-0.77, 0.1)
Fmax_gasmed = 1113;

% gastocmenius lateralis
gaslat_or = [-0.0155 ;-0.3946]; %femur
gaslat_ins = [ 0.0044 ; 0.0310]; %calcaneus
gaslat_wrap = [-0.0254; -0.4018]; % voor kniehoek (-0.77, 0.1)
Fmax_gaslat = 488;


%DE PEZEN WORDEN VOOR NU ALS ONEINDIG STIJF EN MASSALOOS BESCHOUWD


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lengtes spieren
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Posities van origo en insertie in globaal assenstelsel
%biceps femoris short head
% bicfemsh_or = [0.0050; -0.2111]; %femur
% bicfemsh_in = [-0.0101; -0.0725]; %tibia
% Fmax_bicfemsh = 402;


glob_ori_bicfemsh = rH + R_bb*bicfemsh_or;
glob_ins_bicfemsh = rK + R_ob*bicfemsh_in;
l_bicfemsh = norm(glob_ins_bicfemsh-glob_ori_bicfemsh);

% Semimembranosus
glob_ori_semimem = rH + R_romp * semimem_or;
glob_ins_semimem = rK + R_ob * semimem_in;
l_semimem = norm(glob_ins_semimem - glob_ori_semimem);

% Semitendinosus
glob_ori_semiten = rH + R_romp * semiten_or;
glob_ins_semiten = rK + R_ob * semiten_in;
l_semiten = norm(glob_ins_semiten - glob_ori_semiten);

% Biceps femoris long head
glob_ori_bicfemlh = rH + R_romp * bicfemlh_or;
glob_ins_bicfemlh = rK + R_ob * bicfemlh_in;
l_bicfemlh = norm(glob_ins_bicfemlh - glob_ori_bicfemlh);

% Biceps femoris short head (monoarticulair)
glob_ori_bicfemsh = rH + R_bb * bicfemsh_or;
glob_ins_bicfemsh = rK + R_ob * bicfemsh_in;
l_bicfemsh = norm(glob_ins_bicfemsh - glob_ori_bicfemsh);

% glutes maximus superior
glob_ori_gmaxsup = rS + R_romp * gmaxsup_or;
glob_ins_gmaxsup = rH + R_bb * gmaxsup_in;
l_gmaxsup = norm(glob_ins_gmaxsup - glob_ori_gmaxsup);

% glutes maximus middle
glob_ori_gmaxmid = rS + R_romp * gmaxmid_or;
glob_ins_gmaxmid = rH + R_bb * gmaxmid_in;
l_gmaxmid = norm(glob_ins_gmaxmid - glob_ori_gmaxmid);

% glutes maximus inferior
glob_ori_gmaxinf = rS + R_romp * gmaxinf_or;
glob_ins_gmaxinf = rH + R_bb * gmaxinf_in;
l_gmaxinf = norm(glob_ins_gmaxinf - glob_ori_gmaxinf);

% Iliacus
glob_ori_ilia = rS + R_romp * ilia_or;
glob_wrap_ilia = rH + R_bb * ilia_wrap;
glob_ins_ilia = rK + R_bb * ilia_ins;
l_ilia = norm(glob_wrap_ilia - glob_ori_ilia) + norm(glob_ins_ilia - glob_wrap_ilia);

% Psoas
glob_ori_psoas = rS + R_romp * psoas_or;
glob_wrap_psoas = rH + R_bb * psoas_wrap;  
glob_ins_psoas = rK + R_bb * psoas_in;
l_psoas = norm(glob_wrap_psoas - glob_ori_psoas) + norm(glob_ins_psoas - glob_wrap_psoas);

% Rectus femoris
glob_ori_recfem = rS + R_romp * recfem_or;
glob_wrap_recfem = rK + R_bb * recfem_wrap;
glob_ins_recfem = rE + R_ob * recfem_in;
l_recfem = norm(glob_wrap_recfem - glob_ori_recfem) + norm(glob_ins_recfem - glob_wrap_recfem);

% Vastus medialis
glob_ori_vasmed = rH + R_bb * vasmed_or;
glob_wrap_vasmed = rK + R_bb * vasmed_wrap;
glob_ins_vasmed = rE + R_ob * vasmed_in;
l_vasmed = norm(glob_wrap_vasmed - glob_ori_vasmed) + norm(glob_ins_vasmed - glob_wrap_vasmed);

% Vastus intermedius
glob_ori_vasint = rH + R_bb * vasint_or;
glob_wrap_vasint = rK + R_bb * vasint_wrap;
glob_ins_vasint = rE + R_ob * vasint_in;
l_vasint = norm(glob_wrap_vasint - glob_ori_vasint) + norm(glob_ins_vasint - glob_wrap_vasint);

% Vastus lateralis
glob_ori_vaslat = rH + R_bb * vaslat_or;
glob_wrap_vaslat = rK + R_bb * vaslat_wrap;
glob_ins_vaslat = rE + R_ob * vaslat_in;
l_vaslat = norm(glob_wrap_vaslat - glob_ori_vaslat) + norm(glob_ins_vaslat - glob_wrap_vaslat);

% Gemiddelde kop (med + lat)
glob_ori_gas = rH + R_bb * gasmed_or;
glob_wrap_gas = rK + R_ob * [-0.0247; -0.4020]; % gecombineerde wrapping point
glob_ins_gas = rE + R_ob * gasmed_in;
l_gas = norm(glob_wrap_gas - glob_ori_gas) + norm(glob_ins_gas - glob_wrap_gas);

% gastocmenius medialis
glob_ori_gasmed = rH + R_bb * gasmed_or;
glob_wrap_gasmed = rK + R_bb * gasmed_wrap;
glob_ins_gasmed = rE + R_ob * gasmed_in;
l_gasmed = norm(glob_wrap_gasmed - glob_ori_gasmed) + norm(glob_ins_gasmed - glob_wrap_gasmed);

% gastocmenius lateralis lengte
glob_ori_gaslat = rH + R_bb * gaslat_or;
glob_wrap_gaslat = rK + R_bb * gaslat_wrap;
glob_ins_gaslat = rE + R_ob * gaslat_in;
l_gaslat = norm(glob_wrap_gaslat - glob_ori_gaslat) + norm(glob_ins_gaslat - glob_wrap_gaslat);





%snelheden van origo en insertie in globaal assenstelsel
v_ori_bicfemsh_glob = vH + ([ey_lok_bb -ex_lok_bb].*fip(4))*bicfemsh_or;
v_ins_bicfemsh_glob = vK + ([ey_lok_ob -ex_lok_ob].*fip(3))*bicfemsh_in;

%eenheidsvector van insertie t.o.v. oorsprong
%WAT GEBEURT HIER PRECIES???
e_inor_bicfemsh = (glob_ins_bicfemsh - glob_ori_bicfemsh)./(norm(glob_ins_bicfemsh - glob_ori_bicfemsh));
v_bicfemsh = dot(e_inor_bicfemsh,v_ori_bicfemsh_glob) - dot(e_inor_bicfemsh,v_ins_bicfemsh_glob);

%formules voor kracht komen uit syllabus van mechanische analyse

F0_bicfemsh = Fmax_bicfemsh*exp(-2*(fl_bicfemsh/opt_fl_bicfemsh-1)^2);
 

if v_bicfemsh > 0
    F1_bicfemsh = F0_bicfemsh*(0.3125/(v_bicfemsh/vmax_bicfemsh+0.25)-0.25);

else
    F1_bicfemsh = F0_bicfemsh;
end
%Momentsarm en moment berekenen
%vector van origo ten opzichte van het gewricht
r_ori_bicfemsh_knie = glob_ori_bicfemsh - rK;
pot_mom_bicfemsh = cross2d(r_ori_bicfemsh_knie,(F1_bicfemsh*e_inor_bicfemsh));
d_bicfemsh = cross2d(r_ori_bicfemsh_knie,e_inor_bicfemsh);

%wrapping points nodig voor: iliacus, psoas, de quadriceps en gastrocnemius






%lengte van monoarticulaire quadriceps berekenen door afstand van origo tot
%knie 










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: uitwendige krachten en momenten (voor zover niet onbekend) berekenen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fextx = [0 0 0 0 0 0];
Fexty = [-9.81*parms.segparms.m(1) -9.81*parms.segparms.m(2) -9.81*parms.segparms.m(3) -9.81*parms.segparms.m(4) -9.81*parms.segparms.m(5) -9.81*parms.segparms.m(6)];
Mext  = [0 0 0 0 0 0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: (optioneel) tijds- of toestands-afhankelijk deel van u berekenen
% bij constante input zal u gelijk zijn aan u0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = u0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: K en V aanmaken (voor beschrijving: zie segdyn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = ...
    [   1 0 0 0 0 0 1 ...% Frx (n+1)
        1 0 0 0 0 0 1 ...% Fry (n+1)
        1 1 1 1 1 1 1  ...% M (n+1)
        1 1 1 1 1 1  ...% Fextx (n)
        1 1 1 1 1 1  ...% Fexty (n)
        1 1 1 1 1 1  ...% Mext (n)     
        0 0 0 0 0 0 ...% phidd (n)
        0 0]; % basedd (2)

V =...
    [   0 NaN NaN NaN NaN NaN 0 ...% Frx (n+1)
        0 NaN NaN NaN NaN NaN 0 ...% Fry (n+1)
        0 0 0 0 0 0 0 ...% M (n+1)
        Fextx ...% Fextx (n)
        Fexty ...% Fexty (n)
        Mext ...% Mext (n)     
        NaN NaN NaN NaN NaN NaN...% phidd (n)
        NaN NaN]; % basedd (2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: (optioneel) Aconstraints en Bconstraints aanmaken
% (als er geen kinematische constraints zijn mogen Aconstraints en Bconstraints worden
% weggelaten in de aanroep van segdyn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Vastzetten trapas

Aconstraints(1,:) = [zeros(1,45) 1 0];
Bconstraints(1,:) = 0;

Aconstraints(2,:) = [zeros(1,46) 1 ];
Bconstraints(2,:) = 0;

K(1) = 0;
V(1) = NaN;

K(8) = 0;
V(8) = NaN;

%opleggen hoeksnelheid door de versnelling van fi1 gelijk aan 0 te maken
Aconstraints(3,:) = [zeros(1,39) 1 zeros(1,7)];
Bconstraints(3,1) = 0;

K(34) = 0;
V(34) = NaN;

%vastzetten enkelhoek phi2dp-phi3dp = 0

K(17) = 0;
V(17) = NaN;
Aconstraints(4,:) = [zeros(1,40) 1 -1 zeros(1,5)];
Bconstraints(4,1) = 0;

%vastzetten oriëntatie romp
Aconstraints(5,:) = [zeros(1,43) 1 zeros(1,3)];
Bconstraints(5,1) = 0;

K(38) = 0;
V(38) = NaN;



%vastzetten positie heup
Aconstraints(6,:) = zeros(1,47);
Aconstraints(7,:) = zeros(1,47);

Aconstraints(6,40:43) = [-parms.segparms.L(1)*sin(fi(1)) -parms.segparms.L(2)*sin(fi(2)) -parms.segparms.L(3)*sin(fi(3)) -parms.segparms.L(4)*sin(fi(4))]; 
Aconstraints(7,40:43) = [parms.segparms.L(1)*cos(fi(1)) parms.segparms.L(2)*cos(fi(2)) parms.segparms.L(3)*cos(fi(3)) parms.segparms.L(4)*cos(fi(4))]; 

Aconstraints(6,46) = 1;
Aconstraints(7,47) = 1;

Bconstraints(6,1) = parms.segparms.L(1)*cos(fi(1))*(fip(1)^2) + parms.segparms.L(2)*cos(fi(2))*(fip(2)^2) + parms.segparms.L(3)*cos(fi(3))*(fip(3)^2) + parms.segparms.L(4)*cos(fi(4))*(fip(4)^2);%-parms.segparms.L(5)*cos(fi(5))*(fip(5)^2) -parms.segparms.L(6)*cos(fi(6))*(fip(6)^2);
Bconstraints(7,1) = parms.segparms.L(1)*sin(fi(1))*(fip(1)^2) + parms.segparms.L(2)*sin(fi(2))*(fip(2)^2) + parms.segparms.L(3)*sin(fi(3))*(fip(3)^2) + parms.segparms.L(4)*sin(fi(4))*(fip(4)^2);%-parms.segparms.L(5)*sin(fi(5))*(fip(5)^2) -parms.segparms.L(6)*sin(fi(6))*(fip(6)^2);

K(26) = 0;
V(26) = NaN;
K(32) = 0;
V(32) = NaN;



% bij ingewikkelder constraints is het voor de leesbaarheid aan te bevelen om segparms en
% state eerst te ontrafelen in losse variabelen

% Aconstraints(1,:) = ...
%     [   ?? ...% Frx (n+1)
%         ?? ...% Fry (n+1)
%         ?? ...% M (n+1)
%         ?? ...% Fextx (n)
%         ?? ...% Fexty (n)
%         ?? ...% Mext (n)     
%         ?? ...% phidd (n)
%         ??]; % basedd (2)

%Bconstraints(1,1) = ??;

%K(??) = ??
%V(??) = ??

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% segdyn aanroepen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[segdynstatedot, Vnew, succes] = segdyn(segdynstate, parms.segparms, K, V,Aconstraints, Bconstraints);
% als er kinematische constraints zijn dan Aconstraints en Bconstraints toevoegen in bovenstaande aanroep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: statedot maken uit segdynstatedot
% in veel gevallen zullen segdynstatedot en statedot identiek zijn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

statedot = segdynstatedot;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: uitgangen y berekenen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y is een KOLOMvector van willekeurige lengte die alle interessant geachte afhankelijke variabelen
% bevat (bijv reactiekrachten, versnellingen, energietermen, zwaartepuntpositie etc etc)

if parms.calculate_outputs
    y = [Vnew, l_bicfemsh,v_bicfemsh,F1_bicfemsh,pot_mom_bicfemsh,d_bicfemsh];
end

