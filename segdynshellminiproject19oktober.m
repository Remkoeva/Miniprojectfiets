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
segdynstate = state(1:16);
a_bicfemsh = state(17);
a_semimem  = state(18);
a_semiten  = state(19);
a_bicfemlh = state(20);
a_gmaxsup  = state(21);
a_gmaxmid  = state(22);
a_gmaxinf  = state(23);
a_ilia     = state(24);
a_psoas    = state(25);
a_recfem   = state(26);
a_vasmed   = state(27);
a_vasint   = state(28);
a_vaslat   = state(29);
a_gasmed   = state(30);
a_gaslat   = state(31);



fi = state(1:6);
fip = state(7:12);

%Berekenen van locaties van knie en heup

%locatie enkel, knie en heup
rE = [L(1)*cos(fi(1)) + L(2)*cos(fi(2)); L(1)*sin(fi(1)) + L(2)*sin(fi(2)) ];
rK = [L(1)*cos(fi(1)) + L(2)*cos(fi(2)) + L(3)*cos(fi(3)); L(1)*sin(fi(1)) + L(2)*sin(fi(2)) + L(3)*sin(fi(3))];
rH = [L(1)*cos(fi(1)) + L(2)*cos(fi(2)) + L(3)*cos(fi(3)) + L(4)*cos(fi(4)) ; L(1)*sin(fi(1)) + L(2)*sin(fi(2)) + L(3)*sin(fi(3)) + L(4)*sin(fi(4))];
rS = [L(1)*cos(fi(1)) + L(2)*cos(fi(2)) + L(3)*cos(fi(3)) + L(4)*cos(fi(4)) + L(5)*cos(fi(5))  ; L(1)*sin(fi(1)) + L(2)*sin(fi(2)) + L(3)*sin(fi(3)) + L(4)*sin(fi(4))+ L(5)*sin(fi(5))];




% De locatie van de patella wordt t.o.v. het lokale assenstelsel van de
% tibia wordt gedefinieerd; dit komt doordat wij de pezen (en ligamenten)
% momenteel als onvervormbaar beschouwen, en de patellar tendon dus
% hiervoor gebruikt kan worden (zie verslag). 
% Gemiddelde lengte van patella = 34.3 mm
Pat_Hoek = 0.75*(state(3)-state(4));
Pat_ten_angle = deg2rad(90-12);
Pat_ten_length = 16.54 + 0.17*180;
Knie_Tuberosity = [0.039; -0.082];
Tuberosity_Patella = Pat_ten_length*[cos(Pat_ten_angle);sin(Pat_ten_angle)];
Knie_Patella = (Knie_Tuberosity + Tuberosity_Patella)/1000;

Knie_Patella_Distaal = rK + [Knie_Patella(1); Knie_Patella(2)];
rPat = Knie_Patella_Distaal;
Knie_Patella_Proximaal = Knie_Patella_Distaal + 0.0343 * [cos(Pat_Hoek);sin(Pat_Hoek)];

%de hoek van romp - hoek van asis tov heup in assenstelsel van romp =
%0.6523 rad
%de locatie van de oorsprong van het pelvis assenstelsel is de locatie van
%de pelvis + cos/sin(hoek romp+ 0.6523) * lengte van afstand asis tot heup (0.0906) 
firomppel = fi(5) - 0.6523;
Lpelromp = 0.0906;
rPel = [L(1)*cos(fi(1)) + L(2)*cos(fi(2)) + L(3)*cos(fi(3)) + L(4)*cos(fi(4)) + Lpelromp*cos(firomppel); L(1)*sin(fi(1)) + L(2)*sin(fi(2)) + L(3)*sin(fi(3)) + L(4)*sin(fi(4)) + Lpelromp*sin(firomppel)];

%voor de dingen ten opzichte van dit punt gebruiken we gewoon R_romp
%aangezien de oriëntatie van de romp en dit assenstelse vanuit de asis
%hetzelfde is

vE = [-L(1)*sin(fi(1))*fip(1) + -L(2)*sin(fi(2))*fip(2) ; L(1)*cos(fi(1))*fip(1) + L(2)*cos(fi(2))*fip(2)];
vK = [-L(1)*sin(fi(1))*fip(1) + -L(2)*sin(fi(2))*fip(2) + -L(3)*sin(fi(3))*fip(3); L(1)*cos(fi(1))*fip(1) + L(2)*cos(fi(2))*fip(2) + L(3)*cos(fi(3))*fip(3)];
vH = [-L(1)*sin(fi(1))*fip(1) + -L(2)*sin(fi(2))*fip(2) + -L(3)*sin(fi(3))*fip(3) + -L(4)*sin(fi(4))*fip(4); L(1)*cos(fi(1))*fip(1) + L(2)*cos(fi(2))*fip(2) + L(3)*cos(fi(3))*fip(3) +  L(4)*cos(fi(4))*fip(4)];
vS = [-L(1)*sin(fi(1))*fip(1) + -L(2)*sin(fi(2))*fip(2) + -L(3)*sin(fi(3))*fip(3) + -L(4)*sin(fi(4))*fip(4) + -L(5)*sin(fi(5))*fip(5); L(1)*cos(fi(1))*fip(1) + L(2)*cos(fi(2))*fip(2) + L(3)*cos(fi(3))*fip(3) +  L(4)*cos(fi(4))*fip(4) +  L(5)*cos(fi(5))*fip(5)];
vPel = vH + [-Lpelromp*sin(firomppel)*fip(5);  Lpelromp*cos(firomppel)*fip(5)];
vPat = vK ; %DIT MOET NOG AANGEPAST WORDENNNN

%lokale assen onderbeen
ey_lok_ob = (rK-rE)./norm(rK-rE);
ex_lok_ob = [ey_lok_ob(2);-ey_lok_ob(1)];

%lokale assen bovenbeen
ey_lok_bb = (rH-rK)./norm(rH-rK);
ex_lok_bb = [ey_lok_bb(2);-ey_lok_bb(1)];

%lokale assen romp
ey_lok_romp = (rS-rH)./norm(rS-rH);
ex_lok_romp = [ey_lok_romp(2);-ey_lok_romp(1)];

%lokale assen patella
ey_lok_pat = (Knie_Patella_Proximaal-Knie_Patella_Distaal)./norm(Knie_Patella_Proximaal-Knie_Patella_Distaal);
ex_lok_pat = [ey_lok_pat(2);-ey_lok_pat(1)];



%Rotatiematrices
R_ob = [ex_lok_ob ey_lok_ob];
R_bb = [ex_lok_bb ey_lok_bb];
R_romp = [ex_lok_romp ey_lok_romp];
R_pat = [ex_lok_pat ey_lok_pat];



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
semimem_in = [-0.0243; -0.0536]; %tibia
Fmax_semimem = 1030;
pees_sl_semimem = 0.3590;
opt_fl_semimem = 0.0800;
vmax_semimem = vmax*opt_fl_semimem;

%semitendinosus, hier worden meerdere inserties gegeven, er is gekozen voor
%de middelste waarden omdat er drie waren
semiten_or = [-0.1237; -0.1043]; %pelvis
semiten_in = [-0.0113; -0.0193]; %tibia
Fmax_semiten = 328;
pees_sl_semiten = 0.2620;
opt_fl_semiten = 0.2010;
vmax_semiten = vmax*opt_fl_semiten;

%biceps femoris long head
bicfemlh_or = [-0.1244; -0.1001]; %pelvis
bicfemlh_in = [-0.0081; -0.0729]; %tibia
Fmax_bicfemlh = 717;
pees_sl_bicfemlh = 0.3410;
opt_fl_bicfemlh = 0.1090;
vmax_bicfemlh = vmax*opt_fl_bicfemlh;

%biceps femoris short head
bicfemsh_or = [0.0050; -0.2111]; %femur
bicfemsh_in = [-0.0101; -0.0725]; %tibia
Fmax_bicfemsh = 402;
pees_sl_bicfemsh = 0.100;
opt_fl_bicfemsh = 0.1730;
vmax_bicfemsh = vmax*opt_fl_bicfemsh;

%gluteus maximus superior middle en inferior, hier keuze gemaakt uit genoemde
%waarden
gmaxsup_or = [-0.1195; 0.0700]; %pelvis
gmaxsup_in = [-0.0457; -0.0248]; %femur
Fmax_gmaxsup = 382;
pees_sl_gmaxsup = 0.1250;
opt_fl_gmaxsup = 0.1420;
vmax_gmaxsup = vmax*opt_fl_gmaxsup;

gmaxmid_or = [-0.1349; 0.0176];%pelvis
gmaxmid_in = [-0.0426; -0.0530];%femur
Fmax_gmaxmid = 546;
pees_sl_gmaxmid = 0.1270;
opt_fl_gmaxmid = 0.1470;
vmax_gmaxmid = vmax*opt_fl_gmaxmid;

gmaxinf_or = [-0.1556; -0.0314];%pelvis
gmaxinf_in = [-0.0299; -0.1041];%femur
Fmax_gmaxinf = 368;
pees_sl_gmaxinf = 0.1450;
opt_fl_gmaxinf = 0.1440;
vmax_gmaxinf = vmax*opt_fl_gmaxinf;


% gmax_or = [-0.1349; 0.0176];%pelvis
% gmax_in = [-0.0426; -0.0248];%femur
% Fmax_gmax = 382+546+368;


%iliacus en psoas
ilia_or = [-0.0674 ; 0.0365]; %pelvis
ilia_ins = [0.0017 ;-0.0543]; %femur
Fmax_ilia = 429;
pees_sl_ilia = 0.0900;
opt_fl_ilia = 0.1000;
vmax_ilia = vmax*opt_fl_ilia;



psoas_or = [-0.0647 ; 0.0887]; %pelvis
psoas_in = [0.0016 ;-0.0507]; %femur
Fmax_psoas = 371;
pees_sl_psoas = 0.1300;
opt_fl_psoas =0.1040;
vmax_psoas = vmax*opt_fl_psoas;

% ilpso_or = [-0.0647;0.0626 ]; %pelvis
% ilpso_in = [0.0016 ;-0.0507]; %femur
% Fmax_ilpso = 371;

%quadriceps
%Hoe gaan we van patella aanhechting naar aanhechting op onderbeen???
%ten opzichte van wat is die wrapping gedefinieerd???
recfem_or = [-0.0295 ;-0.0311]; %pelvis
recfem_wrap = [0.0334 ;-0.4030]; %femur range KNEEang (-3.0, -1.46)/*radians*/
recfem_in = [0.0121 ; 0.0437]; %patella
Fmax_recfem = 779;
pees_sl_recfem = 0.3460;
opt_fl_recfem = 0.0840;
vmax_recfem = vmax*opt_fl_recfem;

% vastus medialis
vasmed_or   = [0.0356 ; -0.2769];    % origine (femur)
vasmed_in   = [0.0063 ;  0.0445];    % insertie (patella)
vasmed_wrap1 = [0.0370; -0.4048];    % femur bij KNEEang (-3.0, -1.21)
vasmed_wrap2 = [0.0274; -0.4255];    % femur bij KNEEang (-3.0, -1.78)
Fmax_vasmed = 1294;
pees_sl_vasmed = 0.1260;
opt_fl_vasmed = 0.0890;
vmax_vasmed = vmax*opt_fl_vasmed;


% vastus int
vasint_or = [0.0290 ;-0.1924];       %femur
vasint_wrap = [0.0343; -0.4030];     %femur kniehoek (-3.0, -1.42)
vasint_in = [0.0058 ; 0.0480];       %patella 
Fmax_vasint = 1365;
pees_sl_vasint = 0.1360;
opt_fl_vasint = 0.0870;
vmax_vasint = vmax*opt_fl_vasint;

vaslat_or   = [0.0048 ;-0.1854];     % femur
vaslat_wrap1 = [0.0361; -0.4030];    % femur bij KNEEang(-3.0, -1.21)
vaslat_wrap2 = [0.0253; -0.4243];    % femur bij KNEEang(-3.0, -1.92)
vaslat_in   = [0.0103 ; 0.0423];     % patella
Fmax_vaslat = 1871;
pees_sl_vaslat = 0.1570;
opt_fl_vaslat = 0.0840;
vmax_vaslat = vmax*opt_fl_vaslat;


%gastrocnemius medialis
gasmed_or = [-0.0127 ;-0.3929]; %femur
gasmed_in = [0.0044 ; 0.0310]; %calcaneus
gasmed_wrap = [-0.0239; -0.4022]; % voor kniehoek (-0.77, 0.1)
Fmax_gasmed = 1113;
pees_sl_gasmed = 0.4080;
opt_fl_gasmed = 0.0450;
vmax_gasmed = vmax*opt_fl_gasmed;

% gastocmenius lateralis
gaslat_or = [-0.0155 ;-0.3946]; %femur
gaslat_ins = [ 0.0044 ; 0.0310]; %calcaneus
gaslat_wrap = [-0.0254; -0.4018]; % voor kniehoek (-0.77, 0.1)
Fmax_gaslat = 488;
pees_sl_gaslat = 0.3850;
opt_fl_gaslat = 0.0640;
vmax_gaslat = vmax*opt_fl_gaslat;  


% gas_or = [-0.0127 ;-0.3929];%femur
% gas_in = [0.0044 ; 0.0310];%calcaneus
% Fmax_gas = 1113+488;

%DE PEZEN WORDEN VOOR NU ALS ONEINDIG STIJF EN MASSALOOS BESCHOUWD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lengtes spieren
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Biceps femoris short head (femur → tibia) 
glob_ori_bicfemsh = rH + R_bb*bicfemsh_or;
glob_ins_bicfemsh = rK + R_ob*bicfemsh_in;
l_bicfemsh = norm(glob_ins_bicfemsh-glob_ori_bicfemsh);

% Semimembranosus (pelvis → tibia)
glob_ori_semimem = rPel + R_romp * semimem_or;
glob_ins_semimem = rK + R_ob * semimem_in;
l_semimem = norm(glob_ins_semimem - glob_ori_semimem);

% Semitendinosus (pelvis → tibia)
glob_ori_semiten = rPel + R_romp * semiten_or;
glob_ins_semiten = rK + R_ob * semiten_in;
l_semiten = norm(glob_ins_semiten - glob_ori_semiten);

% Biceps femoris long head (pelvis → tibia)
glob_ori_bicfemlh = rPel + R_romp * bicfemlh_or;
glob_ins_bicfemlh = rK + R_ob * bicfemlh_in;
l_bicfemlh = norm(glob_ins_bicfemlh - glob_ori_bicfemlh);

% Gluteus maximus superior (pelvis → femur)
glob_ori_gmaxsup = rPel + R_romp * gmaxsup_or;
glob_ins_gmaxsup = rH + R_bb * gmaxsup_in;
l_gmaxsup = norm(glob_ins_gmaxsup - glob_ori_gmaxsup);

% Gluteus maximus middle
glob_ori_gmaxmid = rPel + R_romp * gmaxmid_or;
glob_ins_gmaxmid = rH + R_bb * gmaxmid_in;
l_gmaxmid = norm(glob_ins_gmaxmid - glob_ori_gmaxmid);

% Gluteus maximus inferior
glob_ori_gmaxinf = rPel + R_romp * gmaxinf_or;
glob_ins_gmaxinf = rH + R_bb * gmaxinf_in;
l_gmaxinf = norm(glob_ins_gmaxinf - glob_ori_gmaxinf);

% Iliacus (pelvis → femur, via wrappunt)
glob_ori_ilia = rPel + R_romp * ilia_or;
glob_ins_ilia = rH + R_bb * ilia_ins;
l_ilia = norm(glob_ins_ilia - glob_ori_ilia);

% Psoas (pelvis → femur, via wrappunt)
glob_ori_psoas = rPel + R_romp * psoas_or;
glob_ins_psoas = rH + R_bb * psoas_in;
l_psoas = norm(glob_ins_psoas - glob_ori_psoas);

% Quadriceps: Rectus femoris (pelvis → patella, via femur wrap)
glob_ori_recfem = rPel + R_romp * recfem_or;
glob_wrap_recfem = rH + R_bb * recfem_wrap;
glob_ins_recfem = rPat + R_pat * recfem_in;
l_recfem = norm(glob_wrap_recfem - glob_ori_recfem) + norm(glob_ins_recfem - glob_wrap_recfem);

% Vastus medialis (femur → patella)
glob_ori_vasmed   = rH + R_bb * vasmed_or;
glob_wrap1_vasmed = rH + R_bb * vasmed_wrap1;
% glob_wrap2_vasmed = rH + R_bb * vasmed_wrap2;
glob_ins_vasmed   = rPat + R_pat * vasmed_in;
l_vasmed = norm(glob_wrap1_vasmed - glob_ori_vasmed) + norm(glob_ins_vasmed - glob_wrap1_vasmed);

% Vastus intermedius
glob_ori_vasint = rH + R_bb * vasint_or;
glob_wrap_vasint = rH + R_bb * vasint_wrap;
glob_ins_vasint = rPat + R_pat * vasint_in;
l_vasint = norm(glob_wrap_vasint - glob_ori_vasint) + norm(glob_ins_vasint - glob_wrap_vasint);

% Vastus lateralis
glob_ori_vaslat   = rH + R_bb * vaslat_or;
glob_wrap1_vaslat = rH + R_bb * vaslat_wrap1;
% glob_wrap2_vaslat = rH + R_bb * vaslat_wrap2;
glob_ins_vaslat   = rPat + R_pat * vaslat_in;
l_vaslat = norm(glob_wrap1_vaslat - glob_ori_vaslat)  + norm(glob_ins_vaslat - glob_wrap1_vaslat);

% Gastrocnemius medialis (femur → enkel)
glob_ori_gasmed = rH + R_bb * gasmed_or;
glob_wrap_gasmed = rH + R_bb * gasmed_wrap;
glob_ins_gasmed = rE; %+ R_ob * gasmed_in; %aangenomen dat het aangrijpt op enkel aangezien alleen effect over knie belangrijk is
l_gasmed = norm(glob_wrap_gasmed - glob_ori_gasmed) + norm(glob_ins_gasmed - glob_wrap_gasmed);

% Gastrocnemius lateralis
glob_ori_gaslat = rH + R_bb * gaslat_or;
glob_wrap_gaslat = rH + R_bb * gaslat_wrap;
glob_ins_gaslat = rE ;%+ R_ob * gaslat_in; %aangenomen dat het aangrijpt op enkel aangezien alleen effect over knie belangrijk is
l_gaslat = norm(glob_wrap_gaslat - glob_ori_gaslat) + norm(glob_ins_gaslat - glob_wrap_gaslat);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contractiesnelheden spieren
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Berekenen van vezellengtes

fl_bicfemsh = l_bicfemsh - pees_sl_bicfemsh;
fl_semimem = l_semimem - pees_sl_semimem;
fl_semiten = l_semiten - pees_sl_semiten;
fl_bicfemlh = l_bicfemlh - pees_sl_bicfemlh;
fl_gmaxsup = l_gmaxsup - pees_sl_gmaxsup;
fl_gmaxmid = l_gmaxmid - pees_sl_gmaxmid;
fl_gmaxinf = l_gmaxinf - pees_sl_gmaxinf;
fl_ilia    = l_ilia - pees_sl_ilia;
fl_psoas   = l_psoas - pees_sl_psoas;
fl_recfem  = l_recfem - pees_sl_recfem;
fl_vasmed  = l_vasmed - pees_sl_vasmed;
fl_vasint  = l_vasint - pees_sl_vasint;
fl_vaslat  = l_vaslat - pees_sl_vaslat;
fl_gasmed  = l_gasmed - pees_sl_gasmed;
fl_gaslat  = l_gaslat - pees_sl_gaslat;

%snelheden van origo en insertie in globaal assenstelsel
v_ori_bicfemsh_glob = vH + ([ey_lok_bb -ex_lok_bb].*fip(4))*bicfemsh_or;
v_ins_bicfemsh_glob = vK + ([ey_lok_ob -ex_lok_ob].*fip(3))*bicfemsh_in;
v_ori_semimem_glob = vPel + ([ey_lok_romp -ex_lok_romp]).*fip(5)*semimem_or;
v_ins_semimem_glob = vK + ([ey_lok_ob -ex_lok_ob]).*fip(3)*semimem_in;
v_ori_semiten_glob = vPel + ([ey_lok_romp -ex_lok_romp]).*fip(5)*semiten_or;
v_ins_semiten_glob = vK + ([ey_lok_ob -ex_lok_ob]).*fip(3)*semiten_in;

% Biceps femoris long head (pelvis → tibia)
v_ori_bicfemlh_glob = vPel + ([ey_lok_romp -ex_lok_romp]).*fip(5)*bicfemlh_or;
v_ins_bicfemlh_glob = vK + ([ey_lok_ob -ex_lok_ob]).*fip(3)*bicfemlh_in;

% Gluteus maximus superior (pelvis → femur)
v_ori_gmaxsup_glob = vPel + ([ey_lok_romp -ex_lok_romp]).*fip(5)*gmaxsup_or;
v_ins_gmaxsup_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*gmaxsup_in;

% Gluteus maximus middle
v_ori_gmaxmid_glob = vPel + ([ey_lok_romp -ex_lok_romp]).*fip(5)*gmaxmid_or;
v_ins_gmaxmid_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*gmaxmid_in;

% Gluteus maximus inferior
v_ori_gmaxinf_glob = vPel + ([ey_lok_romp -ex_lok_romp]).*fip(5)*gmaxinf_or;
v_ins_gmaxinf_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*gmaxinf_in;

% Iliacus (pelvis → femur)
v_ori_ilia_glob = vPel + ([ey_lok_romp -ex_lok_romp]).*fip(5)*ilia_or;
v_ins_ilia_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*ilia_ins;

% Psoas (pelvis → femur)
v_ori_psoas_glob = vPel + ([ey_lok_romp -ex_lok_romp]).*fip(5)*psoas_or;
v_ins_psoas_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*psoas_in;

% Rectus femoris (pelvis → patella via femur wrap)
v_ori_recfem_glob = vPel + ([ey_lok_romp -ex_lok_romp]).*fip(5)*recfem_or;
v_wrap_recfem_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*recfem_wrap;
v_ins_recfem_glob = vPat + ([ey_lok_pat -ex_lok_pat]).*fip(3)*recfem_in;

% Vastus medialis (femur → patella)
v_ori_vasmed_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*vasmed_or;
v_wrap1_vasmed_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*vasmed_wrap1;
% v_wrap2_vasmed_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*vasmed_wrap2;
v_ins_vasmed_glob = vPat + ([ey_lok_pat -ex_lok_pat]).*fip(3)*vasmed_in;

% Vastus intermedius (femur → patella)
v_ori_vasint_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*vasint_or;
v_wrap_vasint_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*vasint_wrap;
v_ins_vasint_glob = vPat + ([ey_lok_pat -ex_lok_pat]).*fip(3)*vasint_in;

% Vastus lateralis (femur → patella)
v_ori_vaslat_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*vaslat_or;
v_wrap1_vaslat_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*vaslat_wrap1;
% v_wrap2_vaslat_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*vaslat_wrap2;
v_ins_vaslat_glob = vPat + ([ey_lok_pat -ex_lok_pat]).*fip(3)*vaslat_in;

% Gastrocnemius medialis (femur → enkel)
v_ori_gasmed_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*gasmed_or;
v_wrap_gasmed_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*gasmed_wrap;
v_ins_gasmed_glob = vE; % enkelpositie

% Gastrocnemius lateralis (femur → enkel)
v_ori_gaslat_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*gaslat_or;
v_wrap_gaslat_glob = vH + ([ey_lok_bb -ex_lok_bb]).*fip(4)*gaslat_wrap;
v_ins_gaslat_glob = vE; % enkelpositie


%eenheidsvector van insertie t.o.v. oorsprong
%en contractiesnelheid
e_inor_bicfemsh = (glob_ins_bicfemsh - glob_ori_bicfemsh)./(norm(glob_ins_bicfemsh - glob_ori_bicfemsh));
v_bicfemsh = dot(e_inor_bicfemsh,v_ori_bicfemsh_glob) - dot(e_inor_bicfemsh,v_ins_bicfemsh_glob);

% --- Geen wrapping ---
[v_bicfemsh, e_inor_bicfemsh] = bereken_spiervector(glob_ori_bicfemsh, glob_ins_bicfemsh, v_ori_bicfemsh_glob, v_ins_bicfemsh_glob);
[v_semimem,  e_inor_semimem]  = bereken_spiervector(glob_ori_semimem,  glob_ins_semimem,  v_ori_semimem_glob,  v_ins_semimem_glob);
[v_semiten,  e_inor_semiten]  = bereken_spiervector(glob_ori_semiten,  glob_ins_semiten,  v_ori_semiten_glob,  v_ins_semiten_glob);
[v_bicfemlh, e_inor_bicfemlh] = bereken_spiervector(glob_ori_bicfemlh, glob_ins_bicfemlh, v_ori_bicfemlh_glob, v_ins_bicfemlh_glob);

[v_gmaxsup,  e_inor_gmaxsup]  = bereken_spiervector(glob_ori_gmaxsup,  glob_ins_gmaxsup,  v_ori_gmaxsup_glob,  v_ins_gmaxsup_glob);
[v_gmaxmid,  e_inor_gmaxmid]  = bereken_spiervector(glob_ori_gmaxmid,  glob_ins_gmaxmid,  v_ori_gmaxmid_glob,  v_ins_gmaxmid_glob);
[v_gmaxinf,  e_inor_gmaxinf]  = bereken_spiervector(glob_ori_gmaxinf,  glob_ins_gmaxinf,  v_ori_gmaxinf_glob,  v_ins_gmaxinf_glob);

[v_ilia,     e_inor_ilia]     = bereken_spiervector(glob_ori_ilia,     glob_ins_ilia,     v_ori_ilia_glob,     v_ins_ilia_glob);
[v_psoas,    e_inor_psoas]    = bereken_spiervector(glob_ori_psoas,    glob_ins_psoas,    v_ori_psoas_glob,    v_ins_psoas_glob);

% --- Met 1 wrap ---
% Rectus femoris
[v1_recfem, e1_recfem] = bereken_spiervector(glob_ori_recfem, glob_wrap_recfem, v_ori_recfem_glob, v_wrap_recfem_glob);
[v2_recfem, e2_recfem] = bereken_spiervector(glob_wrap_recfem, glob_ins_recfem, v_wrap_recfem_glob, v_ins_recfem_glob);
v_recfem = v1_recfem + v2_recfem;

% Vastus intermedius
[v1_vasint, e1_vasint] = bereken_spiervector(glob_ori_vasint, glob_wrap_vasint, v_ori_vasint_glob, v_wrap_vasint_glob);
[v2_vasint, e2_vasint] = bereken_spiervector(glob_wrap_vasint, glob_ins_vasint, v_wrap_vasint_glob, v_ins_vasint_glob);
v_vasint = v1_vasint + v2_vasint;

% Gastrocnemius medialis
[v1_gasmed, e1_gasmed] = bereken_spiervector(glob_ori_gasmed, glob_wrap_gasmed, v_ori_gasmed_glob, v_wrap_gasmed_glob);
[v2_gasmed, e2_gasmed] = bereken_spiervector(glob_wrap_gasmed, glob_ins_gasmed, v_wrap_gasmed_glob, v_ins_gasmed_glob);
v_gasmed = v1_gasmed + v2_gasmed;

% Gastrocnemius lateralis
[v1_gaslat, e1_gaslat] = bereken_spiervector(glob_ori_gaslat, glob_wrap_gaslat, v_ori_gaslat_glob, v_wrap_gaslat_glob);
[v2_gaslat, e2_gaslat] = bereken_spiervector(glob_wrap_gaslat, glob_ins_gaslat, v_wrap_gaslat_glob, v_ins_gaslat_glob);
v_gaslat = v1_gaslat + v2_gaslat;

% --- Met 2 wraps ---
% Vastus medialis
[v1_vasmed, e1_vasmed] = bereken_spiervector(glob_ori_vasmed, glob_wrap1_vasmed, v_ori_vasmed_glob, v_wrap1_vasmed_glob);
% [v2_vasmed, e2_vasmed] = bereken_spiervector(glob_wrap1_vasmed, glob_wrap2_vasmed, v_wrap1_vasmed_glob, v_wrap2_vasmed_glob);
[v3_vasmed, e3_vasmed] = bereken_spiervector(glob_wrap1_vasmed, glob_ins_vasmed, v_wrap1_vasmed_glob, v_ins_vasmed_glob);
v_vasmed = v1_vasmed  + v3_vasmed;

% Vastus lateralis
[v1_vaslat, e1_vaslat] = bereken_spiervector(glob_ori_vaslat, glob_wrap1_vaslat, v_ori_vaslat_glob, v_wrap1_vaslat_glob);
% [v2_vaslat, e2_vaslat] = bereken_spiervector(glob_wrap1_vaslat, glob_wrap2_vaslat, v_wrap1_vaslat_glob, v_wrap2_vaslat_glob);
[v3_vaslat, e3_vaslat] = bereken_spiervector(glob_wrap1_vaslat, glob_ins_vaslat, v_wrap1_vaslat_glob, v_ins_vaslat_glob);
v_vaslat = v1_vaslat  + v3_vaslat;


[adot_bicfemsh] = berekenactivatie(v_bicfemsh,a_bicfemsh);
[adot_semimem]  = berekenactivatie(v_semimem,  a_semimem);
[adot_semiten]  = berekenactivatie(v_semiten,  a_semiten);
[adot_bicfemlh] = berekenactivatie(v_bicfemlh, a_bicfemlh);

[adot_gmaxsup]  = berekenactivatie(v_gmaxsup,  a_gmaxsup);
[adot_gmaxmid]  = berekenactivatie(v_gmaxmid,  a_gmaxmid);
[adot_gmaxinf]  = berekenactivatie(v_gmaxinf,  a_gmaxinf);

[adot_ilia]     = berekenactivatie(v_ilia,     a_ilia);
[adot_psoas]    = berekenactivatie(v_psoas,    a_psoas);

[adot_recfem]   = berekenactivatie(v_recfem,   a_recfem);
[adot_vasmed]   = berekenactivatie(v_vasmed,   a_vasmed);
[adot_vasint]   = berekenactivatie(v_vasint,   a_vasint);
[adot_vaslat]   = berekenactivatie(v_vaslat,   a_vaslat);

[adot_gasmed]   = berekenactivatie(v_gasmed,   a_gasmed);
[adot_gaslat]   = berekenactivatie(v_gaslat,   a_gaslat);

      % equation (1) from De Groote et al 2016, original version as published
      %Muscle excitations were bounded between 0 and 1
      %where u is muscle excitation, a is muscle activation,
% Tact = 0.015 s is the activation time constant,
% Tdeact = 0.060 s is the deactivation time constant, and
% b = 0.1 is a parameter determining transition
% smoothness. 

if v_bicfemsh > 0
u = 1;
else 
    u = 0;
end 

a = a_bicfemsh;
Tact = 0.015;
Tdeact = 0.060;
b = 0.1;
f = 0.5*tanh((u-a)*b);  
% then equation (2) from De Groote et al., 2016
adot = ( 1/Tact/(0.5+1.5*a)*(f+0.5) + (0.5+1.5*a)/Tdeact*(-f+0.5) ) * (u-a);


%Momentsarm en moment berekenen
%vector van insertie ten opzichte van het gewricht
r_ins_bicfemsh_knie = glob_ins_bicfemsh - rK;
r_ins_semimem_knie  = glob_ins_semimem - rK;
r_ori_semimem_heup = glob_ori_semimem - rH;
r_ins_semiten_knie  = glob_ins_semiten - rK;
r_ori_semiten_heup  = glob_ori_semiten - rH;
r_ins_bicfemlh_knie = glob_ins_bicfemlh - rK;
r_ori_bicfemlh_heup = glob_ori_bicfemlh - rH;
r_ins_gmaxsup_heup  = glob_ins_gmaxsup - rH;
r_ins_gmaxmid_heup  = glob_ins_gmaxmid - rH;
r_ins_gmaxinf_heup  = glob_ins_gmaxinf - rH;
r_ins_ilia_heup     = glob_ins_ilia - rH;
r_ins_psoas_knie    = glob_ins_psoas - rH;
r_ins_recfem_knie   = glob_ins_recfem -rK;
r_ori_recfem_heup   = glob_ori_recfem -rH;
r_ins_vasmed_knie   = glob_ins_vasmed - rK;
r_ins_vasint_knie   = glob_ins_vasint - rK;
r_ins_vaslat_knie   = glob_ins_vaslat - rK;
r_ins_gasmed_knie   = glob_ins_gasmed - rK;
r_ins_gaslat_knie   = glob_ins_gaslat - rK;

%Voor biarticulaire spieren ook de eenheidsvector van origo ten opzichte
%van insertie berekenen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Er klopt van de glutei totaal niks
%Van de vastii klopt ook niks volgens mij,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LET OP
%HIER KLOPT NATUURLIJK GEEN FUCK VAN ALS ER WRAPPING POINTS ZIJN
e_orin_semimem = (glob_ori_semimem - glob_ins_semimem)./(norm(glob_ori_semimem - glob_ins_semimem));
e_orin_semiten = (glob_ori_semiten - glob_ins_semiten)./(norm(glob_ori_semiten - glob_ins_semiten));
e_orin_bicfemlh = (glob_ori_bicfemlh - glob_ins_bicfemlh)./(norm(glob_ori_bicfemlh - glob_ins_bicfemlh));


%formules voor kracht komen uit syllabus van mechanische analyse

F0_bicfemsh = Fmax_bicfemsh*exp(-2*(fl_bicfemsh/opt_fl_bicfemsh-1)^2);
 

if v_bicfemsh > 0
    F1_bicfemsh = F0_bicfemsh*(0.3125/(v_bicfemsh/vmax_bicfemsh+0.25)-0.25);

else
    F1_bicfemsh = F0_bicfemsh;
end
F_bicfemsh = F1_bicfemsh*a;
mom_bicfemsh = cross2d(r_ins_bicfemsh_knie,(F_bicfemsh*e_inor_bicfemsh));
d_bicfemsh = cross2d(r_ins_bicfemsh_knie,e_inor_bicfemsh);







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

%INVOEREN VAN BERKENDE MOMENTEN
%HET IS BELANGRIJK DAT ER GOED OPGELET WORDT MET DE CONVENTIE VOOR DE
%MOMENTEN, bijvoorbeeld M of -M afhankelijk van hoe je het moment van de spier
%berekend hebt en of je het moet zien als bijvoorbeeld M12 of -M12

V(18) = V(18) + mom_bicfemsh;



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

[segdynstatedot, Vnew, succes] = segdyn(segdynstate, parms.segparms, K, V, Aconstraints, Bconstraints);
% als er kinematische constraints zijn dan Aconstraints en Bconstraints toevoegen in bovenstaande aanroep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: statedot maken uit segdynstatedot
% in veel gevallen zullen segdynstatedot en statedot identiek zijn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

statedot = segdynstatedot;
statedot = [segdynstatedot; adot_bicfemsh; adot_semimem; adot_semiten; adot_bicfemlh; adot_gmaxsup; adot_gmaxmid; adot_gmaxinf;adot_ilia; adot_psoas; adot_recfem; adot_vasmed; adot_vasint; adot_vaslat; adot_gasmed;  adot_gaslat];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: uitgangen y berekenen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y is een KOLOMvector van willekeurige lengte die alle interessant geachte afhankelijke variabelen
% bevat (bijv reactiekrachten, versnellingen, energietermen, zwaartepuntpositie etc etc)

if parms.calculate_outputs
    y = [Vnew, rPat(1),rPat(2), l_vaslat];%glob_ori_vaslat(1),glob_ori_vaslat(2),glob_wrap1_vaslat(1),glob_wrap1_vaslat(2),glob_ins_vaslat(1),glob_ins_vaslat(2)];%rH(1),rH(2),rPel(1),rPel(2)];
end


% adot_bicfemsh, adot_semimem,adot_semiten, adot_bicfemlh,adot_gmaxsup, adot_gmaxmid, adot_gmaxinf,adot_ilia, adot_psoas, adot_recfem,adot_vasmed,adot_vasint, adot_vaslat,adot_gasmed, adot_gaslat];