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


%momentsarmen uit soest en casius 2000
%momentsarmen knie: 
d_Glut_max_heup = 0.062;
d_Vasti_knie = 0.042;
d_Recfem_knie = 0.042;
d_Recfem_heup = 0.035;
d_Gas_knie = 0.014;
d_Biceps_fem_knie = 0.026;
d_Biceps_fem_heup = 0.077;
d_ilia_heup = 0.050;
% 
% d_Glut_max_heup = 0.062;
% d_Vasti_knie = -0.042;
% d_Recfem_knie = -0.042;
% d_Recfem_heup = -0.035;
% d_Gas_knie = 0.014;
% d_Biceps_fem_knie = 0.026;
% d_Biceps_fem_heup = 0.077;
% d_ilia_heup = -0.050;

delta_phi_heup =  -(fi(5)-fi(4)) ;
delta_phi_knie =  (fi(4)-fi(3)) ;


%Lengtes bij rechtop staan
l0_bicfemsh = 0.293788308140402;
l0_semimem = 0.457841193667097;
l0_semiten = 0.422912915015888;
l0_bicfemlh = 0.480723839398249;
l0_gmaxsup = 0.167855517960221;
l0_gmaxmid = 0.147397625995096;
l0_gmaxinf = 0.161049300957849;
l0_ilia = 0.163408728960102;
l0_psoas = 0.211700922552417;
l0_recfem = 0.464703337159607;
l0_vasmed = 0.149299759544040;
l0_vasint = 0.236747224758782;
l0_vaslat = 0.237233408854345;
l0_gasmed = 0.494952705418321;
l0_gaslat = 0.493112616009206;

l_bicfemsh = l0_bicfemsh - d_Biceps_fem_knie*delta_phi_knie;
l_semimem = l0_semimem - d_Biceps_fem_knie*delta_phi_knie + d_Biceps_fem_heup*delta_phi_heup;
l_semiten = l0_semiten - d_Biceps_fem_knie*delta_phi_knie + d_Biceps_fem_heup*delta_phi_heup ;
l_bicfemlh = l0_bicfemlh - d_Biceps_fem_knie*delta_phi_knie + d_Biceps_fem_heup*delta_phi_heup;
l_gmaxsup = l0_gmaxsup + d_Glut_max_heup*delta_phi_heup;
l_gmaxmid = l0_gmaxmid + d_Glut_max_heup*delta_phi_heup;
l_gmaxinf = l0_gmaxinf + d_Glut_max_heup*delta_phi_heup;
l_ilia = l0_ilia - d_ilia_heup*delta_phi_heup;
l_psoas = l0_psoas - d_ilia_heup*delta_phi_heup;
l_recfem = l0_recfem - d_Recfem_heup*delta_phi_heup + d_Recfem_knie*delta_phi_knie;
l_vasmed = l0_vasmed + d_Vasti_knie*delta_phi_knie;
l_vasint = l0_vasint + d_Vasti_knie*delta_phi_knie;
l_vaslat = l0_vaslat + d_Vasti_knie*delta_phi_knie;
l_gasmed = l0_gasmed - d_Gas_knie*delta_phi_knie;
l_gaslat = l0_gaslat - d_Gas_knie*delta_phi_knie;

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contractiesnelheden spieren
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Berekenen van vezellengtes

fip_heup = -(fip(5) - fip(4));

fip_knie = (fip(4) - fip(3));

v_bicfemsh = d_Biceps_fem_knie*fip_knie;
v_semimem = d_Biceps_fem_knie*fip_knie - d_Biceps_fem_heup*fip_heup;
v_semiten = d_Biceps_fem_knie*fip_knie - d_Biceps_fem_heup*fip_heup;
v_bicfemlh = d_Biceps_fem_knie*fip_knie - d_Biceps_fem_heup*fip_heup;
v_gmaxsup = -d_Glut_max_heup*fip_heup;
v_gmaxmid = -d_Glut_max_heup*fip_heup;
v_gmaxinf = -d_Glut_max_heup*fip_heup;
v_ilia = d_ilia_heup*fip_heup;
v_psoas = d_ilia_heup*fip_heup;
v_recfem = d_Recfem_heup*fip_heup - d_Recfem_knie*fip_knie;
v_vasmed = -d_Vasti_knie*fip_knie;
v_vasint = -d_Vasti_knie*fip_knie;
v_vaslat = -d_Vasti_knie*fip_knie;
v_gasmed = d_Gas_knie*fip_knie;
v_gaslat = d_Gas_knie*fip_knie;

vspier = [v_bicfemsh, v_semimem, v_semiten, v_bicfemlh, ...
          v_gmaxsup, v_gmaxmid, v_gmaxinf, v_ilia, ...
          v_psoas, v_recfem, v_vasmed, v_vasint, ...
          v_vaslat, v_gasmed, v_gaslat];



%ZONDER ACTIVATIEDYNAMICA
[a_bicfemsh] = berekenactivatie2(v_bicfemsh);
%[a_semimem]  = berekenactivatie2(v_semimem);
%[a_semiten]  = berekenactivatie2(v_semiten);
%[a_bicfemlh] = berekenactivatie2(v_bicfemlh);

%[a_gmaxsup]  = berekenactivatie2(v_gmaxsup);%
%[a_gmaxmid]  = berekenactivatie2(v_gmaxmid);
%[a_gmaxinf]  = berekenactivatie2(v_gmaxinf);

%[a_ilia]     = berekenactivatie2(v_ilia);
%[a_psoas]    = berekenactivatie2(v_psoas);

%[a_recfem]   = berekenactivatie2(v_recfem);
[a_vasmed]   = berekenactivatie2(v_vasmed);
[a_vasint]   = berekenactivatie2(v_vasint);
[a_vaslat]   = berekenactivatie2(v_vaslat);

[a_gasmed]   = berekenactivatie2(v_gasmed);
[a_gaslat]   = berekenactivatie2(v_gaslat);
      % equation (1) from De Groote et al 2016, original version as published
      %Muscle excitations were bounded between 0 and 1
      %where u is muscle excitation, a is muscle activation,
% Tact = 0.015 s is the activation time constant,
% Tdeact = 0.060 s is the deactivation time constant, and
% b = 0.1 is a parameter determining transition
% smoothness. 
a_tot_knie = [a_vasmed, a_vasint, a_vaslat, a_gasmed, a_gaslat];

% if v_bicfemsh > 0
% u = 1;
% else 
%     u = 0;
% end 
% 
% a = a_bicfemsh;
% Tact = 0.015;
% Tdeact = 0.060;
% b = 0.1;
% f = 0.5*tanh((u-a)*b);  
% % then equation (2) from De Groote et al., 2016
% adot = ( 1/Tact/(0.5+1.5*a)*(f+0.5) + (0.5+1.5*a)/Tdeact*(-f+0.5) ) * (u-a);



%Voor biarticulaire spieren ook de eenheidsvector van origo ten opzichte
%van insertie berekenen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Er klopt van de glutei totaal niks
%Van de vastii klopt ook niks volgens mij,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LET OP
%HIER KLOPT NATUURLIJK GEEN FUCK VAN ALS ER WRAPPING POINTS ZIJN
% e_orin_bicfemsh = (glob_ori_bicfemsh - glob_ins_bicfemsh)./(norm(glob_ori_bicfemsh - glob_ins_bicfemsh));
% e_orin_semimem = (glob_ori_semimem - glob_ins_semimem)./(norm(glob_ori_semimem - glob_ins_semimem));
% e_orin_semiten = (glob_ori_semiten - glob_ins_semiten)./(norm(glob_ori_semiten - glob_ins_semiten));
% e_orin_bicfemlh = (glob_ori_bicfemlh - glob_ins_bicfemlh)./(norm(glob_ori_bicfemlh - glob_ins_bicfemlh));


%formules voor kracht komen uit syllabus van mechanische analyse









Moment = cos(fi(5))*0.54*4*9.81 + cos(fi(5))*(0.2214)*parms.segparms.m(5)*9.81  ; %-6.547946504007387; %-12.454934381976590; %71.660921167714820;%-6.547946504007387;


    a_bicfemlh = 0;
    a_semimem  = 0;
    a_semiten  = 0;
    a_gmaxsup  = 0;
    a_gmaxmid  = 0;
    a_gmaxinf  = 0;
    a_ilia     = 0;
    a_psoas    = 0;
    a_recfem   = 0;

Fmaxfv_psoas = berekenmaxspierkracht(Fmax_psoas, v_psoas, vmax_psoas, opt_fl_psoas, fl_psoas);
Fmaxfv_recfem = berekenmaxspierkracht(Fmax_recfem, v_recfem, vmax_recfem, opt_fl_recfem, fl_recfem);
Fmaxfv_ilia = berekenmaxspierkracht(Fmax_ilia, v_ilia, vmax_ilia, opt_fl_ilia, fl_ilia);
Fmaxfv_gmaxinf = berekenmaxspierkracht(Fmax_gmaxinf,v_gmaxinf,vmax_gmaxinf,opt_fl_gmaxinf,fl_gmaxinf);
Fmaxfv_gmaxmid = berekenmaxspierkracht(Fmax_gmaxmid,v_gmaxmid,vmax_gmaxmid,opt_fl_gmaxmid,fl_gmaxmid);
Fmaxfv_gmaxsup = berekenmaxspierkracht(Fmax_gmaxsup,v_gmaxsup,vmax_gmaxsup,opt_fl_gmaxsup,fl_gmaxsup);
Fmaxfv_semiten = berekenmaxspierkracht(Fmax_semiten,v_semiten,vmax_semiten,opt_fl_semiten,fl_semiten);
Fmaxfv_semimem = berekenmaxspierkracht(Fmax_semimem,v_semimem,vmax_semimem,opt_fl_semimem,fl_semimem);
Fmaxfv_bicfemlh = berekenmaxspierkracht(Fmax_bicfemlh,v_bicfemlh,vmax_bicfemlh,opt_fl_bicfemlh,fl_bicfemlh);

%kleine waardes toegevoegd zodat er een standaard activatievolgorde is
%vanwege het gelijk zijn van de potentiële vermogens
v_gmax = v_gmaxsup;
P_bicfemlh = 0.000000001+(Moment / d_Biceps_fem_heup) * v_bicfemlh;
P_semimem = 0.00000001+(Moment / d_Biceps_fem_heup) * v_semimem;
P_semiten = 0.000001+(Moment / d_Biceps_fem_heup) * v_semiten;
P_gmaxsup  = 0.00000000001+(Moment / d_Glut_max_heup) * v_gmaxsup;
P_gmaxmid  = 0.00000001+(Moment / d_Glut_max_heup) * v_gmaxmid;
P_gmaxinf  = 0.000001+(Moment / d_Glut_max_heup) * v_gmaxinf;
P_spier = [P_bicfemlh, P_semimem, P_semiten, P_gmaxsup, P_gmaxmid,P_gmaxinf];

fi_doel = 0.4*pi + 0.07*sin(4*pi*t+0.6*pi);
kp_heup = 500;
kd_heup = 50;
M_heup_control = -kp_heup*(fi(5)-fi_doel) -kd_heup*fip(5);
% Moment = M_heup_control + Moment;

%% === Input: moment dat geleverd moet worden ===
Moment_nodig = abs(Moment); % extensiemoment positief genomen

%% === Vermogensschatting per spier ===
P_bicfemlh = 1e-9  + (Moment / d_Biceps_fem_heup) * v_bicfemlh;
P_semimem  = 1e-8  + (Moment / d_Biceps_fem_heup) * v_semimem;
P_semiten  = 1e-6  + (Moment / d_Biceps_fem_heup) * v_semiten;
P_gmaxsup  = 1e-11 + (Moment / d_Glut_max_heup)   * v_gmaxsup;
P_gmaxmid  = 1e-8  + (Moment / d_Glut_max_heup)   * v_gmaxmid;
P_gmaxinf  = 1e-6  + (Moment / d_Glut_max_heup)   * v_gmaxinf;

P_spier = [P_bicfemlh, P_semimem, P_semiten, P_gmaxsup, P_gmaxmid, P_gmaxinf];

%% === Sorteren op vermogen (hoogste eerst) ===
[Ps_sorted, idx_sorted] = sort(P_spier, 'descend');

%% === Bereken maximale kracht per spier ===
Fmaxfv_bicfemlh = berekenmaxspierkracht(Fmax_bicfemlh, v_bicfemlh, vmax_bicfemlh, opt_fl_bicfemlh, fl_bicfemlh);
Fmaxfv_semimem  = berekenmaxspierkracht(Fmax_semimem,  v_semimem,  vmax_semimem,  opt_fl_semimem,  fl_semimem);
Fmaxfv_semiten  = berekenmaxspierkracht(Fmax_semiten,  v_semiten,  vmax_semiten,  opt_fl_semiten,  fl_semiten);
Fmaxfv_gmaxsup  = berekenmaxspierkracht(Fmax_gmaxsup,  v_gmaxsup,  vmax_gmaxsup,  opt_fl_gmaxsup,  fl_gmaxsup);
Fmaxfv_gmaxmid  = berekenmaxspierkracht(Fmax_gmaxmid,  v_gmaxmid,  vmax_gmaxmid,  opt_fl_gmaxmid,  fl_gmaxmid);
Fmaxfv_gmaxinf  = berekenmaxspierkracht(Fmax_gmaxinf,  v_gmaxinf,  vmax_gmaxinf,  opt_fl_gmaxinf,  fl_gmaxinf);

Fmaxfv = [Fmaxfv_bicfemlh, Fmaxfv_semimem, Fmaxfv_semiten, ...
          Fmaxfv_gmaxsup, Fmaxfv_gmaxmid, Fmaxfv_gmaxinf];

d_spier = [d_Biceps_fem_heup, d_Biceps_fem_heup, d_Biceps_fem_heup, ...
            d_Glut_max_heup,  d_Glut_max_heup,  d_Glut_max_heup];

%% === Input: moment dat geleverd moet worden ===
Moment_nodig = abs(Moment); % positief extensiemoment

%% === Vermogensschatting per spier (absolute waarden) ===
P_bicfemlh = 1e-9  + ((Moment / d_Biceps_fem_heup) * v_bicfemlh);
P_semimem  = 1e-8  + ((Moment / d_Biceps_fem_heup) * v_semimem);
P_semiten  = 1e-6  + ((Moment / d_Biceps_fem_heup) * v_semiten);
P_gmaxsup  = 1e-11 + ((Moment / d_Glut_max_heup)   * v_gmaxsup);
P_gmaxmid  = 1e-8  + ((Moment / d_Glut_max_heup)   * v_gmaxmid);
P_gmaxinf  = 1e-6  + ((Moment / d_Glut_max_heup)   * v_gmaxinf);

P_spier = [P_bicfemlh, P_semimem, P_semiten, P_gmaxsup, P_gmaxmid, P_gmaxinf];

%% === Sorteren op vermogen (hoogste eerst) ===
[~, idx_sorted] = sort(P_spier, 'descend');

%% === Bereken maximale kracht per spier ===
Fmaxfv_bicfemlh = berekenmaxspierkracht(Fmax_bicfemlh, v_bicfemlh, vmax_bicfemlh, opt_fl_bicfemlh, fl_bicfemlh);
Fmaxfv_semimem  = berekenmaxspierkracht(Fmax_semimem,  v_semimem,  vmax_semimem,  opt_fl_semimem,  fl_semimem);
Fmaxfv_semiten  = berekenmaxspierkracht(Fmax_semiten,  v_semiten,  vmax_semiten,  opt_fl_semiten,  fl_semiten);
Fmaxfv_gmaxsup  = berekenmaxspierkracht(Fmax_gmaxsup,  v_gmaxsup,  vmax_gmaxsup,  opt_fl_gmaxsup,  fl_gmaxsup);
Fmaxfv_gmaxmid  = berekenmaxspierkracht(Fmax_gmaxmid,  v_gmaxmid,  vmax_gmaxmid,  opt_fl_gmaxmid,  fl_gmaxmid);
Fmaxfv_gmaxinf  = berekenmaxspierkracht(Fmax_gmaxinf,  v_gmaxinf,  vmax_gmaxinf,  opt_fl_gmaxinf,  fl_gmaxinf);

Fmaxfv = [Fmaxfv_bicfemlh, Fmaxfv_semimem, Fmaxfv_semiten, ...
          Fmaxfv_gmaxsup, Fmaxfv_gmaxmid, Fmaxfv_gmaxinf];

d_spier = abs([d_Biceps_fem_heup, d_Biceps_fem_heup, d_Biceps_fem_heup, ...
               d_Glut_max_heup,  d_Glut_max_heup,  d_Glut_max_heup]);

%% === Initialisatie ===
a = zeros(1,6);
RemainingMoment = Moment_nodig;

%% === Cumulatieve rekrutering ===
for i = 1:length(idx_sorted)
    spier_idx = idx_sorted(i);
    Fmax_i = Fmaxfv(spier_idx);
    d_i = d_spier(spier_idx);

    % Maximale momentcapaciteit van deze spier
    Mmax_i = abs(Fmax_i * d_i);

    % Benodigde activatie
    a_needed = RemainingMoment / Mmax_i;

    if a_needed <= 1
        a(spier_idx) = a_needed;
        RemainingMoment = 0;
        break;
    else
        a(spier_idx) = 1;
        RemainingMoment = RemainingMoment - Mmax_i;
    end

    if RemainingMoment <= 1e-6
        break;
    end
end

%% === Controle ===
a = max(a, 0); % geen negatieve activaties
a(a > 1) = 1;

%% === Toewijzing ===
a_bicfemlh = a(1);
a_semimem  = a(2);
a_semiten  = a(3);
a_gmaxsup  = a(4);
a_gmaxmid  = a(5);
a_gmaxinf  = a(6);
a_ilia   = 0;
a_psoas  = 0;
a_recfem = 0;

%% === Controleoutput ===
Moment_gelev = sum(a .* Fmaxfv .* d_spier);
fprintf('Moment gevraagd: %.2f Nm, geleverd: %.2f Nm\n', Moment_nodig, Moment_gelev);
fprintf('Activaties (positief, 0–1):\n');
disp(a);





    a_all = [a_bicfemlh, a_semimem, a_semiten, a_gmaxsup, ...
         a_gmaxmid, a_gmaxinf, a_ilia, a_psoas, a_recfem];

    a_heup = [a_semimem, a_semiten, a_bicfemlh, a_gmaxsup, a_gmaxmid, a_gmaxinf];
a_knie = [a_bicfemsh, a_semimem, a_semiten, a_bicfemlh, a_recfem,a_vasmed,a_vasint,a_vaslat,a_gasmed,a_gaslat];


[F_bicfemsh] = berekenspierkracht(Fmax_bicfemsh, a_bicfemsh, v_bicfemsh, vmax_bicfemsh, opt_fl_bicfemsh, fl_bicfemsh);
[F_semimem] = berekenspierkracht(Fmax_semimem, a_semimem, v_semimem, vmax_semimem, opt_fl_semimem, fl_semimem);
[F_semiten] = berekenspierkracht(Fmax_semiten, a_semiten, v_semiten, vmax_semiten, opt_fl_semiten, fl_semiten);
[F_bicfemlh] = berekenspierkracht(Fmax_bicfemlh, a_bicfemlh, v_bicfemlh, vmax_bicfemlh, opt_fl_bicfemlh, fl_bicfemlh);
[F_gmaxsup] = berekenspierkracht(Fmax_gmaxsup, a_gmaxsup, v_gmaxsup, vmax_gmaxsup, opt_fl_gmaxsup, fl_gmaxsup);
[F_gmaxmid] = berekenspierkracht(Fmax_gmaxmid, a_gmaxmid, v_gmaxmid, vmax_gmaxmid, opt_fl_gmaxmid, fl_gmaxmid);
[F_gmaxinf] = berekenspierkracht(Fmax_gmaxinf, a_gmaxinf, v_gmaxinf, vmax_gmaxinf, opt_fl_gmaxinf, fl_gmaxinf);
[F_ilia] = berekenspierkracht(Fmax_ilia, a_ilia, v_ilia, vmax_ilia, opt_fl_ilia, fl_ilia);
[F_psoas] = berekenspierkracht(Fmax_psoas, a_psoas, v_psoas, vmax_psoas, opt_fl_psoas, fl_psoas);
[F_recfem] = berekenspierkracht(Fmax_recfem, a_recfem, v_recfem, vmax_recfem, opt_fl_recfem, fl_recfem);
[F_vasmed] = berekenspierkracht(Fmax_vasmed, a_vasmed, v_vasmed, vmax_vasmed, opt_fl_vasmed, fl_vasmed);
[F_vasint] = berekenspierkracht(Fmax_vasint, a_vasint, v_vasint, vmax_vasint, opt_fl_vasint, fl_vasint);
[F_vaslat] = berekenspierkracht(Fmax_vaslat, a_vaslat, v_vaslat, vmax_vaslat, opt_fl_vaslat, fl_vaslat);
[F_gasmed] = berekenspierkracht(Fmax_gasmed, a_gasmed, v_gasmed, vmax_gasmed, opt_fl_gasmed, fl_gasmed);
[F_gaslat] = berekenspierkracht(Fmax_gaslat, a_gaslat, v_gaslat, vmax_gaslat, opt_fl_gaslat, fl_gaslat);




%Momenten knie berekenen
M_bicfemsh_knie = F_bicfemsh*d_Biceps_fem_knie;
M_semimem_knie = F_semimem*d_Biceps_fem_knie;
M_semiten_knie = F_semiten*d_Biceps_fem_knie;
M_bicfemlh_knie = F_bicfemlh*d_Biceps_fem_knie;
M_recfem_knie = -F_recfem*d_Recfem_knie;
M_vasmed_knie = -F_vasmed*d_Vasti_knie;
M_vasint_knie = -F_vasint*d_Vasti_knie;
M_vaslat_knie = -F_vaslat*d_Vasti_knie;
M_gasmed_knie = F_gasmed*d_Gas_knie;
M_gaslat_knie = F_gaslat*d_Gas_knie;

M_knie_tot = M_bicfemsh_knie + M_semimem_knie + M_semiten_knie + ...
             M_bicfemlh_knie + M_recfem_knie + M_vasmed_knie + ...
             M_vasint_knie + M_vaslat_knie + M_gasmed_knie + ...
             M_gaslat_knie;


M_knie = [ M_bicfemsh_knie, M_semimem_knie, M_semiten_knie, ...
           M_bicfemlh_knie, M_recfem_knie, M_vasmed_knie, ...
           M_vasint_knie, M_vaslat_knie, M_gasmed_knie, ...
           M_gaslat_knie ];

%Momenten heup berekeken
M_bicfemlh_heup = F_bicfemlh*d_Biceps_fem_heup;
M_semimem_heup = F_semimem*d_Biceps_fem_heup;
M_semiten_heup = F_semiten*d_Biceps_fem_heup;
M_gmaxsup_heup = F_gmaxsup*d_Glut_max_heup;
M_gmaxmid_heup = F_gmaxmid*d_Glut_max_heup;
M_gmaxinf_heup = F_gmaxinf*d_Glut_max_heup;
M_ilia_heup = -F_ilia*d_ilia_heup;
M_psoas_heup = -F_psoas*d_ilia_heup;
M_recfem_heup = -F_recfem*d_Recfem_heup;

M_heup_tot = (M_bicfemlh_heup + M_semimem_heup + M_semiten_heup + ...
            M_gmaxsup_heup + M_gmaxmid_heup + M_gmaxinf_heup + ...
            M_ilia_heup + M_psoas_heup + M_recfem_heup);

M_spier_heup = [ ...
    M_bicfemlh_heup, ...
    M_semimem_heup, ...
    M_semiten_heup, ...
    M_gmaxsup_heup, ...
    M_gmaxmid_heup, ...
    M_gmaxinf_heup, ...
    M_ilia_heup, ...
    M_psoas_heup, ...
    M_recfem_heup ...
];




Fmaxfv_psoas = berekenmaxspierkracht(Fmax_psoas, v_psoas, vmax_psoas, opt_fl_psoas, fl_psoas);
Fmaxfv_recfem = berekenmaxspierkracht(Fmax_recfem, v_recfem, vmax_recfem, opt_fl_recfem, fl_recfem);
Fmaxfv_ilia = berekenmaxspierkracht(Fmax_ilia, v_ilia, vmax_ilia, opt_fl_ilia, fl_ilia);
Fmaxfv_gmaxinf = berekenmaxspierkracht(Fmax_gmaxinf,v_gmaxinf,vmax_gmaxinf,opt_fl_gmaxinf,fl_gmaxinf);
Fmaxfv_gmaxmid = berekenmaxspierkracht(Fmax_gmaxmid,v_gmaxmid,vmax_gmaxmid,opt_fl_gmaxmid,fl_gmaxmid);
Fmaxfv_gmaxsup = berekenmaxspierkracht(Fmax_gmaxsup,v_gmaxsup,vmax_gmaxsup,opt_fl_gmaxsup,fl_gmaxsup);
Fmaxfv_semiten = berekenmaxspierkracht(Fmax_semiten,v_semiten,vmax_semiten,opt_fl_semiten,fl_semiten);
Fmaxfv_semimem = berekenmaxspierkracht(Fmax_semimem,v_semimem,vmax_semimem,opt_fl_semimem,fl_semimem);
Fmaxfv_bicfemlh = berekenmaxspierkracht(Fmax_bicfemlh,v_bicfemlh,vmax_bicfemlh,opt_fl_bicfemlh,fl_bicfemlh);


potmomflex = -(Fmaxfv_psoas*d_ilia_heup + Fmaxfv_ilia*d_ilia_heup + Fmaxfv_recfem*d_Recfem_heup);
potmomext = -(Fmaxfv_gmaxinf*d_Glut_max_heup + Fmaxfv_gmaxmid*d_Glut_max_heup + Fmaxfv_gmaxsup*d_Glut_max_heup + Fmaxfv_semiten*d_Biceps_fem_heup + Fmaxfv_semimem*d_Biceps_fem_heup + Fmaxfv_bicfemlh*d_Biceps_fem_heup);


Mom_x_V = [
    v_bicfemsh  * F_bicfemsh;
    v_semimem   * F_semimem;
    v_semiten   * F_semiten;
    v_bicfemlh  * F_bicfemlh;
    v_gmaxsup   * F_gmaxsup;
    v_gmaxmid   * F_gmaxmid;
    v_gmaxinf   * F_gmaxinf;
    v_ilia      * F_ilia;
    v_psoas     * F_psoas;
    v_recfem    * F_recfem;
    v_vasmed    * F_vasmed;
    v_vasint    * F_vasint;
    v_vaslat    * F_vaslat;
    v_gasmed    * F_gasmed;
    v_gaslat    * F_gaslat
]';

Ptot = sum(Mom_x_V);


Precfem =  v_recfem * F_recfem;
Psemimem =    v_semimem   * F_semimem;
Psemiten =    v_semiten   * F_semiten;
Pbicfemlh =     v_bicfemlh  * F_bicfemlh;

Pbiarti = [Precfem, Psemimem, Psemiten, Pbicfemlh];

Pgmaxsup =     v_gmaxsup   * F_gmaxsup;
Pgmaxmid =   v_gmaxmid   * F_gmaxmid;
Pgmaxinf =     v_gmaxinf   * F_gmaxinf;


Pheuprest = [Pgmaxsup,Pgmaxmid,Pgmaxinf];









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
V(18) = M_knie_tot;
V(19) = M_heup_tot;

Pknie = fip_knie*V(18);
Pheup = -fip_heup*V(19);

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





%vastzetten positie heup
Aconstraints(5,:) = zeros(1,47);
Aconstraints(6,:) = zeros(1,47);

Aconstraints(5,40:43) = [-parms.segparms.L(1)*sin(fi(1)) -parms.segparms.L(2)*sin(fi(2)) -parms.segparms.L(3)*sin(fi(3)) -parms.segparms.L(4)*sin(fi(4))]; 
Aconstraints(6,40:43) = [parms.segparms.L(1)*cos(fi(1)) parms.segparms.L(2)*cos(fi(2)) parms.segparms.L(3)*cos(fi(3)) parms.segparms.L(4)*cos(fi(4))]; 

Aconstraints(5,46) = 1;
Aconstraints(6,47) = 1;

Bconstraints(5,1) = parms.segparms.L(1)*cos(fi(1))*(fip(1)^2) + parms.segparms.L(2)*cos(fi(2))*(fip(2)^2) + parms.segparms.L(3)*cos(fi(3))*(fip(3)^2) + parms.segparms.L(4)*cos(fi(4))*(fip(4)^2);%-parms.segparms.L(5)*cos(fi(5))*(fip(5)^2) -parms.segparms.L(6)*cos(fi(6))*(fip(6)^2);
Bconstraints(6,1) = parms.segparms.L(1)*sin(fi(1))*(fip(1)^2) + parms.segparms.L(2)*sin(fi(2))*(fip(2)^2) + parms.segparms.L(3)*sin(fi(3))*(fip(3)^2) + parms.segparms.L(4)*sin(fi(4))*(fip(4)^2);%-parms.segparms.L(5)*sin(fi(5))*(fip(5)^2) -parms.segparms.L(6)*sin(fi(6))*(fip(6)^2);

K(26) = 0;
V(26) = NaN;
K(32) = 0;
V(32) = NaN;



%Moment van externe krachten op romp toevoegen vanuit aangrijpingspunt in
%heup
Aconstraints(7,:) = zeros(1,47);
Aconstraints(7,26) = sin(fi(5))*0.2214;
Aconstraints(7,32) = -cos(fi(5))*0.2214;
Aconstraints(7,38) = -1;
Bconstraints(7,1) = 0;
K(38) = 0;
V(38) = NaN;



%vastzetten positie hand
% Aconstraints(7,:) = zeros(1,47);
% Aconstraints(8,:) = zeros(1,47);
% 
% Aconstraints(7,40:45) = [-parms.segparms.L(1)*sin(fi(1)) -parms.segparms.L(2)*sin(fi(2)) -parms.segparms.L(3)*sin(fi(3)) -parms.segparms.L(4)*sin(fi(4)) -parms.segparms.L(6)*sin(fi(6)) -parms.segparms.L(6)*sin(fi(6))]; 
% Aconstraints(8,40:45) = [parms.segparms.L(1)*cos(fi(1)) parms.segparms.L(2)*cos(fi(2)) parms.segparms.L(3)*cos(fi(3)) parms.segparms.L(4)*cos(fi(4)) parms.segparms.L(5)*cos(fi(5)) parms.segparms.L(6)*cos(fi(6))]; 
% 
% Aconstraints(7,46) = 1;
% Aconstraints(8,47) = 1;
% 
% Bconstraints(7,1) = parms.segparms.L(1)*cos(fi(1))*(fip(1)^2) + parms.segparms.L(2)*cos(fi(2))*(fip(2)^2) + parms.segparms.L(3)*cos(fi(3))*(fip(3)^2) + parms.segparms.L(4)*cos(fi(4))*(fip(4)^2) + parms.segparms.L(5)*cos(fi(5))*(fip(5)^2) + parms.segparms.L(6)*cos(fi(6))*(fip(6)^2);
% Bconstraints(8,1) = parms.segparms.L(1)*sin(fi(1))*(fip(1)^2) + parms.segparms.L(2)*sin(fi(2))*(fip(2)^2) + parms.segparms.L(3)*sin(fi(3))*(fip(3)^2) + parms.segparms.L(4)*sin(fi(4))*(fip(4)^2) + parms.segparms.L(5)*sin(fi(5))*(fip(5)^2) + parms.segparms.L(6)*sin(fi(6))*(fip(6)^2);
% 
% K(7) = 0;
% V(7) = NaN;
% K(14) = 0;
% V(14) = NaN;

%vastzetten oriëntatie romp
% Aconstraints(7,:) = [zeros(1,43) 1 zeros(1,3)];
% Bconstraints(7,1) = 0;
% 
% K(38) = 0;
% V(38) = NaN;


%INVOEREN VAN BERKENDE MOMENTEN
%HET IS BELANGRIJK DAT ER GOED OPGELET WORDT MET DE CONVENTIE VOOR DE
%MOMENTEN, bijvoorbeeld M of -M afhankelijk van hoe je het moment van de spier
%berekend hebt en of je het moet zien als bijvoorbeeld M12 of -M12





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
%statedot = [segdynstatedot; adot_bicfemsh; adot_semimem; adot_semiten; adot_bicfemlh; adot_gmaxsup; adot_gmaxmid; adot_gmaxinf;adot_ilia; adot_psoas; adot_recfem; adot_vasmed; adot_vasint; adot_vaslat; adot_gasmed;  adot_gaslat];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% applicatiespecifiek: uitgangen y berekenen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y is een KOLOMvector van willekeurige lengte die alle interessant geachte afhankelijke variabelen
% bevat (bijv reactiekrachten, versnellingen, energietermen, zwaartepuntpositie etc etc)

if parms.calculate_outputs
    y = [Vnew,vspier,a_all,M_knie,M_spier_heup,a_heup,P_spier,Ptot,Pknie,Pheup,Pbiarti,Pheuprest,fip_heup, fip_knie,v_semimem,v_gmax,Ptot];%,v_gmaxsup,a_all];
end


% adot_bicfemsh, adot_semimem,adot_semiten, adot_bicfemlh,adot_gmaxsup, adot_gmaxmid, adot_gmaxinf,adot_ilia, adot_psoas, adot_recfem,adot_vasmed,adot_vasint, adot_vaslat,adot_gasmed, adot_gaslat];