function [adot_spier] = berekenactivatie(v_spier,a_spier)

if v_spier > 0
u = 1;
else 
    u = 0;
end 

a = a_spier;
Tact = 0.015;
Tdeact = 0.060;
b = 0.1;
f = 0.5*tanh((u-a)*b);  
% then equation (2) from De Groote et al., 2016
adot_spier = ( 1/Tact/(0.5+1.5*a)*(f+0.5) + (0.5+1.5*a)/Tdeact*(-f+0.5) ) * (u-a);