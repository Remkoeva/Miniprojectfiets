function [F1] = berekenmaxspierkracht(Fmax,v_spier,vmax_spier,opt_fl,fl)

%formule om de huidige kracht in de spier te berekenen op basis van de
%maximale kracht, verkortingssnelheid, maximale
%verkortingssnelheid, optimale vezellengte en huidige vezellengte

%formules voor kracht komen uit syllabus van mechanische analyse (2024/2025)
F0 = Fmax*exp(-2*(fl/opt_fl-1)^2);
if v_spier > 0
    F1 = F0*(0.3125/(v_spier/vmax_spier+0.25)-0.25);

else
    F1 = F0;
end




