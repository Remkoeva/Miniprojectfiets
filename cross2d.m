function [M] = cross2d(r,F)
%cross2d Berekend het crossproduct van de 2D vectoren r en F (met een x en
%y componenten).
%   input:
%       - (1) r: 2D vector met typisch de plaats van het aangrijpingspunte
%       ten opzichte van het draaipunt
%       - (2) F: 2D vector met typisch de krachvector

% controleer of r en F 1x2 vectoren zijn. Indien niet, voeg bewerking uit
% op alle rijen/kolommen.
[nr,nc] = size(r);
if nr == 2
    % berekening
    M = r(1,:) .* F(2,:) - r(2,:).*F(1,:);
elseif nc == 2
    % berekening
    M = r(:,1) .* F(:,2) - r(:,2).*F(:,1);
else
    % berekening
    M = r(1) * F(2) - r(2)*F(1);
end


end