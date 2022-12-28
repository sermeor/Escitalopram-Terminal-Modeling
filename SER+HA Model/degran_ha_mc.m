%% Function of mast cell degranulation.
%s is neuroinflammation factor, or scale of neuroinflammation.
% f is release histamine/time in uM/hour.  
function f = degran_ha_mc(s)
f = 3000*s;
end