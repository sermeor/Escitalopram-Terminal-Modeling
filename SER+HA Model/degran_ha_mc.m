%% Function of mast cell degranulation.
%s is neuroinflammation factor, or scale of neuroinflammation.
% f is release rate in hour-1.  
function f = degran_ha_mc(s)
f = 3.*s;
end