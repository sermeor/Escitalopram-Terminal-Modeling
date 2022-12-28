%% Inihibition of serotonin release from 5-HT1B due to serotonin.
% Dependent on Serotonin-activated g-protein (b) and the  equilibrium  value
% of  the  activated  G-protein (c).
% UNITS IN uM. 
function a = inhibR5HTto5HT(b, c)
min_a = 0;
a = 1 - (1).*(b - c);  %3.5
a(a<min_a) = min_a;
end 