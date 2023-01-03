%% Function of serotonin inhibition of histamine release.
% Dependent on Serotonin-activated g-protein (b) and the  equilibrium  value
% of  the  activated  G-protein (c).
% Units in uM. 
function a = inhibR5HTtoHA(b, c)
min_a = 0;
a = 1 + (3).*(b - c); %3
a(a<min_a) = min_a;
end
