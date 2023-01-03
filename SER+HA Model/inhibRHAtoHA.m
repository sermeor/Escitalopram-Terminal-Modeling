%% Function of histamine inhibition of histamine firing.
% Dependent on Serotonin-activated g-protein (b) and the  equilibrium  value
% of  the  activated  G-protein (c).
% Units in uM. 
function a = inhibRHAtoHA(b, c)
min_a = 0; 
a = 1 - (2).*(b - c);    
a(a<min_a) = min_a;
end 
%Fitted