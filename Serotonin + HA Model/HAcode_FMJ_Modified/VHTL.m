%% Function of histidine transporter. 
% Histidine transport from blood to cytosol.
% UNITS in uM/h. 
% b = bht.

function a = VHTL(b,sc)

km = 1000;  % (lobster)(6.2-19 muM conrad05)

k1 = 4680;

a = (1).*k1.*(b./(km + b)).*sc;

