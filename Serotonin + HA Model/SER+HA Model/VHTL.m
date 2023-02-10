%% Function of histidine transporter. 
% Histidine transport from blood to cytosol.
% UNITS in uM/h. 
% b = bht.
function a = VHTL(b)
km = 1000;  % (lobster)(6.2-19 muM conrad05)
k1 = 4680;
a = k1.*(b./(km + b));

end