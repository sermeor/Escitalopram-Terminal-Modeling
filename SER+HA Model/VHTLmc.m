%% Function of histidine transporter. 
% Histidine transport from blood to mast cell.
% UNITS in uM/h. 
% b = bht.
function a = VHTLmc(b)
km = 1000;  % (lobster)(6.2-19 muM conrad05)
vmax = (1)*4680;
a = vmax.*(b./(km + b));
end