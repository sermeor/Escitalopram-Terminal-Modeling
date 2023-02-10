%% Function of  the putative HA transporter in glia. 
% UNITS in uM/h. 
% b = eha
function a = VHATg(b)
k = 10; 
V = 13500;
a = (V.*b./(k + b));
end
