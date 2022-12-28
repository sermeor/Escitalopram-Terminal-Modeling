%% Function of  the putative HA transporter in glia. 
% UNITS in uM/h. 
% b = eha
function a = VHATg(b)
k = 10; 
V = (1)*(2.5)*(90)*60;
a = (V.*b./(k + b));
end
