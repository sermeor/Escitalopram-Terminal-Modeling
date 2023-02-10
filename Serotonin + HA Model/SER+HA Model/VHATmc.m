%% Function of  the putative HA transporter in mast cells. 
% UNITS in uM/h. 
% b = eha
function a = VHATmc(b)
k = 10; 
V =  3375;
a = (V.*b./(k + b));
end
