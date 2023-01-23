%% Function of  the putative HA transporter in HA neuron. 
% UNITS in uM/h. 
% b = eha
function a = VHAT(b)
k = 10;       
V = 4128.3;             
a = V.*(b./(k + b));
end
