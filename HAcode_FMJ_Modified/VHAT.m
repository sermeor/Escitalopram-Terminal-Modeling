%% Function of  the putative HA transporter in HA neuron. 
% UNITS in uM/h. 

% b = eha

function a = VHAT(b,sc)

k = 10;  
          
V = 1.35*6116;   % 6513;            

a = (.5).*V.*(b./(k + b)).*sc;

