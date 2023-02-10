%% Rate of production of c5ht from 5htp. 
%Follows conventional michaelis menten mechanism. 
% b = 5htp
%UNITS in uM and uM/h. 
function a = VAADC(b, sc);

km  = 160; 
vmax = (1)*400.*sc;  %Vmax (vitamin B6 is a cofactor and will affect Vmax)

a = vmax.*b./(km + b);

%4.1.1.28 AADC

%Km = 160 muM for 5HTp  Chico06 (wrong in BRENDA)

      %Km = 45 muM  for 5HTP  gilbert95.
      %Km = .88 muM for 5HTP  Roberge80
      %Km = 20 muM  for 5HTP  

