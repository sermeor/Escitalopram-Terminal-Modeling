%% Function of histamine methyltransferase, 
% Histamine metabolisis in cytosol.
% UNITS in uM/h. 
% b = cha

function a = VHNMT(b,sc)

k = 4.2;  % Francis 80

V = (1)*(3.5)*(.53)*100;
      
a = (V.*b./(k + b)).*sc;

