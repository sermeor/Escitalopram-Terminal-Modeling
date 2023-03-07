%% Function of histamine methyltransferase, 
% Histamine metabolisis in glia.
% UNITS in uM/h. 
% b = gha.
function a = VHNMTg(b)
k = 4.2;  % Francis 80
V = 212;
     
a = (V.*b./(k + b));

end
