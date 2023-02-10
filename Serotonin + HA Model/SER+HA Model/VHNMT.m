%% Function of histamine methyltransferase, 
% Histamine metabolisis in cytosol.
% UNITS in uM/h. 
% b = cha
function a = VHNMT(b)
k = 4.2;  % Francis 80
V = 185.5;
a = (V.*b./(k + b));
end