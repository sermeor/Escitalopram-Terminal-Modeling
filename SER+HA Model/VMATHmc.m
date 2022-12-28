%% Function of monoamine transporter in vesicles.
% Transport of histamine from cytosol to vesicles.
% b = cha
% c = vha
function a = VMATHmc(b, c)
k = 24;         
V =  (1)*21104; 
a = (V.*(b./(k + b)) - 5.*c);
end