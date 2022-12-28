%% Function of monoamine transporter in vesicles.
% Transport of histamine from cytosol to vesicles.
% b = cha
% c = vha
function a = VMATH(b,c)
k = 24;         
V =  1*10552; %31500

a = (V.*(b./(k + b)) - 5.*c);
%merickel95   reports 24 micromolar and also 3 micromolar from previous authors. 
% erickson96  reports 200 micromolar
%travis00(caron,wightman)     reports 3.06 + or - 1 micromolar