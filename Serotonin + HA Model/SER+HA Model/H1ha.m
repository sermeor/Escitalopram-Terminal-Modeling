%% Factor of glia histamine reuptake.
% Units in uM. 
% a = eha. 
function f=H1ha(a)
if a < 19
      f = 0.025; 
elseif  a < 29
     f = (0.1)*(a-19);   
else
    f = 1;
end