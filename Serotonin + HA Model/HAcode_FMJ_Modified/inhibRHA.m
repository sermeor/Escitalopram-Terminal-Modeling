%% Function of histamine inhibition of firing.
% Units in uM. 

function a = inhibRHA(b)
%b = gstar;

a = 1 - 3.5.*(b - .6945);     %17 is the magic number for control .9282

%a = 1;
end 