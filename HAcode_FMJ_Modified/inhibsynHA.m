%% Function of histamine inhibition of synthesis.
% Units in uM. 

function a = inhibsynHA(b)
a = (1 - 3.5.*(b - .6945));
end