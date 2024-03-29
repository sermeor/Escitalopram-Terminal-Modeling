%% Function of histamine inhibition of synthesis of histamine.
%autoreceptors and depends on levels of g-coupled proteins activated (b) by
%histamine and basal levels of g-coupled histamine (c). 
% Units in uM. 
function a = inhibsynHAtoHA(b, c)
min_a = 0;
a = 1 - (0.1).*(b - c);
a(a<min_a) = min_a;
end
