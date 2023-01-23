%% Histamine H3 inhibition of serotonin release.
%autoreceptors and depends on levels of g-coupled proteins activated (b) by
%histamine and basal levels of g-coupled histamine (c). 
%UNITS IN ratio. 
function a = inhibRHAto5HT(b, c)
min_a = 0;
%b = gha*
a = 1 - (3).*(b - c);
a(a<min_a) = min_a;
end 
