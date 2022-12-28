%% 5HT catabolism in the cytosol
%Conventional Michaelis Menten rate.
% b = cytosolic 5ht in the terminal, sc is the scaling factor for time. 
% UNITS IN uM and uM/h. 
function a = TCcatab(b)
km = 95;  %Gottowik93 and Fowler94
vmax = 4000;
a = (vmax.*b./(km + b));
end