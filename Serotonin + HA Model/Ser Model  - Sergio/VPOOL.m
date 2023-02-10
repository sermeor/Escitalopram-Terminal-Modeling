
%% Rate of storage of tryptophan in the tryptophan pool from cytosolic tryptophan. 
% b = trp
% c = pool
% UNITS IN uM and uM/h. 

function a = VPOOL(b,c,sc);

k1 = 9; %to pool
k2 = .6; %from pool
 
a = (k1.*b - k2.*c).*sc;

