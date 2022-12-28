%% Rate of reuptake from serotonin transporters in serotonin terminals.
% b = e5ht
%UNITS IN uM and uM/h. 
function a = VSERT(b, sert_density, ssri)
k = (1)*0.060;  
vmax = (1)*250;
ki = (1)*27.6/1000; %constant of inhibition of ESCIT in uM. doi:10.1017/S1461145706007486. 
vmax_app = (1)*vmax*sert_density;
k_app = k * (1 + ssri/ki);
a = (1)*vmax_app.*b./(k_app + b); 
end

