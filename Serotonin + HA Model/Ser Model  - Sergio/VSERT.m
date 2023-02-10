%% Rate of reuptake from serotonin transporters in serotonin terminals.
% b = e5ht
%UNITS IN uM and uM/h. 
function a = VSERT(b, sert_density, ssri, sert_binding, sc)
k = (1).*.060;  
vmax = 250;

ki = (1)*27.6/1000; %constant of inhibition of ESCIT in uM. doi:10.1017/S1461145706007486. It was 1.1 nM in vitro.

vmax_app = (1)*vmax*sert_density;
k_app = k * (1 + sert_binding.*ssri/ki);

a = (1)*vmax_app.*b./(k_app + b).*sc; 
end

