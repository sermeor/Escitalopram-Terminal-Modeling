function fitting()


function f = fitting_function(x)

% Extra parameters
y = [0.099536676200271 0.900463323799727 20.165129111390300 1.611006959181445 ...
    0.037652145123716 67.494409421869349 0.060055870432186 1.584725391415195 1.134288512515704e+02 ...
    0.863896510449192 1.006770277834519 0.975698298241236 0.0000 0.696152764831885 1.27547664778701 ...
    0.957723862365818 1 0.1 0 0 0 0 0 3.144939308740220 1.441700144015534e+02 1.394002658962972 ...
    0.756745984753609 99.736992949085135 2.500592520752401e+02 3.125740650940501e+02 0.696152764831885 ... 
    1.27547664778701 0.957723862365818 0.863896510449192 1.006770277834519 0.975698298241236];
t = 0;
protein_binding = 0.15;
sert_binding = 0.15;
v2 = 0.41;
molecular_weight = 324.392; % g/mol, or ug/umol.

q0 = y(20); %Peritoneum concentration in ug.
q1 = y(21); %Blood concentration in ug.
q2 = y(22); %Brain concentration in ug.
q3 = y(23); %Periphery concentration in ug.




%% Pharmacokinetic model
% Rates between comparments. 
k01p = 0.6;
k10p = 3;
k12p = 9.9;
k21p = 2910;
k13p = 6;
k31p = 0.6;





%% Serotonin Terminal Model
NADP = 26;  % NADP concentration in uM.
NADPH = 330;    % NADPH concentration in uM. 
a9 = 20;    %bound auto produce G5ht*                
a10 = 200;  %T5ht* reverses G5ht* to G5ht.                  
a11 = 30;   %G5ht* produces T5ht*                     
a12 = 200;  %decay of T5ht*.                
a13 = 36000;% rate eht bounding to receptors.                       
a14 = 20000;% rate of 5ht unbonding from autoreceptor.                        
g0 = 10; % total g-protein of serotonin autoreceptors.                           
t0 = 10; % total T regulary protein of serotonin autoreceptors.                           
b0 = 10; % total serotonin autoreceptors.  
a15 = 1; %diffusion from cytosol to extracell.
a16 = 1;  %diffusion from glia to extracellular space
a17 = 1;  %catabolism of hiaa
a18 = 0.001;   %0.05 UP2 multiplier
a19 = 2;  %removal of trp
a20 = 1;   %removal of pool
a21 = 40;  %eht removal rate.
a22 = x(1);  %release per action potential.
a23 = 2; %Coefficient of SERT reuptake. 
a24 = x(2); % Coefficient of VMAT transport. 
k11 = 100;  %bound auto produce Gha*.
k12 = 961.094;   %Tha* reverses Gha* to Gha.
k13 = 20;  %Gha* produces THA*.
k14 = 66.2992;  %decay of Tha*.
k15 = 1000000;  %eha binds to auto
k16 = 13123578;   %eha dissociates from auto
gh0 = 10; %total g-protein of histamine autoreceptors. 
th0 =  10;  %total T regulary protein of histamine autoreceptors. 
bh0 =  10; %total histamine autoreceptors
ssri = (q2/v2)*1000/(molecular_weight); % SSRI concentration from compartmental model in uM -> umol/L. 
%Parameters for SERT membrane, inactivity and pool transport.
k_ps = (100) .* k_5ht1ab_rel(y(10));
k_sp = (10);
k_si = (10) .*  k_ssri_reupt(ssri);
k_is = (1);
%Equation parameters. 
% y(1) = bh2
% y(2) = bh4
% y(3) = trp
% y(4) = htp
% y(5) = cht
% y(6) = vht
% y(7) = eht
% y(8) = hia
% y(9) = trppool
% y(10) = gstar
% y(11) = tstar
% y(12) = bound
% y(13) =  glialht
%  y(14) =  Gha*
%  y(15) = Tha*
%  y(16)  =  bound ha
% y(17) = SERTs_surface
% y(18) = SERT_pool
% y(19) = SERT inactive

f(1) = VAADC(y(4)) - a24.*VMAT(y(5),y(6)) + a23.*VSERT(y(7), y(17), ssri, sert_binding) - TCcatab(y(5)) -a15.*y(5);
f(2) = a24.*VMAT(y(5),y(6)) - a22.*inhibR5HTto5HT(y(10)).*inhibRHAto5HT(y(14)).*fireht(t).*y(6); 
f(3) = a22.*inhibR5HTto5HT(y(10)).*inhibRHAto5HT(y(14)).*fireht(t).*y(6) - a23.*VSERT(y(7), y(17), ssri, sert_binding) - a18.*H1ht(y(7)).*VUP2(y(7)) - a21.*y(7) - a13.*y(7).*(b0 - y(12)) + a14.*y(12) + a15.*y(5) + a16.*y(13); 
end

x0 = [1,1];

display(fitting_function(x0));

display(fsolve(@fitting_function, x0));

end