function dy=msc(t, y, protein_binding, sert_binding, molecular_weight, v2, sc)

dy=zeros(23,1);

%% Pharmacokinetic model
% Rates between comparments. 
k01p = 0.6*sc;
k10p = 3*sc;
k12p = 9.9*sc;
k21p = 2910*sc;
k13p = 6*sc;
k31p = 0.6*sc;

q0 = y(20); %Peritoneum concentration in ug.
q1 = y(21); %Blood concentration in ug.
q2 = y(22); %Brain concentration in ug.
q3 = y(23); %Periphery concentration in ug. 

% Diff. equations. 
dy(20) = -k01p*(q0);
dy(21) = k01p*(q0) - (k10p + k12p)*(q1*(1-protein_binding)) + k21p*(q2*(1-sert_binding)) - k13p*(q1*(1-protein_binding)) + k31p*(q3);
dy(22) = k12p*(q1*(1-protein_binding)) - k21p*(q2*(1-sert_binding));
dy(23) = k13p*(q1*(1-protein_binding)) - k31p*(q3);

% End of pharmacokinetic model

%% Serotonin Terminal Model
NADP = 26;  % NADP concentration in uM.
NADPH = 330;    % NADPH concentration in uM. 
a9 = 20;   %20  %bound auto produce G5ht*                
a10 = 200; %200 %T5ht* reverses G5ht* to G5ht.                  
a11 = 30;  %30  %G5ht* produces T5ht*                     
a12 = 200; %200 %decay of T5ht*.                
a13 = 36000; %36000 % rate eht bounding to receptors.                       
a14 = 20000; %20000 % rate of 5ht unbonding from autoreceptor.                        
g0 = 10; % total g-protein of serotonin autoreceptors.                           
t0 = 10; % total T regulary protein of serotonin autoreceptors.                           
b0 = 10; % total serotonin autoreceptors.  
a15 = 1; %diffusion from cytosol to extracell.
a16 = 1;  %diffusion from glia to extracellular space
a17 = 1;  %catabolism of hiaa
a18 = 0.05;   %UP2 multiplier
a19 = 2;  %removal of trp
a20 = 1;   %removal of pool
a21 = 40; %40 % eht removal rate. 
a22 = 1;  %release per action potential.
k11 = 4.32;  %bound auto produce Gha*
k12 = 1.296;  %Tha* reverses Gha* to Gha.
k13 = 14.4;   %Gha* produces THA*
k14 =  25.92;  %decay of Tha*
k15  =  432;  %eha binds to auto
k16 =  1440;  %eha dissociates from auto
gh0 = 1; %total g-protein of histamine autoreceptors. 
th0 =  60;  %total T regulary protein of histamine autoreceptors. 
bh0 =  10; %total histamine autoreceptors
eha = 1.39; % Histamine concentration in uM. 
ssri = (q2/v2)*1000/(molecular_weight); % SSRI concentration from compartmental model in uM -> umol/L. 
%Parameters for SERT membrane and pool transport.
k_ps = (1) .* k_5ht1ab_rel(y(11));
k_sp = (1) .*  k_ssri_reupt(ssri);

dy(1) = inhibsyn5HT(y(11),1).*VTPH(y(3),y(2),sc) - VDRR(y(1),NADPH,y(2),NADP,sc); 
dy(2) = VDRR(y(1),NADPH,y(2),NADP,sc) - inhibsyn5HT(y(11),1).*VTPH(y(3),y(2),sc);
dy(3) = VTRPin(btrp(t),sc) - inhibsyn5HT(y(11),1).*VTPH(y(3),y(2),sc) - VPOOL(y(3),y(9),sc) - a19*y(3).*sc;
dy(4) = inhibsyn5HT(y(11),1).*VTPH(y(3),y(2),sc) - VAADC(y(4),sc);
dy(5) = VAADC(y(4),sc) - VMAT(y(5),y(6),sc) + VSERT(y(7), y(18), ssri, sert_binding, sc) - TCcatab(y(5),sc) -a15.*y(5).*sc;
dy(6) = VMAT(y(5),y(6),sc) - a22.*inhibR5HT(y(11),1).*inhibHA(y(15)).*fireht(t).*y(6).*sc; 
dy(7) = a22.*inhibR5HT(y(11),1).*inhibHA(y(15)).*fireht(t).*y(6).*sc - VSERT(y(7), y(18), ssri, sert_binding, sc) - a18*H1ht(y(7)).*VUP2(y(7),sc) - a21.*y(7).*sc - (1).*a13.*y(7).*(b0 - y(13)).*sc + (1).*a14.*y(13).*sc + a15.*y(5).*sc + a16.*y(14).*sc; 
dy(8) = TCcatab(y(5),sc) +  TCcatab(y(14),sc) - a17.*y(8).*sc;
dy(9) = VPOOL(y(3),y(9),sc) - a20.*y(9).*sc;
dy(10) = sin(t);
dy(11)  = (a9.*y(13).^2.*(g0 - y(11)) - a10.*y(12).*y(11)).*sc;
dy(12) = (a11.*y(11).^2.*(t0 - y(12))  - a12.*y(12)).*sc;
dy(13) = (a13.*y(7).*(b0 - y(13)).*sc  - a14.*y(13)).*sc;
dy(14) = a18*H1ht(y(7)).*VUP2(y(7),sc)  - TCcatab(y(14),sc) -a16.*y(14).*sc;
dy(15) = (k11.*y(17).^2.*(gh0 - y(15)) - k12.*y(16).*y(15)).*sc;
dy(16) = (k13.*y(15).^2.*(th0 - y(16))  - k14.*y(16)).*sc;
dy(17) = (k15.*eha.*(bh0 - y(17))  - k16.*y(17)).*sc;
dy(18) = k_ps .* y(19) -  k_sp.*y(18);
dy(19) = k_sp.*y(18) - k_ps .* y(19);

%END OF TESTING

% y(1) = bh2
% y(2) = bh4
% y(3) = trp
% y(4) = htp
% y(5) = ht
% y(6) = vht
% y(7) = eht
% y(8) = hia
% y(9) = trppool
% y(10) = dummy
% y(11) = gstar
% y(12) = tstar
% y(13) = bound
% y(14) =  glialht
%  y(15) =  Gha*
%  y(16) = Tha*
%  y(17)  =  bha
% y(18) = SERTs_surface
% y(19) = SERT_pool


end
