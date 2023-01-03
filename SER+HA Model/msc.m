%% Function with the system of differential equations for 5-HT+HA varicosities.
function dy=msc(t, y, molecular_weight, v2, SSRI_start_time, mc_switch, mc_start_time, btrp0, eht_basal, gstar_5ht_basal, gstar_ha_basal, bht0,vht_basal)
dy=zeros(48,1);

%% Pharmacokinetic model
% Rates between comparments. 
k01p = (1)*0.6;
k10p = (1)*3;
k12p = (1)*9.9;
k21p = (1)*2910;
k13p = (1)*6;
k31p = (1)*0.6;


%Parameters.
protein_binding = 0.56;
sert_binding = 0.15;


q0 = y(23).*start_SSRI(t, SSRI_start_time); %Peritoneum concentration in ug.
q1 = y(24); %Blood concentration in ug.
q2 = y(25); %Brain concentration in ug.
q3 = y(26); %Periphery concentration in ug. 

% Diff. equations. 
dy(23) = - k01p*(q0);
dy(24) = k01p*(q0) - (k10p + k12p)*(q1*(1-protein_binding)) + k21p*(q2*(1-sert_binding)) - k13p*(q1*(1-protein_binding)) + k31p*(q3);
dy(25) = k12p*(q1*(1-protein_binding)) - k21p*(q2*(1-sert_binding));
dy(26) = k13p*(q1*(1-protein_binding)) - k31p*(q3);

% End of pharmacokinetic model

%% Serotonin Terminal Model
NADP = 26;  % NADP concentration in uM.
NADPH = 330;    % NADPH concentration in uM. 
a8 = 5; %Strength of btrp stabilization.  
a9 = 20;    %bound auto produce G5ht*                
a10 = 200;  %T5ht* reverses G5ht* to G5ht.                  
a11 = 30;   %G5ht* produces T5ht*                     
a12 = 200;  %decay of T5ht*.                
a13 = (1)*36;% rate eht bounding to receptors.                       
a14 = (1)*20;% rate of 5ht unbonding from autoreceptor. 
g0 = 10; % total g-protein of serotonin autoreceptors.                           
t0 = 10; % total Twas your regulary protein of serotonin autoreceptors.                           
b0 = 10; % total serotonin autoreceptors.  
a15 = 1; %5-HT leakage diffusion from cytosol to extracellular space.
a16 = 1;  %5-HT leakage diffusion from glia to extracellular space.
a17 = 1;  %catabolism of hiaa
a18 = 0.01;   %UP2 multiplier
a19 = 2;  %removal of trp
a20 = 1;   %removal of pool
a21 = 40;  %eht removal rate.
a22 = 1; % rate of vht_reserve moving to vht.
a23 = (1)*1.89; %Release constant. 
a24 = (1)*1; % VMAT constant.  
a25 = (1)*1; % VSERT constant.
k11 = 100;  %bound auto produce Gha*.
k12 = 961.094;   %Tha* reverses Gha* to Gha.
k13 = 20;  %Gha* produces THA*.
k14 = 66.2992;  %decay of Tha*.
k15 = 5;  %eha binds to auto
k16 = 65.6179;   %eha dissociates from auto
gh0 = 10; %total g-protein of histamine autoreceptors. 
th0 =  10;  %total T regulary protein of histamine autoreceptors. 
bh0 =  10; %total histamine autoreceptors
ssri = (q2/v2)*1000/(molecular_weight); % SSRI concentration from compartmental model in uM -> umol/L. 

%Parameters for SERT membrane, inactivity and pool transport.
k_ps = 10 .* k_5ht1ab_rel_ps(y(12), gstar_5ht_basal);
k_sp = 10 .* k_5ht1ab_rel_sp(y(12), gstar_5ht_basal);
k_si = (1)*7.5 .*  k_ssri_reupt(ssri);
k_is = (1)*0.75;


%Equation parameters. 
% y(1) = btrp
% y(2) = bh2
% y(3) = bh4
% y(4) = trp
% y(5) = htp
% y(6) = cht
% y(7) = vht
% y(8) = vht_reserve
% y(9) = eht
% y(10) = hia
% y(11) = trppool
% y(12) = gstar
% y(13) = tstar
% y(14) = bound
% y(15) =  glialht
%  y(16) =  Gha*
%  y(17) = Tha*
%  y(18) = bound ha
% y(20) = SERTs_surface
% y(19) = SERT_surface_phospho
% y(21) = SERT_pool
% y(22) = SERT inactive
dy(1) = TRPin(t) - VTRPin(y(1)) - a8.*(y(1) - btrp0); 
dy(2) = inhibsyn5HTto5HT(y(12), gstar_5ht_basal).*VTPH(y(4),y(3)) - VDRR(y(2), NADPH, y(3), NADP); 
dy(3) = VDRR(y(2),NADPH,y(3),NADP) - inhibsyn5HTto5HT(y(12), gstar_5ht_basal).*VTPH(y(4),y(3));
dy(4) = VTRPin(y(1)) - inhibsyn5HTto5HT(y(12), gstar_5ht_basal).*VTPH(y(4),y(3)) - VPOOL(y(4),y(11)) - a19*y(4);
dy(5) = inhibsyn5HTto5HT(y(12), gstar_5ht_basal).*VTPH(y(4),y(3)) - VAADC(y(5));
dy(6) = VAADC(y(5)) - a24.*VMAT(y(6),y(7)) + a25.*VSERT(y(9), y(20), ssri) - TCcatab(y(6)) - a15.*(y(6) - y(9));
dy(7) = a24.*VMAT(y(6),y(7)) - a23.*fireht(t, inhibR5HTto5HT(y(12), gstar_5ht_basal).*inhibRHAto5HT(y(16), gstar_ha_basal)).*y(7) + vht_trafficking(y(7), vht_basal);
dy(8) = a24.*VMAT(y(6),y(8)) - a22.*vht_trafficking(y(7), vht_basal);
dy(9) = a23.*fireht(t, inhibR5HTto5HT(y(12), gstar_5ht_basal).*inhibRHAto5HT(y(16), gstar_ha_basal)).*y(7) - a25.*VSERT(y(9), y(20), ssri) - a18.*H1ht(y(9), eht_basal).*VUP2(y(9)) - a21.*y(9) + a15.*(y(6) - y(9)) + a16.*(y(15) - y(9));
dy(10) = TCcatab(y(6)) +  TCcatab(y(15)) - a17.*y(10);
dy(11) = VPOOL(y(4),y(11)) - a20.*y(11);
dy(12)  = (a9.*y(14).^2.*(g0 - y(12)) - a10.*y(13).*y(12));
dy(13) = (a11.*y(12).^2.*(t0 - y(13))  - a12.*y(13));
dy(14) = (a13.*y(9).*(b0 - y(14))  - a14.*y(14));
dy(15) = a18*H1ht(y(9), eht_basal).*VUP2(y(9))  - TCcatab(y(15)) - a16.*(y(15) - y(9));
dy(16) = (k11.*y(18).^2.*(gh0 - y(16)) - k12.*y(17).*y(16));
dy(17) = (k13.*y(16).^2.*(th0 - y(17))  - k14.*y(17));
dy(18) = (k15.*y(29).*(bh0 - y(18))  - k16.*y(18));
dy(19) = k_ps .* y(21) - k_sp .*y(19);
dy(20) = dy(19) - k_si .* y(20) + k_is .* y(22);
dy(21) = k_sp .* y(19)  - k_ps .* y(21);
dy(22) = k_si .* y(20)  - k_is .* y(22);

%% Histamine Terminal Model. 
b1 = 15;  %HA leakage from the cytosol to the extracellular space. 
b2 = 3.5;%5 %HA release per action potential. 
b3 = 15; %HA leakage from glia to the extracellular space.
b4 = 0.05; %0.001 %HA removal from the extracellular space
b5 = 0.25;  %Strength of stabilization of blood HT near 100μM. 
b6 = 2.5; %From cHT to HTpool.
b7 = 1; %From HTpool to cHT. 
b8 = 1; %Other uses of HT remove HT. 
b9 = 1;  %From gHT to gHTpool. 
b10 = 1; %From gHTpool to gHT.  
b11 = 1; %Removal of gHT or use somewhere else. 
b12 = 200; %Factor of mast cell activation of glia histamine production. 
c9 = 100; %Bound autoreceptors produceG∗. 
c10 = 961.094; %T∗ facilitates the reversion of G∗ to G. 
c11 = 20; %G∗ produces T∗. 
c12 = 66.2992; %decay coefficient of T∗
c13 = 5;  %eHA binds to autoreceptors. 
c14 = 65.6179; %eHA dissociates from autoreceptors
d1 = 20;  %bound auto produce G5ht* 
d2 = 200;  %T5ht* reverses G5ht* to G5ht.
d3 = 30;   %G5ht* produces T5ht* 
d4 =  200;  %decay of T5ht*.
d5 =  36;  %rate eht bounding to receptors.
d6 =  20;  %rate of 5ht unbonding from autoreceptor. 
g0HH = 10;  %Total gstar for H3 on HA neuron
t0HH = 10; %Total tstar for H3 on HA neuron
b0HH = 10;  %Total H3 receptors on HA neuron
g05ht = 10; % total g-protein of serotonin receptors in histamine varicosity.                           
t05ht = 10; % total T regulary protein of serotonin receptors in histamine varicosity.                           
b05ht = 10; % total serotonin receptors in histamine varicosities.

% y(27) = cha
% y(28) = vha 
% y(29) = eha
% y(30) = gha 
% y(31) = bht 
% y(32) = cht 
% y(33) = chtpool 
% y(34) = gstar 
% y(35) = tstar 
% y(36) = bound 
% y(37) = gstar5ht 
% y(38) = tstar5ht 
% y(39) = bound5ht 
% y(40) = ght
% y(41) = ghtpool

dy(27) = inhibsynHAtoHA(y(34), gstar_ha_basal) .* VHTDC(y(32))  - VMATH(y(27),y(28)) -  VHNMT(y(27)) - b1*(y(27) - y(29)) + VHAT(y(29));
dy(28) = VMATH(y(27),y(28)) - fireha(t, inhibRHAtoHA(y(34), gstar_ha_basal).*inhibR5HTtoHA(y(37), gstar_5ht_basal)).*b2.*y(28);
dy(29) = fireha(t, inhibRHAtoHA(y(34), gstar_ha_basal).*inhibR5HTtoHA(y(37), gstar_5ht_basal)).*b2.*y(28) - VHAT(y(29)) + b3.*(y(30) - y(29)) + b1.*(y(27) - y(29)) - H1ha(y(29)).*VHATg(y(29)) - b4.*y(29) - mc_activation(t, mc_switch, mc_start_time) .* VHATmc(y(29)) + inhibRHAtoHA(y(46), gstar_ha_basal).*degran_ha_mc(mc_activation(t, mc_switch, mc_start_time));
dy(30) = H1ha(y(29)).*VHATg(y(29)) - b3.*(y(30) - y(29)) - VHNMTg(y(30)) + (1 + b12*mc_activation(t, mc_switch, mc_start_time))*VHTDCg(y(40));
dy(31) = HTin(t) - VHTL(y(31)) - VHTLg(y(31)) - b5.*(y(31) - bht0) - mc_activation(t, mc_switch, mc_start_time).*VHTLmc(y(31)); 
dy(32) = VHTL(y(31)) - inhibsynHAtoHA(y(34), gstar_ha_basal) .* VHTDC(y(32)) - b6.*y(32) + b7.*y(33);
dy(33) = (b6.*y(32) - b7.*y(33) - b8.*y(33));
dy(34)  = (c9.*y(36).^2.*(g0HH - y(34)) - c10.*y(35).*y(34));
dy(35) = (c11.*y(34).^2.*(t0HH - y(35))  - c12.*y(35));
dy(36) = (c13.*y(29).*(b0HH - y(36))  - c14.*y(36));
dy(37) = (d1.*y(39).^2.*(g05ht - y(37)) - d2.*y(38).*y(37));
dy(38) = (d3.*y(37).^2.*(t05ht - y(38))  - d4.*y(38));
dy(39) = (d5.*y(9).*(b05ht - y(39))  - d6.*y(39));
dy(40) = VHTLg(y(31)) - (1 + b12*mc_activation(t, mc_switch, mc_start_time))*VHTDCg(y(40)) - b9.*y(40) + b10.*y(41);
dy(41) = b9.*y(40) - b10.*y(41) - b11 .* y(41);



%% Mast Cell Model
% y(42) = cht. 
% y(43) = chtpool.
% y(44) = cha. 
% y(45) = vha. 
% y(46) =  Gha*.
% y(47) = Tha*.
% y(48)  =  bound ha.

e1 = 1; %From cHT to HTpool.
e2 = 1; %From HTpool to cHT. 
e3 = 1; %Removal of gHT or use somewhere else. 
e4 = 100; % Bound autoreceptors produce g*. 
e5 = 961.094; %T∗ facilitates the reversion of G∗ to G.
e6 = 20; %G∗ produces T∗.
e7 = 66.2992; %decay coefficient of T∗
e8 = 5;  %eHA binds to autoreceptors. 
e9 = 65.6179; %eHA dissociates from autoreceptors
g0Hmc = 10;  %Total gstar for H3 on mast cell.
t0Hmc = 10; %Total tstar for H3 on mast cell.
b0Hmc = 10;  %Total H3 receptors on mast cell.

dy(42) = mc_activation(t, mc_switch, mc_start_time).*VHTLmc(y(31)) - inhibsynHAtoHA(y(46), gstar_ha_basal).*VHTDCmc(y(42)) - e1.*(y(42)) + e2.*(y(43));
dy(43) = e1.*(y(42)) - e2.*(y(43)) - e3.*(y(43));
dy(44) = inhibsynHAtoHA(y(46), gstar_ha_basal).*VHTDCmc(y(42)) - VMATHmc(y(44), y(45)) - VHNMTmc(y(44)) + mc_activation(t, mc_switch, mc_start_time) .* VHATmc(y(29));
dy(45) = VMATHmc(y(44), y(45)) - inhibRHAtoHA(y(46), gstar_ha_basal).*degran_ha_mc(mc_activation(t, mc_switch, mc_start_time));
dy(46) = e4.*y(48).^2.*(g0Hmc - y(46)) - e5.*y(47).*y(46);
dy(47) = (e6.*y(46).^2.*(t0Hmc - y(47))  - e7.*y(47));
dy(48) = (e8.*y(29).*(b0Hmc - y(48)) - e9.*y(48));
end