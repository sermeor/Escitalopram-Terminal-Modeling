%% Function with the system of differential equations for 5-HT+HA varicosities.
function dy=msc(t, y, molecular_weight, v2, SSRI_start_time, SSRI_repeat_time, q_inj, mc_switch, mc_start_time, btrp0, eht_basal, gstar_5ht_basal, gstar_ha_basal, bht0, vht_basal, vha_basal)
dy=zeros(49,1);

%% Pharmacokinetic model
% Rates between comparments. 
k01p = 0.6;
k10p = 3;
k12p = 9.9;
k21p = 2910;
k13p = 6;
k31p = 0.6;

%Parameters.
protein_binding = 0.56;
protein_brain_binding = 0.15;
q0 = y(23); %Peritoneum concentration in ug.
q1 = y(24); %Blood concentration in ug.
q2 = y(25); %Brain concentration in ug.
q3 = y(26); %Periphery concentration in ug. 

% Diff. equations. 
dy(23) = SSRI_inj(t, SSRI_start_time, SSRI_repeat_time, q_inj) - k01p*(q0);
dy(24) = k01p*(q0) - (k10p + k12p)*(q1*(1-protein_binding)) + k21p*(q2*(1-protein_brain_binding)) - k13p*(q1*(1-protein_binding)) + k31p*(q3);
dy(25) = k12p*(q1*(1-protein_binding)) - k21p*(q2*(1-protein_brain_binding));
dy(26) = k13p*(q1*(1-protein_binding)) - k31p*(q3);

%% Serotonin Terminal Model
NADP = 26;  % NADP concentration in uM.
NADPH = 330;    % NADPH concentration in uM.
TRPin = 157.6; %addition of tryptophan into blood (uM/h).
a1 = 5; %Strength of btrp stabilization.  
a2 = 20;    %bound 5ht to autoreceptors produce G5ht*                
a3 = 200;  %T5ht* reverses G5ht* to G5ht.                  
a4 = 30;   %G5ht* produces T5ht*                     
a5 = 200;  %decay of T5ht*.                
a6 = 36;% rate eht bounding to autoreceptors.                       
a7 = 20;% rate of 5ht unbonding from autoreceptors.   
a8 = 1; %5-HT leakage diffusion from cytosol to extracellular space.
a9 = 1;  %5-HT leakage diffusion from glia to extracellular space.
a10 = 1;  %catabolism of hiaa
a11 = 0.01;   %UP2 multiplier
a12 = 2;  %removal of trp
a13 = 1;   %removal of pool
a14 = 40;  %eht removal rate.
a15 = 1; % rate of vht_reserve moving to vht.
a16 = 1.89; %Factor of release per firing event. 
a17 = 100;  %bound eha to heteroreceptors produce Gha*.
a18 = 961.094;   %Tha* reverses Gha* to Gha.
a19 = 20;  %Gha* produces THA*.
a20 = 66.2992;  %decay of Tha*.
a21 = 5;  %eha binds to heteroreceptors. 
a22 = 65.6179;   %eha dissociates from heteroreceptors.
g0 = 10; % total g-protein of serotonin autoreceptors.                           
t0 = 10; % total T protein of serotonin autoreceptors.                           
b0 = 10; % total serotonin autoreceptors.
gh0 = 10; %total g-protein of histamine heteroreceptors. 
th0 =  10;  %total T regulary protein of histamine heteroreceptors. 
bh0 =  10; %total histamine heteroreceptors
ssri = (q2/v2)*1000/(molecular_weight); % SSRI concentration from compartmental model in uM -> umol/L. 

%Parameters for SERT membrane, inactivity and pool transport.
k_ps = 10 .* k_5ht1ab_rel_ps(y(12), gstar_5ht_basal);
k_sp = 10 .* k_5ht1ab_rel_sp(y(12), gstar_5ht_basal);
k_si = 7.5 .*  k_ssri_reupt(ssri);
k_is = 0.75;

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
% y(16) =  Gha*
% y(17) = Tha*
% y(18) = bound ha
% y(19) = SERT_surface_phospho
% y(20) = SERTs_surface
% y(21) = SERT_pool
% y(22) = SERT inactive
dy(1) = TRPin - VTRPin(y(1)) - a1.*(y(1) - btrp0); 
dy(2) = inhibsyn5HTto5HT(y(12), gstar_5ht_basal).*VTPH(y(4),y(3)) - VDRR(y(2), NADPH, y(3), NADP); 
dy(3) = VDRR(y(2),NADPH,y(3),NADP) - inhibsyn5HTto5HT(y(12), gstar_5ht_basal).*VTPH(y(4),y(3));
dy(4) = VTRPin(y(1)) - inhibsyn5HTto5HT(y(12), gstar_5ht_basal).*VTPH(y(4),y(3)) - VPOOL(y(4),y(11)) - a12*y(4);
dy(5) = inhibsyn5HTto5HT(y(12), gstar_5ht_basal).*VTPH(y(4),y(3)) - VAADC(y(5));
dy(6) = VAADC(y(5)) - VMAT(y(6), y(7)) - VMAT(y(6),y(8)) + VSERT(y(9), y(20), ssri, allo_ssri_ki(ssri)) - TCcatab(y(6)) - a15.*(y(6) - y(9));
dy(7) = VMAT(y(6),y(7)) - a16.*fireht(t, inhibR5HTto5HT(y(12), gstar_5ht_basal).*inhibRHAto5HT(y(16), gstar_ha_basal)).*y(7) + vht_trafficking(y(7), vht_basal);
dy(8) = VMAT(y(6),y(8)) - a15.*vht_trafficking(y(7), vht_basal);
dy(9) = a16.*fireht(t, inhibR5HTto5HT(y(12), gstar_5ht_basal).*inhibRHAto5HT(y(16), gstar_ha_basal)).*y(7) - VSERT(y(9), y(20), ssri, allo_ssri_ki(ssri)) - a11.*H1ht(y(9), eht_basal).*VUP2(y(9)) - a14.*y(9) + a8.*(y(6) - y(9)) + a9.*(y(15) - y(9));
dy(10) = TCcatab(y(6)) +  TCcatab(y(15)) - a10.*y(10);
dy(11) = VPOOL(y(4),y(11)) - a13.*y(11);
dy(12)  = (a2.*y(14).^2.*(g0 - y(12)) - a3.*y(13).*y(12));
dy(13) = (a4.*y(12).^2.*(t0 - y(13))  - a5.*y(13));
dy(14) = (a6.*y(9).*(b0 - y(14))  - a7.*y(14));
dy(15) = a11*H1ht(y(9), eht_basal).*VUP2(y(9))  - TCcatab(y(15)) - a9.*(y(15) - y(9));
dy(16) = (a17.*y(18).^2.*(gh0 - y(16)) - a18.*y(17).*y(16));
dy(17) = (a19.*y(16).^2.*(th0 - y(17))  - a20.*y(17));
dy(18) = (a21.*y(30).*(bh0 - y(18))  - a22.*y(18));
dy(19) = k_ps .* y(21) - k_sp .*y(19);
dy(20) = dy(19) - k_si .* y(20) + k_is .* y(22);
dy(21) = k_sp .* y(19)  - k_ps .* y(21);
dy(22) = k_si .* y(20)  - k_is .* y(22);

%% Histamine Terminal Model. 
b1 = 15;  %HA leakage from the cytosol to the extracellular space. 
b2 = 3.5; %HA release per action potential. 
b3 = 15; %HA leakage from glia to the extracellular space.
b4 = 0.05; %HA removal from the extracellular space
b5 = 0.25;  %Strength of stabilization of blood HT near 100μM. 
b6 = 2.5; %From cHT to HTpool.
b7 = 1; %From HTpool to cHT. 
b8 = 1; %Other uses of HT remove HT. 
b9 = 1;  %From gHT to gHTpool. 
b10 = 1; %From gHTpool to gHT.  
b11 = 1; %Removal of gHT or use somewhere else. 
b12 = 10; %Factor of activation of glia histamine production. 
b13 = 100; %Bound eha to autoreceptors produce G∗. 
b14 = 961.094; %T∗ facilitates the reversion of G∗ to G. 
b15 = 20; %G∗ produces T∗. 
b16 = 66.2992; %decay coefficient of T∗
b17 = 5;  %eha binds to autoreceptors. 
b18 = 65.6179; %eha dissociates from autoreceptors.
b19 = 20;  %bound e5ht to heteroreceptors to produce G5ht* 
b20 = 200;  %T5ht* reverses G5ht* to G5ht.
b21 = 30;   %G5ht* produces T5ht* 
b22 =  200;  %decay of T5ht*.
b23 =  36;  %rate eht bounding to heteroreceptors.
b24 =  20;  %rate of 5ht unbonding from heteroreceptors. 
g0HH = 10;  %Total gstar for H3 on HA neuron
t0HH = 10; %Total tstar for H3 autoreceptors on HA neuron
b0HH = 10;  %Total H3 autoreceptors on HA neuron
g05ht = 10; % total g-protein of serotonin heteroreceptors in histamine varicosity.                           
t05ht = 10; % total T protein serotonin heteroreceptors in histamine varicosity.                           
b05ht = 10; % total serotonin heteroreceptors in histamine varicosities.
HTin = 636.5570; % Histidine input to blood histidine uM/h. 
% y(27) = cha
% y(28) = vha 
% y(29) = vha_reserve
% y(30) = eha
% y(31) = gha 
% y(32) = bht 
% y(33) = cht 
% y(34) = chtpool 
% y(35) = gstar 
% y(36) = tstar 
% y(37) = bound 
% y(38) = gstar5ht 
% y(39) = tstar5ht 
% y(40) = bound5ht 
% y(41) = ght
% y(42) = ghtpool

dy(27) = inhibsynHAtoHA(y(35), gstar_ha_basal) .* VHTDC(y(33))  - VMATH(y(27),y(28))  -  VHNMT(y(27)) - b1*(y(27) - y(30)) + VHAT(y(30)) - VMATH(y(27), y(29));
dy(28) = VMATH(y(27),y(28)) - fireha(t, inhibRHAtoHA(y(35), gstar_ha_basal).*activR5HTtoHA(y(38), gstar_5ht_basal)).*b2.*y(28) + vha_trafficking(y(28), vha_basal);
dy(29) = VMATH(y(27), y(29)) - vha_trafficking(y(28), vha_basal); 
dy(30) = fireha(t, inhibRHAtoHA(y(35), gstar_ha_basal).*activR5HTtoHA(y(38), gstar_5ht_basal)).*b2.*y(28) - VHAT(y(30)) + b3.*(y(31) - y(30)) + b1.*(y(27) - y(30)) - H1ha(y(30)).*VHATg(y(30)) - b4.*y(30) - mc_activation(t, mc_switch, mc_start_time) .* VHATmc(y(30)) + inhibRHAtoHA(y(47), gstar_ha_basal).*degran_ha_mc(mc_activation(t, mc_switch, mc_start_time));
dy(31) = H1ha(y(30)).*VHATg(y(30)) - b3.*(y(31) - y(30)) - VHNMTg(y(31)) + (1 + b12*mc_activation(t, mc_switch, mc_start_time))*VHTDCg(y(41));
dy(32) = HTin - VHTL(y(32)) - VHTLg(y(32)) - b5.*(y(32) - bht0) - mc_activation(t, mc_switch, mc_start_time).*VHTLmc(y(32)); 
dy(33) = VHTL(y(32)) - inhibsynHAtoHA(y(35), gstar_ha_basal) .* VHTDC(y(33)) - b6.*y(33) + b7.*y(34);
dy(34) = (b6.*y(33) - b7.*y(34) - b8.*y(34));
dy(35)  = (b13.*y(37).^2.*(g0HH - y(35)) - b14.*y(36).*y(35));
dy(36) = (b15.*y(35).^2.*(t0HH - y(36))  - b16.*y(36));
dy(37) = (b17.*y(30).*(b0HH - y(37))  - b18.*y(37));
dy(38) = (b19.*y(40).^2.*(g05ht - y(38)) - b20.*y(39).*y(38));
dy(39) = (b21.*y(38).^2.*(t05ht - y(39))  - b22.*y(39));
dy(40) = (b23.*y(9).*(b05ht - y(40))  - b24.*y(40));
dy(41) = VHTLg(y(32)) - (1 + b12*mc_activation(t, mc_switch, mc_start_time))*VHTDCg(y(41)) - b9.*y(41) + b10.*y(42);
dy(42) = b9.*y(41) - b10.*y(42) - b11 .* y(42);

%% Mast Cell Model
% y(43) = cht. 
% y(44) = chtpool.
% y(45) = cha. 
% y(46) = vha. 
% y(47) =  Gha*.
% y(48) = Tha*.
% y(49)  =  bound ha.

c1 = 1; %From cHT to HTpool.
c2 = 1; %From HTpool to cHT. 
c3 = 1; %Removal of cHT or use somewhere else. 
c4 = 100; % Bound autoreceptors produce g*. 
c5 = 961.094; %T∗ facilitates the reversion of G∗ to G.
c6 = 20; %G∗ produces T∗.
c7 = 66.2992; %decay coefficient of T∗
c8 = 5;  %eHA binds to autoreceptors. 
c9 = 65.6179; %eHA dissociates from autoreceptors
g0Hmc = 10;  %Total gstar for H3 on mast cell.
t0Hmc = 10; %Total tstar for H3 on mast cell.
b0Hmc = 10;  %Total H3 receptors on mast cell.

dy(43) = mc_activation(t, mc_switch, mc_start_time).*VHTLmc(y(32)) - inhibsynHAtoHA(y(47), gstar_ha_basal).*VHTDCmc(y(43)) - c1.*(y(43)) + c2.*(y(44));
dy(44) = c1.*(y(43)) - c2.*(y(44)) - c3.*(y(44));
dy(45) = inhibsynHAtoHA(y(47), gstar_ha_basal).*VHTDCmc(y(43)) - VMATHmc(y(45), y(46)) - VHNMTmc(y(45)) + mc_activation(t, mc_switch, mc_start_time) .* VHATmc(y(30));
dy(46) = VMATHmc(y(45), y(46)) - inhibRHAtoHA(y(47), gstar_ha_basal).*degran_ha_mc(mc_activation(t, mc_switch, mc_start_time));
dy(47) = c4.*y(49).^2.*(g0Hmc - y(47)) - c5.*y(48).*y(47);
dy(48) = (c6.*y(47).^2.*(t0Hmc - y(48))  - c7.*y(48));
dy(49) = (c8.*y(30).*(b0Hmc - y(49)) - c9.*y(49));
end