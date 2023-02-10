function dy=msc(t,y,sc)

dy=zeros(10,1);

b1 = 12;  %HA leakage from the cytosol to the extracellular space. 
b2 = 3.5;%5 %HA release per action potential. 
b3 = 12; %HA leakage from glia to the extracellular space.
b4 = 0.05; %0.001 %HA removal from the extracellular space
b5 = .25;  %Strength of stabilization of blood HT near 100μM. 
b6 = 2.5; %From cHT to HTpool.
b7 = 1; %From HTpool to cHT. 
b8 = 1; %Other uses of HT remove HT. 
c9 = 4.32; %Bound autoreceptors produceG∗. 
c10 = 1.296; %T∗ facilitates the reversion of G∗ to G. 
c11 = 14.4; %G∗ produces T∗. 
c12 = 25.92; %decay coefficient of T∗
c13 = 432;  %eHA binds to autoreceptors. 
c14 = 1440; %eHA dissociates from autoreceptors
g0HH = 1;  %Total gstar for H3 on HA neuron
t0HH = 60; %Total tstar for H3 on HA neuron
b0HH = 10;  %Total H3 receptors on HA neuron

% Initial parameters. 
bht0 = 100; 

% y(1) = cha
% y(2) = vha
% y(3) = eha   
% y(4) = gha
% y(5) = bht
% y(6) = cht
% y(7) = chtpool
% y(8) = gstar
% y(9) = tstar
% y(10) = bound

dy(1) = inhibsynHA(y(8)) .* VHTDC(y(6),sc)  - VMATH(y(1),y(2),sc) -  VHNMT(y(1),sc) - b1*y(1).*sc + VHAT(y(3),sc);  
dy(2) = VMATH(y(1),y(2),sc) - inhibRHA(y(8)).*fireha(t).*b2.*y(2).*sc;
dy(3) = inhibRHA(y(8)).*fireha(t).*b2.*y(2).*sc - VHAT(y(3),sc) + b3.*y(4).*sc + b1.*y(1).*sc - H1ha(y(3)).*VHATg(y(3),sc) - b4.*y(3).*sc;
dy(4) = H1ha(y(3)).*VHATg(y(3),sc) - b3.*y(4).*sc - VHNMTg(y(4),sc);
dy(5) = HTin(t).*sc - VHTL(y(5),sc) - b5.*(y(5) - bht0).*sc;  
dy(6) = VHTL(y(5),sc) - inhibsynHA(y(8)) .* VHTDC(y(6),sc) - b6.*y(6).*sc + b7.*y(7).*sc;
dy(7) = (b6.*y(6) - b7.*y(7) - b8.*y(7)).*sc;
dy(8)  = (c9.*y(10).^2.*(g0HH - y(8)) - c10.*y(9).*y(8)).*sc;
dy(9) = (c11.*y(8).^2.*(t0HH - y(9))  - c12.*y(9)).*sc;
dy(10) = (c13.*y(3).*(b0HH - y(10))  - c14.*y(10)).*sc;
