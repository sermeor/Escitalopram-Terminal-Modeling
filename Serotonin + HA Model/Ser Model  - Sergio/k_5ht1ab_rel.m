%% Autoreceptor 5-HT1B control of SERT density. 
%Likely mediated by protein kinases. 5ht1b activation increases g-coupled
%protein which downregulates PKA, which reduces the phosporilation of SERTs
%(removal) and increases the density of SERTs. 
%Sigmoid function. 

function  ratio = k_5ht1ab_rel(g)
gbasal = 0.8639; %Steady state value of serotonin.
gdiff = g - gbasal;
max_r = 2; % Max increase of ratio.
min_r = 0; %Min value of ratio.
b = 500; % Strength.
c = 1; %Parameter for f(eht_basal) = 1) 
ratio = max_r ./ (1 + c*exp(-b.*gdiff)) + min_r;
end




