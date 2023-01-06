%% Function of SSRI decrease of SERT ratio. 
% units in uM. 
function ratio = k_ssri_reupt(ssri)
max_r = 2; % Max increase of speed.
b = 2.5; % Strength.
ratio =  b * ssri;
ratio(ratio > max_r) = max_r;
end







