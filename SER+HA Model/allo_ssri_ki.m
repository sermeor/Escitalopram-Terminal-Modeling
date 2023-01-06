%% Function that models slow down in escit dissociation from SERTs due to allosteric binding
%Units in uM.
function ki = allo_ssri_ki(ssri)
min_ki = 0.001;
b = 4; %Strength.
ki = 0.05 *exp(-b*ssri); %constant of inhibition of ESCIT in uM. doi:10.1017/S1461145706007486.
ki(ki<min_ki) = min_ki;
end