function ratio = k_ssri_reupt(ssri)
b = 50; %Strength
c = 1.75; %Amplitude
d = 5; %Plateau size.
ratio = 1 + c.*(ssri.^(2.*d))./(((1./b).^(2.*d)) +(ssri.^(2.*d)));
end