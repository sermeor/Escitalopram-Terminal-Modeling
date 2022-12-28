%% Uptake 2 rate of reuptake into glial terminals. 
% b = e5ht
% UNITS IN uM and uM/h. 

function a = VUP2(b,sc)

k = .17; %Km Wightman-bunin
V = (1)*1400; % Vmax Wightman-Bunin .78*1800

a = (V.*b./(k + b)).*sc;  % (4700)
% Km is from Bunin98, see also Daws05