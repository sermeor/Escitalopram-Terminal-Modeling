%% Function of histidine decarboxylase in glia. 
% Histamine synthesis in glia.
% UNITS in uM/h. 
% b = gha.
function a = VHTDCg(b)
% b = cht
%  c = G*
km = 270;  %homosapiens BRENDA  (other wild type =100, Komori2012)
v = 61.4250;
a =  v.*(b./(km + b)); 

