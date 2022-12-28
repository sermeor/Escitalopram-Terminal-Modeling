%% Function of histidine decarboxylase
% Histamine synthesis in cytosol.
% UNITS in uM/h. 
% b = gha.

function a = VHTDC(b, sc)

% b = cht
%  c = G*

km = 270;  %homosapiens BRENDA  (other wild type =100, Komori2012)

k1 = 1*234;
 
a =  k1.*(b./(km + b)).*sc; 

