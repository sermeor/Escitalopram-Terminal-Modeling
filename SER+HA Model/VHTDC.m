%% Function of histidine decarboxylase
% Histamine synthesis in cytosol.
% UNITS in uM/h. 
function a = VHTDC(b)
% b = cht
%  c = G*
km = 270;  %homosapiens BRENDA  (other wild type =100, Komori2012)
v = 1*234;
a =  v.*(b./(km + b)); 

