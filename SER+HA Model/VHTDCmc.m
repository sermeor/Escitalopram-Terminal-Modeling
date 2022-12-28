%% Function of histidine decarboxylase in mast cells. 
% Histamine synthesis in mast cells.
% UNITS in uM/h. 
% b = gha.
function a = VHTDCmc(b)
% b = cht
%  c = G*
km = 270;  %homosapiens BRENDA  (other wild type =100, Komori2012)
v = (1)*877.50;
a =  v.*(b./(km + b)); 

