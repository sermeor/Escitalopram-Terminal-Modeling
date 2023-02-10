%% Full script with compartmental model and Serotonin Terminal Model.
close all;
clear all;

%% Parameters for computation
sc = 1;   %scaling factor for time: sc=1 for hours
time = 25; % In hours. 
sampling_rate = 100; % number of samples per hour.
time_array = linspace(0, time, time * sampling_rate + 1); 

%% Compartmental Model of Escitalopram parameters.
%Dose parameters. 
weight = 20;                            %Mouse weight in g
dose_factor = 0;                       %mg/kg of body weight. 
dose = (dose_factor*1e6)*(weight/1000) * 0.001; % In ug. 

volume_factor = 5;                                                      %ml/kg of body weight.
volume_injection = volume_factor*(weight/1000);                         % in ml.
bioavailability = 0.8;
protein_binding = 0.56;
sert_binding = 0.15;

IP_volume = 2; % in ml.
plasma_volume = 2; %in ml.
brain_volume = 0.41; % in ml.
peripheral_volume = 25; % in ml

% Volumes in mL.
v0 = IP_volume+volume_injection;
v1 = plasma_volume;
v2 = brain_volume;
v3 = peripheral_volume;

%Molecular weight of escitalopram
molecular_weight = 324.392; % g/mol, or ug/umol.


%% Serotonin Terminal Model. 
[T,Y] = ode15s(@msc, time_array, [0.0995 0.9005 20.1651 1.6110 0.0376 67.4267 0.0601 1.5846 113.4289 0.0000 0.8639 1.0068 0.9757 0.0000 0.6947 12.6860 2.9428 1 1 dose*bioavailability 0 0 0],[], protein_binding, sert_binding, molecular_weight, v2, sc);

%% Extracting and calculating parameters. 
%Getting the ssri array in concentration.
ssri_array = (Y(:,22)/v2)*1000/molecular_weight; % uM. 


%% Plotting of results
figure;
plot(T, 1000.*Y(:,7),'g','LineWidth',3);
leg1=legend('eht');
set(leg1,'FontSize',14);
xlabel('Time (h)');
ylabel('serotonin (nM)');

figure;
plot(T, fireht(T), 'r','LineWidth',3);
leg1=legend('fire');
set(leg1,'FontSize',14);
xlabel('Time (h)');
ylabel('fireht(t) (/s)');

figure;
plot(T, ssri_array*1000, '','LineWidth',3);
leg1=legend('Escit');
set(leg1,'FontSize',14);
xlabel('Time (h)');
ylabel('Escit (nM)');

figure;
plot(T, Y(:,18), 'magenta','LineWidth',3);
leg1 = legend('SERTs ratio');
set(leg1,'FontSize',14);
xlabel('Time (h)');
ylabel('SERT ratio');

figure;
plot(T, Y(:,11), 'r','LineWidth',3);
leg1 = legend('g');
set(leg1,'FontSize',14);
xlabel('Time (h)');
ylabel('g');

