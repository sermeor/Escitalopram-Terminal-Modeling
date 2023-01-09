%% Full script with compartmental model and Serotonin Terminal Model.
close all;
%clear all;

%% Parameters for computation
t_factor = 3600; % Time factor for graphs.
time = 2*3600/t_factor; % Time of simulation depending on t_factor.
sampling_rate = 1*t_factor; % number of samples per time factor units.
time_array = linspace(0, time, time * sampling_rate + 1);

%% Compartmental Model of Escitalopram parameters.
%Dose parameters. 
weight = 20;                            % Mouse weight in g
dose_factor = 10;                        %mg/kg of body weight. 
SSRI_start_time = 1*3600/t_factor;           % Starting time of SSRI dose in same units as t_factor.
dose = (dose_factor*1e6)*(weight/1000) * 0.001; % In ug. 
SSRI_repeat_time = 10000*3600/t_factor; %Time for repeat of dose. 


volume_factor = 5;                                                      % ml/kg of body weight.
volume_injection = volume_factor*(weight/1000);                         % in ml.
bioavailability = 0.8;

IP_volume = 2; % in ml.
plasma_volume = 2; %in ml.
brain_volume = 0.41; % in ml.
peripheral_volume = 25; % in ml.

% Volumes in mL.
v0 = IP_volume+volume_injection;
v1 = plasma_volume;
v2 = brain_volume;
v3 = peripheral_volume;

%Molecular weight of escitalopram
molecular_weight = 324.392; % g/mol, or ug/umol.

%% Mast cell model of neuroinflammation. 
mc_start_time = 0.5*3600/t_factor; %Time to start neuroinflammation effects with mast cells.
mc_switch = 0; %Switch that turns on an off all effects of mast cell presence.

%% Basal parameters. 
btrp0 = 96; %Blood tryptophan equilibrium value. 
eht_basal = 0.06; %Steady state basal concentration of serotonin.
gstar_5ht_basal = 0.8561; %Equilibrium concentration of g* serotonin.
gstar_ha_basal =  0.7484; %Equilibrium concentration of g* histamine. 
bht0 = 100; % Blood histidine equilibrium value. 
vht_basal = 63.0457; %Basal vesicular 5ht. 
vha_basal = 136.3639; %Basal vesicular ha.

%% Model Solving. 
[T,Y] = ode45(@msc, time_array, [95.9766 0.0994	0.9006	20.1618	1.6094	0.0373	63.0383 280.0048	0.0603	1.6824	113.4099 0.8660	1.0112	0.9791	0.0027	0.7114	1.3245	0.9874	0.2666 1.0203 	0.2297	0	0	0	0	0	3.1074	136.3639 241.9217	1.4378	2.0126	99.7316	249.3265	311.6581	0.7114	1.3245	0.9874	0.8660	1.0112	0.9791	354.6656	177.3328	350	150	3	140	0.7205	1.3539	1.0051],[], molecular_weight, v2, SSRI_start_time, SSRI_repeat_time, dose*bioavailability, mc_switch, mc_start_time, btrp0, eht_basal, gstar_5ht_basal, gstar_ha_basal, bht0, vht_basal, vha_basal);

%% Extracting and calculating parameters. 
%Getting the ssri array in concentration.
ssri_array = (Y(:,25)/v2)*1000/molecular_weight; % uM. 

%% Plotting of results. 
figure;
plot(T.*t_factor, 1000.*Y(:,9),'LineWidth',3);
leg1=legend('eht');
set(leg1,'FontSize',14);
xlabel('Time');
ylabel('serotonin (nM)');

figure;
plot(T.*t_factor, fireht(time_array, inhibR5HTto5HT(Y(:, 12), gstar_5ht_basal).*inhibRHAto5HT(Y(:, 16), gstar_ha_basal)),'g','LineWidth',3);
leg1=legend('fire');
set(leg1,'FontSize',14);
xlabel('Time');
ylabel('fire');

figure;
plot(T.*t_factor, ssri_array.*1000,'LineWidth',3);
leg1 = legend('escit');
set(leg1,'FontSize',14);
xlabel('Time');
ylabel('escit (nM)');


figure;
plot(T.*t_factor, Y(:,20), 'magenta','LineWidth',3);
leg1 = legend('SERTs ratio');
set(leg1,'FontSize',14);
xlabel('Time');
ylabel('SERT ratio');

figure;
plot(T.*t_factor, Y(:,30), 'magenta','LineWidth',3);
leg1 = legend('eha');
set(leg1,'FontSize',14);
xlabel('Time');
ylabel('eha');

%Copy results to CSV.
%csvwrite(strcat('ha_stim_control.csv'), horzcat(T,Y));
%csvwrite('datacontrolstimha', horzcat(T,Y));