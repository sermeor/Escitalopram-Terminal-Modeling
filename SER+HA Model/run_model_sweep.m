function run_model_sweep(z)
close all;
%% Parameters for computation
t_factor = 1; % Time factor for graphs.
time = 100/t_factor; % Time of simulation depending on t_factor.
sampling_rate = 100*t_factor; % number of samples per time factor units.
time_array = linspace(0, time, time * sampling_rate + 1);

%% Compartmental Model of Escitalopram parameters.
%Dose parameters. 
weight = 20;                            % Mouse weight in g
dose_factor = z;                        % mg/kg of body weight. 
SSRI_start_time = 0/t_factor;           % Starting time of SSRI dose in same units as t_factor.
dose = (dose_factor*1e6)*(weight/1000) * 0.001; % In ug. 

volume_factor = 5;                                                     % ml/kg of body weight.
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

%% Mast cell model of neuroinflammation. 
mc_start_time = 0.5/t_factor; %Time to start neuroinflammation effects with mast cells.
mc_switch = 0; %Switch that turns on an off all effects of mast cell presence.

%% Basal parameters. 
btrp0 = 96; %Blood tryptophan equilibrium value. 
eht_basal = 0.06; %Steady state basal concentration of serotonin.
gstar_5ht_basal = 0.8561; %Equilibrium concentration of g* serotonin.
gstar_ha_basal =  0.7484; %Equilibrium concentration of g* histamine. 
bht0 = 100; % Blood histidine equilibrium value. 

%% Model Solving. 
[T,Y] = ode45(@msc, time_array, [95.9766	0.0993	0.9007	20.166	1.6065	0.0367	63.0457	0.0593	1.6028	113.4370	0.8573	0.9931	0.9650	0.0014	0.7221	1.3593	1.0084	1.0037	0.2463	0	dose*bioavailability	0	0	0	3.1968	140.3708	1.4717	2.0490	99.7316	247.6260	309.5325	0.7221	1.3593	1.00584	0.8573	0.9931	0.9650	354.6656	177.3328	350	150	3	140	0.7205	1.3539	1.0051],[], protein_binding, sert_binding, molecular_weight, v2, SSRI_start_time, mc_switch, mc_start_time, btrp0, eht_basal, gstar_5ht_basal, gstar_ha_basal, bht0);

%% Extracting and calculating parameters. 
%Getting the ssri array in concentration.
ssri_array = (Y(:,23)/v2)*1000/molecular_weight; % uM. 

%% Plotting of results. 
figure;
plot(T.*t_factor, 1000.*Y(:,8),'g','LineWidth',3);
leg1=legend('eht');
set(leg1,'FontSize',14);
xlabel('Time');
ylabel('serotonin (nM)');

% figure;
% plot(T.*t_factor, fireht(time_array),'g','LineWidth',3);
% leg1=legend('fire');
% set(leg1,'FontSize',14);
% xlabel('Time');
% ylabel('fire');

% figure;
% plot(T.*t_factor, ssri_array*1000, '','LineWidth',3);
% leg1=legend('Escit');
% set(leg1,'FontSize',14);
% xlabel('Time');
% ylabel('Escit (nM)');

figure;
plot(T.*t_factor, Y(:,18), 'magenta','LineWidth',3);
leg1 = legend('SERTs ratio');
set(leg1,'FontSize',14);
xlabel('Time');
ylabel('SERT ratio');

% figure;
% plot(T.*t_factor, fireha(time_array), 'magenta','LineWidth',3);
% leg1 = legend('fireha');
% set(leg1,'FontSize',14);
% xlabel('Time');
% ylabel('fireha');

figure;
plot(T.*t_factor, Y(:,27), 'magenta','LineWidth',3);
leg1 = legend('eha');
set(leg1,'FontSize',14);
xlabel('Time');
ylabel('eha');

%Copy results to CSV.
csvwrite(strcat('data',num2str(dose_factor),'mgkg.csv'), horzcat(T,Y));
%csvwrite('datacontrolstimha', horzcat(T,Y));