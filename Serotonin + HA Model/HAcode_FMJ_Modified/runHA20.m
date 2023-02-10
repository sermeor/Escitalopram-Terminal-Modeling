% HA including synthesis
close all;
clear all;

%% Parameters for computation
sc = 1;   %scaling factor for time: sc=1 for hours
time = 10; % In hours. 
sampling_rate = 100; % number of samples per hour.
time_array = linspace(0, time, time * sampling_rate + 1); 

[T,Y] = ode15s(@msc, time_array, [3.144939308740220 1.441700144015534e+02 1.399984798740298 0.756745984753609 99.736992949085135 2.500592520752401e+02 3.125740650940501e+02 0.696152764831885 12.727547664778701 2.957723862365818]  ,[],sc);


figure; plot(time_array, Y(:, 1), 'LineWidth', 3); title('cha');
figure; plot(time_array, Y(:, 2), 'LineWidth', 3); title('vha');
figure; plot(time_array, Y(:, 3), 'LineWidth', 3); title('eha');
figure; plot(time_array, Y(:, 4), 'LineWidth', 3); title('gha');

figure; plot(time_array, Y(:, 5), 'r', 'LineWidth', 3); title('bht');
figure; plot(time_array, Y(:, 6), 'r','LineWidth', 3); title('cht');
figure; plot(time_array, Y(:, 7), 'r','LineWidth', 3); title('chtpool');
figure; plot(time_array, Y(:, 8), 'r','LineWidth', 3); title('gstar');
figure; plot(time_array, Y(:, 9), 'r','LineWidth', 3); title('tstar');
figure; plot(time_array, Y(:, 10), 'r','LineWidth', 3); title('boundha');


