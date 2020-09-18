%% Test file to experiment with data visualization
% This file is to experiment with the data visualization of the network
% output
%% Load the necessary data files
addpath('test_visualization_data')
output = load('G_8_4.mat');
single_run = load('9.mat');

%% parameters

% set epoch
epoch = 1;

%% Validation plots
% Network with Q = 3 and G = 9 second

% spike times
tspike = single_run.spikes{epoch};

% network and target signals
network = single_run.network{epoch};
target = single_run.target{epoch};

% reservoir spikes validation trial
figure(1)
subplot(2,1,1)
plot(tspike(:,2),tspike(:,1),'k.')
ylabel('Neuron')
xlabel('time/ ms')
ylim([0 200])
title('Network spikes')

% network output vs target function
subplot(2,1,2)
plot(network)
hold on
plot(target)
hold off

%% G-Q accuracy plot

% define G and Q (get this data by saving the par_comb)
G = [7,7,7,8,8,8,9,9,9];
Q = [1,2,3,1,2,3,1,2,3];

% network accuracy
run_acc = output.run_acc(:, epoch);

figure(2)
scatter3(Q, G, run_acc, 'filled');
xlabel('Q')
ylabel('G')
zlabel('accuracy')

%% Firing activity - CV plot
% Network with Q = 3 and G = 9 second

% define Cv and run accuracy
Cv = output.Cv{9, epoch}{1};
run_acc = output.run_acc(9, epoch);

% calculate firing rate
neurons = sort(single_run.spikes{epoch}(:,1).');
edges = 1:401;
fire = histcounts(neurons, edges);
avg_fire = fire / (T*10^-3);


figure(3)
scatter(Cv, avg_fire, 'filled')
xlabel('Cv')
ylabel('Firing rate Hz')
