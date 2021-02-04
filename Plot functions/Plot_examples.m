%% Plot examples
% this script runs multiple plot examples

%% Load the network output file
% For example: 'Win_0.5G_6Q_1Winp_1Pexc_0.mat'
file = load('Win_0.5G_6Q_1Winp_1Pexc_0.mat');
training_output = file.training_output;
scale_param = file.scale_param;

%% Network parameters
N = 2000; 
dt = 0.05;
epoch = 2;
trial = 1;

test_output = training_output(epoch).test_output;

%% Single trial network output and activity

% Output vs target
output = test_output.Z_out{trial};
target = test_output.Zx{trial};

figure
output_vs_target_plot(output, target)

% Network activity

% convert the spiketimes
tspikes = test_output.tspikes{trial};
tmax = length(test_output.Zx{trial});
spikes = convert_spiketimes(tspikes, dt, N, tmax, trial);

% rasterplot
figure
subplot(3, 1, [1 2])
samplesize = 100;
raster_plot(spikes, N, dt, samplesize)

% PSTH
subplot(3, 1, 3)
PSTH_plot(spikes, dt)

%% Network features

% prepare data
stats = training_output(epoch).test_output.stats;
feedback = training_output(epoch).weights.feedback;

% calculate mean firing rate and CV
[mean_fire_rate, mean_CV] = calc_mean_stats(stats, N);

% Firing rate vs CV
figure
firing_rate_vs_CV_plot(mean_fire_rate, mean_CV, feedback)

%% Network weights

% Weight matrices

% get the matrices
static = training_output(epoch).weights.static;
feedback = training_output(epoch).weights.feedback;
output = training_output(epoch).weights.output; 

% set the window and make the figure
window = 100;
figure
weight_matrices_plot(static, feedback, output, N, window)


%% Helper functions
function spikes = convert_spiketimes(tspikes, dt, N, tmax, trial)
% Converts the spiketimes from an array of spikecounts and spike times to a
% boolean array of the spikes over time

% filter the zeros out
tspikes = tspikes(tspikes(:, 1) ~= 0, :);

% prepare the spikes matrix
time = 0:dt:tmax;
Ntrial = length(trial);
Ntot = length(time);
spikes = zeros(Ntrial, N, Ntot);

% convert the spiketimes to the spikemat
for neuron = 1:N
    index = find(tspikes(:, 1) == neuron);
    index_times = round(tspikes(index, 2)./dt);
    
    spikes(Ntrial, neuron, index_times) = 1;
end

end

function [mean_fire_rate, mean_Cv] = calc_mean_stats(stats, N)
% Calculates the mean firing rate and CV over multiple trials for each
% neuron

mean_fire_rate = zeros(N, 1);
mean_Cv = zeros(N, 1);
% loop through the test trials 
for i = 1:length(stats)
    stat = stats{i};
    mean_fire_rate = mean_fire_rate + stat.A_t';
    mean_Cv = mean_Cv + stat.Cv;
end

mean_fire_rate = mean_fire_rate/length(stats);
mean_Cv = mean_Cv/length(stats);

end

%% Plot Functions
function output_vs_target_plot(output, target)
%Plots the output signal vs the target signal

% Plot the output vs target
hold on
plot(output, 'LineWidth', 1.5)
plot(target, 'LineWidth', 1.5)
ylim([-2.5 2.5])
xlim([-50 3500])
legend('network output', 'target')
title('Output vs target')
hold off
end

function raster_plot(spikes, N, dt, samplesize)
%Creates a reaster plot of a samplesize of random neurons in the network.
%The spikes should be given as an boolean array with the neuron spikes over time 

% select random neurons
index = randperm(N, samplesize);

hold on
title('Reservoir spikes')
for i = 1:length(index)
    nn = index(i);
    spikeindices = find(squeeze(spikes(1,nn,:))==1);
    plot(dt*spikeindices, i*ones(1,length(spikeindices)), '.k')
end
xlim([-50 3500])
hold off

end

function PSTH_plot(spikes, dt)
%Makes a PSTH plot of the entire network activity.
%The spikes should be given as an boolean array with the neuron spikes over time 

% create PSTH over network of neurons
PSTH = squeeze(sum(spikes,2));

% squeezed network PSTH for the plot
temp_signal = PSTH(1:(end-1))';
binwidth = 1/dt;
PSTH_squeeze = squeeze(sum(reshape(temp_signal, size(temp_signal, 1),binwidth,[]),2));

% plot the PSTH
p = plot(PSTH_squeeze, 'LineWidth', 1.5);
ylabel('spikes', 'FontWeight', 'bold', 'FontSize', 10)
xlabel('time in ms', 'FontWeight', 'bold', 'FontSize', 10) 
ylim([0 200])
xlim([-50 3500])

end

function firing_rate_vs_CV_plot(mean_fire_rate, mean_CV, feedback)
% Creates a scatter plot of the mean firing rate vs CV of neurons

scatter3(mean_CV, mean_fire_rate, feedback, 15, feedback, 'filled')
colormap(cool)
hcb = colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = 'Feedback';
set(colorTitleHandle ,'String',titleString, 'FontWeight', 'bold', 'FontSize', 11);
title('Average firing rate vs CV', 'FontSize', 11)
xlabel('CV', 'FontWeight', 'bold', 'FontSize', 11)
%ylim([0 3])
%xlim([0 6])
ylabel('firing rate in Hz', 'FontWeight', 'bold', 'FontSize', 11)
xlim([0 2.5])
ylim([0 35])

end

function weight_matrices_plot(static, feedback, output, N, window)
% Plots the static, feedback and learnt weight matrices as well as the
% summed synaptic input of a neuron

learned = feedback * output';
final = static + learned;

% sample neurons from the matrices
lower = 0.5*N - 0.5*window;
upper = 0.5*N + 0.5*window;
filter_static = static(lower:upper, lower:upper);
filter_learned = learned(lower:upper, lower:upper);
filter_final = final(lower:upper, lower:upper);

% plot the matrices
subplot(2,3,1)
imagesc(filter_static)
title('Static weight matrix')
xlabel('neurons', 'FontWeight', 'bold') 
ylabel('neurons', 'FontWeight', 'bold')
caxis([-4 4])

subplot(2,3,2)
imagesc(filter_learned)
title('Learnt weight matrix')
caxis([-4 4])

subplot(2,3,3)
imagesc(filter_final)
title('Learnt adjacency matrix')
caxis([-4 4])
ylabel('Neuron', 'FontWeight', 'bold', 'FontSize', 11)
xlabel('Neuron', 'FontWeight', 'bold', 'FontSize', 11)
colorbar


% plot the synaptic input over the neurons
subplot(2,3,5)
sum_synaptic_input = sum(final, 2);
%nbins = 20;
histogram(sum_synaptic_input, 'DisplayStyle', 'stairs', 'LineWidth', 2)
title('Distribution of synaptic input (without thalamus input)')
ylabel('# neurons')
xlabel('synaptic input')
end
