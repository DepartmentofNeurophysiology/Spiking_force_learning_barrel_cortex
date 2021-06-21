%% Plot examples
% this script runs multiple plot examples
%   * single trial network output
%   * sequence of trials with network output 
%   * network activity
%   * network features
%   * network weights
%   * synaptic input
%   * network weight changes
% Helper functions:
%   * convert_spiketimes
%   * calc_mean_stats
%   * create_trial_seq
% Plot functions:
%   * output_vs_target_plot
%   * trial_sequence_plot
%   * raster_plot
%   * PSTH_plot
%   * firing_rate_vs_CV_plot
%   * weight_matrices_plot
%   * synaptic_input_plot
%   * weight_change_plot
%% Load the network output file
% for example: 'Win_0.5G_10Q_1Winp_1Pexc_0.mat'
file = load('Win_0.5G_3Q_1Winp_1Pexc_0.mat');
training_output = file.training_output;
scale_param = file.scale_param;

%% Network parameters
N = 2000;    % number of neurons
N_train = 600;  % number of training trials
dt = 0.05;   % integration time constant (ms)
epoch = 1;   % number of epoch (of N_total)
N_test = 50;  % validation trial number (of N_test)

test_output = training_output(epoch).test_output;

%% Single trial network output
% output vs target
output = test_output.Z_out{N_test};
target = test_output.Zx{N_test};

figure
output_vs_target_plot(output, target)
%% Sequence of trials with network output
% Select a trainings epoch
epoch = 2;

% prepare data
test_output = training_output(epoch).test_output;

% get a sample of test trials 
sample_size = 30;
trials_num = length(test_output.Z_out);
[Z_out_seq, Zx_seq] = create_trial_seq(trials_num, sample_size, test_output);

figure
trial_sequence_plot(Z_out_seq, Zx_seq);

sgtitle(['Network output and target for a sample of ' num2str(sample_size)...
    ' trials'])

%% network activity
% convert the spiketimes
tspikes = test_output.tspikes{N_test};
tmax = length(test_output.Zx{N_test});
spikes = convert_spiketimes(tspikes, dt, N, tmax, N_test);

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

% firing rate vs CV
figure
firing_rate_vs_CV_plot(mean_fire_rate, mean_CV, feedback)

%% Network weights
% weight matrices

% get the matrices
static = training_output(epoch).weights.static;
feedback = training_output(epoch).weights.feedback;
output = training_output(epoch).weights.output; 

% set the window and make the figure
window = 100;
figure
weight_matrices_plot(static, feedback, output, N, window)

%% Synaptic input

static = training_output(epoch).weights.static;
feedback = training_output(epoch).weights.feedback;
output = training_output(epoch).weights.output; 

figure 
synaptic_input_plot(static, feedback, output)

%% Network weight changes
weight_change = training_output(1).weight_change;

figure
weight_change_plot(weight_change, N_train)

%% Helper functions
function spikes = convert_spiketimes(tspikes, dt, N, tmax, trial)
% CONVERT_SPIKETIMES converts the spiketimes from an array of spikecounts 
% and spike times to a boolean array of the spikes over time
% Input:
%   * tspikes = spike times
%   * dt = integration time constant (ms)
%   * N = number of neurons 
%   * tmax = length of target function
%   * trial = number of trials
% Output:
%   * spikes = struct with spike times converted to the spikemat

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
% CALC_MEAN_STATS calculates the mean firing rate and CV 
% over multiple trials for each neuron
% Input:
%   * stats = statistics of spike data
%   * N = number of neurons 
% Output:
%   * mean_fire_rate = mean firing rate for each neuron
%   * mean_Cv = mean coefficient of variation for each neuron

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

function [Z_out_seq, Zx_seq] = create_trial_seq(trials_num, sample_size, test_output)
% creates a sequence of test trials (output vs target signal)
% INPUT
% trials_num : total number of output trials
% sample_size : the sample size
% test_output : the test_output file
% OUTPUT
% Z_out : sequence of network output
% Zx_out : sequence of corresponding target signals

% select random trials of the sample size
sample_trials = randi([1 trials_num], 1,sample_size);

Z_out_seq = [];
Zx_seq = [];

for i = 1 : sample_size
    sample_trial = sample_trials(i);
    
    % append the network signals
    new_Z_out = test_output.Z_out{sample_trial}';
    Z_out_seq = [Z_out_seq new_Z_out];

    % append the target signals
    new_Zx = test_output.Zx{sample_trial};
    Zx_seq = [Zx_seq new_Zx];
end

end

%% Plot functions
function output_vs_target_plot(output, target)
% OUTPUT_VS_TARGET_PLOT plots the output signal vs the target signal

% plot the output vs target
hold on
plot(output, 'LineWidth', 1.5)
plot(target, 'LineWidth', 1.5)
ylim([-2.5 2.5])
xlim([-50 3500])
legend('network output', 'target')
title('Output vs target')
hold off
end

function trial_sequence_plot(Z_out_seq, Zx_seq)
% Plots a sequence of target and output signals
% INPUT
% Z_out_seq : sequence of output trials
% Zx_seq : sequence of corresponding target signals

% Plot the output vs target
hold on
plot(Z_out_seq)
plot(Zx_seq)
hold off
ylim([-2.5 2.5])
legend('network output', 'target')
xlabel('time in ms')
ylabel('network output')
end

function raster_plot(spikes, N, dt, samplesize)
% RASTER_PLOT creates a raster plot of a samplesize of random neurons 
% in the network. The spikes should be given as an boolean array 
% with the neuron spikes over time 
% Input:
%   * spikes = converted spikes
%   * N = number of neurons 
%   * dt = integration time constant (ms)
%   * samplesize = samplesize of random neurons in the network

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
% PSTH_PLOT makes a PSTH plot of the entire network activity. The spikes 
% should be given as an boolean array with the neuron spikes over time 
% Input:
%   * spikes = converted spikes
%   * dt = integration time constant (ms)

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
% FIRING_RATE_VS_CV_PLOT creates a scatter plot of the mean firing rate 
% vs CV of neurons
% Input:
%   * mean_fire_rate = mean firing rate for each neuron
%   * mean_Cv = mean coefficient of variation for each neuron
%   * feedback = feedback weights matrix

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
% WEIGHT_MATRICES_PLOT plots the static, feedback and 
% learnt weight matrices as well as the summed synaptic input of a neuron
% Input:
%   * static = static weight matrix
%   * feedback = feeback weight matrix
%   * ouput = network output
%   * N = number of neurons
%   * window = window size

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

function synaptic_input_plot(static, feedback, output)
learned = feedback * output';
final = static + learned;
sum_synaptic_input = sum(final, 2);

histogram(sum_synaptic_input, 'DisplayStyle', 'stairs', 'LineWidth', 2)
title('Distribution of synaptic input (without thalamus input)')
ylabel('# neurons')
xlabel('synaptic input')
end


function weight_change_plot(weight_change_1, N_train)
% plots the learnt output weight changes during the training trials
% Input:
%   * weight_change = N_trainx1 matrix of learnt output weight differences
%   * N_train = number of training trials

% plot the weight changes
hold on
plot(weight_change_1)
% plot(weight_change_2)
% plot(weight_change_3)
% plot(weight_change_4)

hold off
xlim([0 N_train])
xlabel('number of training trials')
ylabel('weight change')
% legend('epoch = 1', 'epoch = 2')
title('Changes of the learnt output weights during training')

end

