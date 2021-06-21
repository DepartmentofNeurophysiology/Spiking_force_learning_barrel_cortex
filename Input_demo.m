%% Input demo
% Runs short simulations of the networks with different inputs.
% Input conditions:
%   * Filtered whisker traces (Convtrace)
%   * Whisker trace PSTH
%   * Thalamus spikes
% Helper functions
%   * convert spiketimes
% Plot functions
%   * mean_input_plot
%   * raster_plot
%   * PSTH_plot
%   * angle_curve_plot

disp('Remember to adjust the training files to save the network input')
pause;
%% Add paths
f = filesep;
main_dir = 'Input_demo_scaled';
mkdir(['Output' f main_dir])

%% Save files
savename = 'filter';
savefolder_filter = ['.' f 'Output' f main_dir f savename];

savename = 'PSTH';
savefolder_PSTH = ['.' f 'Output' f main_dir f savename];

savename = 'thalamus';
savefolder_thalamus = ['.' f 'Output' f main_dir f savename];
%% Set the parameters for all three simulations
% input parameters
N = 2000;        % number of neurons in the reservoir 
N_th = 200;      % number of thalamus neurons
N_train = 600;   % number of training trials
N_test = 100;    % number of validation trials
N_total = 1;     % number of epochs 

% weight scaling parameters
% Win = 0.5;     % the input weights
G = 10;          % the static weights
Q = 1;           % the feedback weights
Winp = 1;        % network sparsity
alpha = 0.05;    % learning rate

% logicals
FORCE = true;       % if TRUE; apply FORCE learning during trials
makespikes = true;  % if TRUE; make the trial spiking structures 

% Dale's law
Pexc = 0;   % percentage of excitatory neurons; set to 0 to ignore Dale's law restrictions

% Input 
input_type = 'spikes';  % options: 'ConvTrace', 'PSTH', 'spikes'
%% Filtered trace
Win = 0.2; % the input weight
run_filter = run_sim(N, N_th, N_train, N_test, N_total, Win, G,...
    Q, Winp, alpha, Pexc, FORCE, makespikes, input_type, savefolder_filter);

%% PSTH of trace
Win = 0.1; % the input weight
run_trace = run_sim(N, N_th, N_train, N_test, N_total, Win, G,...
    Q, Winp, alpha, Pexc, FORCE, makespikes, input_type, savefolder_PSTH);

%% Thalamus spikes
Win = 0.05; % the input weight
run_thalamus = run_sim(N, N_th, N_train, N_test, N_total, Win, G, Q,...
    Winp, alpha, Pexc, FORCE, makespikes, input_type, savefolder_thalamus);

%% Plot the different inputs and raster plots
epoch = 1;  % number of epochs
trial = 1;  % number of trials
dt = 0.05;  % integration time constant (ms)

% thalamus input
figure
angle_curve_plot(main_dir, trial);

subplot(7, 1, 3)
mean_input_plot(savefolder_thalamus, trial, 'neuron_input');

subplot(6, 1, 4)
mean_input_plot(savefolder_thalamus, trial, 'syn_input');

subplot(7, 1, [5 6])
spikes = convert_spiketimes(savefolder_thalamus, N, dt, trial, epoch);
raster_plot(spikes, N, dt);

subplot(7, 1, 7)
PSTH_plot(spikes, dt);

sgtitle('Thalamus input')

% thalamus input simulation (no filter)
figure
angle_curve_plot(main_dir, trial);

subplot(7, 1, 3)
mean_input_plot(savefolder_filter, trial, 'neuron_input');

subplot(7, 1, [4 5 6])
spikes = convert_spiketimes(savefolder_filter, N, dt, trial, epoch);
raster_plot(spikes, N, dt);

subplot(7, 1, 7)
PSTH_plot(spikes, dt);

sgtitle('Convtrace input')

% curve and angle trace simulations
figure
angle_curve_plot(main_dir, trial);

subplot(7, 1, 3)
mean_input_plot(savefolder_PSTH, trial, 'neuron_input');

subplot(7, 1, [4 5 6])
spikes = convert_spiketimes(savefolder_PSTH, N, dt, trial, epoch);
raster_plot(spikes, N, dt);

subplot(7, 1, 7)
PSTH_plot(spikes, dt);

sgtitle('PSTH input')

%% Helper functions
function spikes = convert_spiketimes(input_folder, N, dt, trial, epoch)
% CONVERT_SPIKETIMES converts the spiketimes from an array of spikecounts 
% and spike times to a boolean array of the spikes over time 
% Input:
% * input_folder: struct with spike times
% * N: number of neurons
% * dt: integration time constant (ms)
% * trial: number of trials
% * epoch: number of epochs
% Output:
% * spikes: struct with spiketimes converted to the spikemat

%% Load everything from the input folder
f = filesep;
structs = dir(fullfile(input_folder, '*.mat'));

for i = 1 : length(structs)
    struct = structs(i).name;
    
    load([input_folder f struct]);
end

test_output = training_output(epoch).test_output;
tmax = length(test_output.Zx{trial});

% get the spike times
tspikes = test_output.tspikes{trial};

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

%% Plot Functions
function mean_input_plot(input_folder, trial, input_type)
% MEAN_INPUT_PLOT plots the mean of the input
% Input:
%   * input_folder = struct with input data
%   * trial = number of trials
%   * input_type = 'neuron_input' or 'syn_input'

f = filesep;
structs = dir(fullfile(input_folder, '*.mat'));

for i = 1 : length(structs)
    struct = structs(i).name;
    
    load([input_folder f struct]);
end

input_mat = getfield(input_save{trial}, input_type);

mean_input = mean(input_mat,1);

if strcmp(input_type, 'syn_input')
    plottitle = 'Mean network input after synapse filter';
else
    plottitle = 'Mean network input';
end

plot(mean_input)
title(plottitle)
xlim([-100 4000])

end

function raster_plot(spikes, N, dt)
% RASTER_PLOT creates a raster plot of a samplesize of random neurons 
% in the network. The spikes should be given as an boolean array 
% with the neuron spikes over time 
% Input:
%   * spikes = converted spikes
%   * N = number of neurons 
%   * dt = integration time constant (ms)

% plot rasterplot of 200 neurons
hold on
title('Network raster plot and PSTH')
mid = ceil(N/2);
window = 100;
for nn=(mid-window):(mid+window)
    spikeindices = find(squeeze(spikes(1,nn,:))==1);
    plot(dt*spikeindices, nn*ones(1,length(spikeindices)), '.k')
end
xlim([-100 4000])
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
plot(PSTH_squeeze);
ylabel('PSTH')
xlabel('time (ms)') 
ylim([0 300])
xlim([-100 4000])

end


function angle_curve_plot(main_dir, trial)
% ANGLE_CURVE_PLOT plots the angle and curvature of the wisking data
% Input: 
%   * main_dir = main directory of input folder
%   * trial = number of trials

f = filesep;
% Load everything from the input folder
filename = ['Output' f main_dir f 'angle_whisk.mat'];

load(filename);

curve = input_save{trial}.curve;
angle = input_save{trial}.angle;

subplot(7, 1, 1)
plot(curve, 'LineWidth', 1.5)
title('Curvature')
ylim([-0.01 0.01])
xlim([-100 4000])
xlabel('time in ms')
ylabel('Kappa')

subplot(7, 1, 2)
plot(angle, 'LineWidth', 1.5)
title('Angle')
%ylim([-0.5 0.5])
xlim([-100 4000])
%xlabel('time in ms')
ylabel('Theta')

end


