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
main_dir = 'Input_demo_var';
mkdir(['Output' f main_dir])

%% Save files
savename = 'filter';
savefolder_filter = ['.' f 'Output' f main_dir f savename];

savename = 'PSTH';
savefolder_PSTH = ['.' f 'Output' f main_dir f savename];

savename = 'thalamus';
savefolder_thalamus = ['.' f 'Output' f main_dir f savename];
%% Set the parameters for all three simulations
N = 2000;        % number of neurons
N_th = 200;     % number of thalamus neurons
N_train = 2;    % number of training trials
N_test = 2;     % number of validation trials
N_total = 1;    % number of epochs (only functions properly wiht one)

% varying parameters
Win = [0.5];       % scales the input weights
G = [10];         % scales the static weights
Q = [0];         % scales the feedback weights
Winp = [1];      % network sparsity
alpha = 0.05;    % learning rate

% logicals
FORCE = false;   % apply FORCE learning during trials
makespikes = false; % make the trial spiking structures 

% Dale's law
Pexc = 0; 

%% Filtered trace
Win = 100;
run_filter = run_sim_filter(N, N_th, N_train, N_test, N_total, Win, G,...
    Q, Winp, alpha, Pexc, FORCE, makespikes, savefolder_filter);

%% PSTH of trace
Win = 0.1;
run_trace = run_sim_PSTH(N, N_th, N_train, N_test, N_total, Win, G,...
    Q, Winp, alpha, Pexc, FORCE, makespikes, savefolder_PSTH);

%% Thalamus spikes
Win = 0.5;
run_thalamus = run_sim(N, N_th, N_train, N_test, N_total, Win, G, Q,...
    Winp, alpha, Pexc, FORCE, makespikes, savefolder_thalamus);

%% Plot the different inputs and raster plots
epoch = 1;
trial = 1;
dt = 0.05;

% Thalamus input
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

% Thalamus input simulation (no filter)
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

% Curve and angle trace simulations
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
%NETWORK_OUTPUT_PLOT 
%   Detailed explanation goes here

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
% input_type is 'neuron_input' or 'syn_input'

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
%RASTER_PLOT
%   Detailed explanation goes here

% Plot rasterplot of 200 neurons
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
%PSTH_PLOT 
%   Detailed explanation goes here

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


function angle_curve_plot(main_dir ,trial)
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


