function [reservoir_input, target_function ] = noisy_reservoir_input(SpikeTrainStruct, Ein, N, N_th, dat, rate, input_type)
% Function that takes the nth thalamic spike times and makes neuron spiking
% input for the reservoir and returns this. WITH NOISY BACKGROUND

% SpikeTrainStruct: structure containing PSTH, SpikeTimes and SpikeCount of
% a number of trials.
% n: selects the trials in the struct for which you want to make the
% reservoir input and the target function.
% Ein: is the input weights matrix.
% dat: contains the pole location.
% rate: contains the rate of the intermediate poisson firing. 

% output: is a structure containing the spiking reservoir input and the
% target function.

%% Target pulse
trial_len = length(SpikeTrainStruct{1}.PSTH{1});
pulse_length = 1000; % Length of pulse
amp = 2; % Amplitude of pulse
decay = 200; % decay of exponential pulse
constant = 100; % How long the pulse is kept constant
reset = 1500; % Time inbetween the trials
dt = 0.001;                                     % 1 msec                          
T_vec = 0:dt:reset*dt;                         % a vector with each time step	

% if SpikeTrainStruct{1, 1}.first_touch == 0
%     first_touch = 1000;
% elseif SpikeTrainStruct{1, 1}.first_touch < 0
%     first_touch = 0;
% else
%     first_touch = SpikeTrainStruct{1, 1}.first_touch;
% end
first_touch = 500;

z_all = [];
start_early = 500; % How long to start before the end of the input

%% Make the reservoir input and the target function
l_t = trial_len + reset;

z_t = zeros(l_t,1);

% from 500ms (first toucht|) to the end of the input signal + a constant of
% 100ms
z_t(first_touch:trial_len+constant) = -dat(1)*amp*ones(length(first_touch:trial_len+constant),1);
% from the end of the input signal + a constant of 100ms until the end of
% the pulse length + constant of 100ms
z_t(trial_len+constant:1:trial_len+pulse_length+constant) =  -dat(1)*amp*exp(-(0:1:pulse_length)./decay);

z_all = [z_all z_t'];
target_function = z_all;
%% Background noise 
 % make the poisson input
 rng('shuffle')
 trial_time = trial_len + reset;
 trial_time_vec = 0:dt:trial_time*dt;
 for n = 1:N_th
     vt = rand(size(trial_time_vec) - [0 1]);
     spikes = (rate*dt) > vt;
     thalamus_poisson(n).spike_times = find( spikes == 1);
 end
 poisson_noise = make_poisson_spikes_weighted(N, Ein, trial_time, thalamus_poisson);

%% Neuron Input
% load the scale values file if necessary
%f = filesep;
%addpath(['..' f 'Helper data'])
if ~strcmp('spikes', input_type)
    filename = 'scale_values.mat';
    if ~exist(filename, 'file')
        error([filename ' is not in the helper data folder'])
    else
        file = load(filename);
        scale_values = file.ScaleVal;
   end 
end

% ConvTrace condition reservoir input
if strcmp('ConvTrace', input_type)
    ConvTrace = zeros(N_th, trial_len);
    for n = 1:N_th
        ConvTrace(n, :) = SpikeTrainStruct{1}.ConvTrace{n};
    end
    neuron_input = Ein*ConvTrace;
    
    % scale input signal
    scaled_neuron_input = scale_input(neuron_input, scale_values, 'ConvTrace');
    mu = scale_values.spikes.mu;
    
    %[scaled_neuron_input, mu] = scale_input_alt(SpikeTrainStruct, neuron_input, Ein, N, N_th);
    
   
    reservoir_input = [scaled_neuron_input ones(N, reset)*mu];
    reservoir_input = reservoir_input + poisson_noise;
    %disp('ConvTrace')

% PSTH condition reservoir input
elseif strcmp('PSTH', input_type)
    PSTH = zeros(N_th, trial_len);
    for n = 1:N_th
        PSTH(n, :) = SpikeTrainStruct{1}.PSTH{n};
    end
    neuron_input = Ein*PSTH;
    scaled_neuron_input = scale_input(neuron_input, scale_values, 'PSTH');
    mu = scale_values.spikes.mu;
   
    %[scaled_neuron_input, mu] = scale_input_alt(SpikeTrainStruct, neuron_input, Ein, N, N_th);
    
    reservoir_input = [scaled_neuron_input ones(N, reset)*mu];
    reservoir_input = reservoir_input + poisson_noise;
    %disp('PSTH')

% spike condition reservoir input
elseif strcmp('spikes', input_type)
    spikes_mat = zeros(N_th, trial_len);
    for n = 1:N_th
        spiketimes = SpikeTrainStruct{1}.SpikeTimes{n};
        spikes_mat(n, spiketimes) = ones(length(spiketimes) ,1);
    end
    neuron_input = Ein*spikes_mat;
    reservoir_input = poisson_noise;
    reservoir_input(:, 1:trial_len) =  reservoir_input(:, 1:trial_len) + neuron_input;
    %disp('Spikes')
end

end

%% Helper functions
function scaled = scale_input(neuron_input, scale_values, name)
% scales the input of a given signal to the spike signal

temp = getfield(scale_values, name);
spikes = scale_values.spikes;
scaled = (neuron_input - temp.mu)./temp.sigma;
scaled = (scaled*spikes.sigma) + spikes.mu;
end

function [scaled, mu] = scale_input_alt(SpikeTrainStruct, neuron_input, Ein, N, N_th)
trial_len = length(SpikeTrainStruct{1}.PSTH{1});
spikes_mat = zeros(N_th, trial_len);
for n = 1:N_th
    spiketimes = SpikeTrainStruct{1}.SpikeTimes{n};
    spikes_mat(n, spiketimes) = ones(length(spiketimes) ,1);
end
spikes_mat = Ein*spikes_mat;
filt_spikes_mat = double_exponential_filter(spikes_mat,N);
mu = mean(filt_spikes_mat, 'all');
sigma = std(filt_spikes_mat, 0, 'all');

scaled = (neuron_input - mean(neuron_input, 'all'))./std(neuron_input, 0, 'all');
scaled = (scaled*sigma) + mu;
%{
figure 
subplot(2,1,1)
plot(mean(filt_spikes_mat))
title('spikes')
subplot(2,1,2)
plot(mean(scaled))
title('scaled')
pause;
%}
end

function f2_t = double_exponential_filter(neuron_input, N)
%DOUBLE_EXPONENTIAL_FILTER Summary of this function goes here
%   Detailed explanation goes here
 % run the spikes through a double exponential filter
 dt = 0.05;
 tau_r = 2;
 tau_d = 50;
 %T = length(SpikingStruct{1}.PSTH{1});
 T = size(neuron_input,2);
 nt = T/dt;
 
 %{
 spikearray = zeros(N_th, T);   
 for j = 1:N_th
     spikearray(j ,SpikingStruct{1}.SpikeTimes{j}) = ones(length(SpikingStruct{1}.SpikeTimes{j}) ,1);
 end
 spikearray = Ein*spikearray;
 %}   
    
 input = zeros(N, nt);
 for t = 1:nt 
     if mod(t, 1/dt) == 0
         input(:, t) = neuron_input(:, t*dt); 
     end
 end
 
 f = zeros(N,1);
 f2 = zeros(N,1);
 f2_t = zeros(size(input));
    
 for j = 1:nt
     f2_t(:, j) = f2;
     f2 = f2 * exp(-dt/tau_r) + f*dt;
     f = f * exp(-dt/tau_d) + input(:, j)/(tau_r*tau_d);
 end
end
