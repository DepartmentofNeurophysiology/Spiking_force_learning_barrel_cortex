function [reservoir_input, target_function] = reservoir_input(SpikeTrainStruct, Ein, N, N_th, dat, rate, input_type)
% RESERVOIR_INPUT takes the nth thalamic spike times and makes neuron spiking
% input for the reservoir and returns this. 
% Input:
%   * SpikeTrainStruct = structure containing PSTH, SpikeTimes 
% and SpikeCount of a number of trials
%   * Ein = the input weights matrix
%   * N = number of trials
%   * N_th = selects the trials in the struct for which you want 
% to make the reservoir input and the target function
%   * dat = contains the pole location
%   * rate = contains the rate of the intermediate poisson firing
% Output: 
%   * reservoir_input = structure containing the neuron spiking input for the reservoir 
%   * target_funciton = structure containing the target function
% Helper functions:
%   * scale_input

%% Target pulse
trial_len = length(SpikeTrainStruct{1}.PSTH{1}); % length of trial
pulse_length = 1000;    % length of pulse (ms)
amp = 2;                % amplitude of pulse (V)
decay = 200;            % decay of exponential pulse (ms)
constant = 100;         % how long the pulse is kept constant (ms)
reset = 1500;           % time inbetween the trials (ms)
dt = 0.001;             % integration time constant (ms)                    
T_vec = 0:dt:reset*dt;  % vector with each time step	

% if SpikeTrainStruct{1, 1}.first_touch == 0
%     first_touch = 1000;
% elseif SpikeTrainStruct{1, 1}.first_touch < 0
%     first_touch = 0;
% else
%     first_touch = SpikeTrainStruct{1, 1}.first_touch;
% end

first_touch = 500;  % start of first touch (ms)

z_all = [];
start_early = 500;  % how long to start before the end of the input (ms)

%% Make the reservoir input and the target function
l_t = trial_len + reset;
z_t = zeros(l_t,1);

% from 500ms (first toucht|) to the end of the input signal +
% a constant of 100ms
z_t(first_touch:trial_len+constant) = -dat(1)*amp*ones(length(first_touch:trial_len+constant),1);
% from the end of the input signal + a constant of 100ms until the end of
% the pulse length + constant of 100ms
z_t(trial_len+constant:1:trial_len+pulse_length+constant) =  -dat(1)*amp*exp(-(0:1:pulse_length)./decay);

z_all = [z_all z_t'];
target_function = z_all;
%% Neuron input
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
   
    reservoir_input = [scaled_neuron_input ones(N, reset)*mu];
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
   
    reservoir_input = [scaled_neuron_input ones(N, reset)*mu];

% spike condition reservoir input
elseif strcmp('spikes', input_type)
    spikes_mat = zeros(N_th, trial_len);
    for n = 1:N_th
        spiketimes = SpikeTrainStruct{1}.SpikeTimes{n};
        spikes_mat(n, spiketimes) = ones(length(spiketimes) ,1);
    end
    neuron_input = Ein*spikes_mat;
    
    % make the poisson input
    % makes the poisson spikes random
    rng('shuffle') 
    
    % makes the poisson spikes fixed
    %rng(0) 
    
    for n = 1:N_th
        vt = rand(size(T_vec) - [0 1]);
        spikes = (rate*dt) > vt;
        thalamus_poisson(n).spike_times = find(spikes == 1);
    end
    poisson_input = make_poisson_spikes_weighted(N, Ein, reset, thalamus_poisson);
    reservoir_input = [neuron_input poisson_input];
end

end

%% Helper functions
function scaled = scale_input(neuron_input, scale_values, name)
% SCALE_INPUT scales the input of a given signal to the spike signal
    temp = getfield(scale_values, name);
    spikes = scale_values.spikes;
    scaled = (neuron_input - temp.mu)./temp.sigma;
    scaled = (scaled*spikes.sigma) + spikes.mu;
end


