function [ neuron_input] = make_poisson_spikes_weighted(N, Ein, T, thalamus_poisson)
% MAKE_POISSON_SPIKES_WEIGHTED makes weighted poisson spikes by taking 200
% spikes as input over a time T and multiplies it with the N x 200 weights
% matrix Ein to output a N x T array of 'weighted' spikes which can be 
% used as input for the neural network.
% Input:
%   * N = number of trials
%   * Ein = the input weights matrix
%   * T = time inbetween the trials (ms)
%   * thalamus_poisson = thalamus poisson spikes
% Output:
%   * neuron_input = N x T array of weighted poisson spikes

spike_array{1}.trial = zeros(200, T);

for i = 1:200
    spike_array{1}.trial(i,thalamus_poisson(i).spike_times) = ones(length(thalamus_poisson(i).spike_times) ,1);
end

spike_array{1,1}.neuron_input = zeros(N,length(spike_array{1,1}.trial));
for t=1:length(spike_array{1,1}.trial)
    index_input = find(spike_array{1,1}.trial(:,t) == 1);  % Find thalamic input neurons that have spiked
    if length(index_input) ~= 0
        spike_array{1,1}.neuron_input(:,t) = sum(Ein(:,index_input),2);
    end
end

neuron_input = spike_array{1}.neuron_input;

end
