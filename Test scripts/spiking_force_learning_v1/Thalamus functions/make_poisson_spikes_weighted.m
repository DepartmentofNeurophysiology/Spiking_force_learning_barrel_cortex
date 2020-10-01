function [ neuron_input] = make_poisson_spikes_weighted(N, Ein, T, thalamus_poisson)
% Takes 200 spikes as input over a time T and multiplies it with the N x
% 200 weights matrix Ein to ouput a N x T array of 'weighted' spikes which
% can be used as input for the neural network.

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
