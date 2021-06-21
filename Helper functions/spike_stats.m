function [ A_t, ISI, Cv ] = spike_stats( tspike, trial_length , N )
% SPIKE_STATS calculates the average activity, interspike intervals
% and the coefficient of variation of all neurons in the reservoir
% Input:
%   * tspikes = spike times (s)
%   * trial_length = length of trial (s)
%   * N = number of neurons in the reservoir
% Output:
%   * A_t = average activity / average firing rate 
%   * ISI = interspike intervals
%   * Cv = coefficient of variantion

% remove the zero values
not_null = tspike(:, 1) ~= 0;
tspike = tspike(not_null, :);

% calculate the spikecounts per neuron
neurons = tspike(:, 1);
[A_t, ~] = histcounts(neurons, 1:1:N + 1);
A_t = A_t / (trial_length*10^-3);

% ISI
ISI = [];
for n=1:N
    index = find(tspike(:,1)==n);
    spike_neuron(n).time = tspike(index,2);
    spike_neuron(n).ISI = diff(spike_neuron(n).time);
    ISI = [ISI diff(spike_neuron(n).time)'];
end

% Coefficient of variation
Cv = zeros(N,1);
for n = 1:N
    Cv(n,1) = std(spike_neuron(n).ISI)/mean(spike_neuron(n).ISI);
end


end

