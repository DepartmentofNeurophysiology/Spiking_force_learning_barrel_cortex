function voltage_trace(vtrace, T)
%VOLTAGE_TRACE Summary of this function goes here
%   Detailed explanation goes here
% Neuron potential trace
x = 1: T; 
plot(x, vtrace)
xlim([0 T])
legend
xlabel('t in ms')
ylabel('membrane potential in mV')
end

