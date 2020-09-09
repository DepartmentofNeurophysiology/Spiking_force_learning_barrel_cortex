%% Test and visualize the a single spiking neuron
% Plot the membrane potential of a neuron that receives an input current
%% Set the parameters
vrest = -70; % rest potential
tm = 10; % membrane time constant
T = 200; % total time in ms
tstart = 50; % start of the current pulse
tend = 150; % end of the current pulse
Ia = 15; % amplitude of the current pulse in nA
Ibias = - 70; % current bias
dt = 1; % time steps in ms

Vtrace = zeros(1,T);

%%
v = vrest;

for i = 1:dt:T 
    % calculate the current
    I = ((i >= tstart & i <= tend) * Ia) + Ibias; 
    
    % calculate the change in voltage
    dv = (-v + I)/tm; 
    v = v + dt*dv;
    
    % store the voltage and current change
    Vtrace(i) = v;
    Itrace(i) = I - Ibias;
end

%% Plot of the neuron
figure
x = linspace(0,T,T);
y1 = Vtrace;
y2 = Itrace;

% voltage trace
subplot(2,1,1)
plot(x, y1)
axis([0 200 -75 -50])
ylabel('membrane potential in mV')
xlabel('t in ms')
title('simulation of a single spiking neuron')

% current trace
subplot(2,1,2)
plot(x, y2)
axis([0 200 0 10])
ylabel('current pulse in nA')
xlabel('t in ms')