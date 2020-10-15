function [ error, BPhi, zx, z_out, tspike ] = spiking_network( param, weights, input_struct, FORCE, iter )
% This is the main spiking neural network function. As input it receives
% some essential network paramters, connections weights for the reservoir,
% output, and feedback. It also receives the thalamic input and the target
% function. FORCE is a parameter that can be turned on or off, which
% enables the network to be trained or not. 

% Parameters that are returned are the error, the output weights (BPhi), the
% network output (zx), the target output (z_out), and the spike times
% (tspike).

N = param.N; 
alpha = param.alpha; % learning rate
BIAS = param.BIAS;
step = param.step; 
dt = param.dt; % integration step time
td = param.td; % synaptic decay

%% Network parameters
tref = 2; %Refractory time constant in milliseconds 
tm = 10; %Membrane time constant 
vreset = -65; %Voltage reset 
vpeak = -40; %Voltage threshold 
rng(1);
tr = 2; % synaptic rise time
%% Storage parameters and some others
Pinv = eye(N)*alpha; %initialize the correlation weight matrix for RLMS
IPSC = zeros(N,1); %post synaptic current storage variable 
h = zeros(N,1); %Storage variable for filtered firing rates
r = zeros(N,1); %second storage variable for filtered rates 
hr = zeros(N,1); %Third variable for filtered rates 
JD = 0*IPSC; %storage variable required for each spike time  
v = vreset + rand(N,1)*(30-vreset); %Initialize neuronal voltage with random distribtuions
tlast = zeros(N,1); %This vector is used to set the refractory times 
z = 0; % Initial z value
%% Define weights
BPhi = weights.output_weights; % learnt output weights
OMEGA = weights.static_weights; % static weights
E = weights.feedback_weights; % feedback weights
%% Define input and time
zx = input_struct.target_function; % target function
T = length(zx);
nt = T/dt;
% Define some storage arrays
input = zeros(N,nt);
z_out = zeros(T,1);
% Adjust the time resolution of input to integration step of the
% network.
for t = 1:nt
    if mod(t,1/dt) == 0
        input(:,t) = input_struct.reservoir_input(:,t*dt);
    end
end
in = 1;
tspike = zeros(1,2); %Storage variable for spike times
ns = 0;
%% MAIN NETWORK LOOP
for i = 1:1:nt % for every integration step time T = length(zx), nt = T/dt
    if mod(i,1/dt) == 0 % this loop makes sure that only every 1 ms a new data point is presented
        in=in+1;
    end
    if in > T
        break
    end
    I = IPSC + E*z + BIAS ; %Neuronal Current
    dv = (dt*i>tlast + tref).*(-v+I)/tm; %Voltage equation with refractory period
    v = v + dt*(dv);
    index = find(v>=vpeak);  %Find the neurons that have spiked
    %Store spike times, and get the weight matrix column sum of spikers
    if length(index)>0
        JD = sum(OMEGA(:,index),2); %compute the increase in current due to spiking
        if FORCE == 0 
            %tspike(ns+1:ns+length(index),:) = [index,0*index+dt*i]; % Store spike times
            spikes = [index,0*index+dt*i];
            tspike = [tspike; spikes ];
            %ns = ns + length(index);  % total number of spikes so far
        end
        
    end
        
    %Implement RLMS with the FORCE method
    z = BPhi'*r; %approximant
    err = z - zx(:,in); %error
    %Store the output of the network
    z_out(in,1) = z;
    % RLS
    if FORCE
        if mod(i,step) == 1 % every step iterations (in this case 20) the output weights are updated
            if zx(in) ~= 0 % Makes sure there is only RLS when the target is non-zero
                cd = Pinv*r;
                BPhi = BPhi - (cd*err');
                Pinv = Pinv -((cd)*(cd'))/( 1 + (r')*(cd));
            end
        end
    end
    tlast = tlast + (dt*i -tlast).*(v>=vpeak);  %Used to set the refractory period of LIF neurons
    
    % filter the spikes that go through the network
    IPSC = IPSC*exp(-dt/tr) + h*dt;
    h = h*exp(-dt/td) + JD*(length(index)>0)/(tr*td) + input(:,i)/(tr*td);  % THE LAST TERM ARE THE THALAMIC SPIKES
    
    % filter the spikes of the synaptic output
    r = r*exp(-dt/tr) + hr*dt; 
    hr = hr*exp(-dt/td) + (v>=vpeak)/(tr*td);
        
    v = v + (30 - v).*(v>=vpeak); % Sets the voltage to spike 
    v = v + (vreset - v).*(v>=vpeak); %reset with spike time interpolant implemented.
end
%%
% Mean Square Error between network output and target
s_t = length(zx) - 800 - 500; % only calculate error when output is trained
error = immse(z_out(s_t:end)', zx(s_t:end));

end

