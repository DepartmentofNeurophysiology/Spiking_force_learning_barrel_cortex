function [error, output_weights, Zx, Z_out, tspikes] = LIF_spiking_network_no_filt(param, weights, thalamus_input, target, FORCE)
% LIF_SPIKING_NETWORK_NO_FILT computes the dynamics of the spiking 
% neural network without synaps filter, and applies FORCE learning, 
% only for input types 'ConvTrace' or ' PSTH'
% Input:
%   * param = network parameters
%   * weights = weights parameters
%   * thalamus_input = thalamic input
%   * target = target function
%   * FORCE = 0 or 1, apply FORCE learning or not
% Output:
%   * error = mean square error between network output and target
%   * output_weights = network output weights
%   * Zx = target function
%   * Z_out = network output
%   * tspikes = spike times

disp('no synapse filter')
%% Network parameters
% input parameters
N = param.N;
alpha = param.alpha;
Ibias = param.Ibias;
step = param.step;
dt = param.dt;
tau_d = param.tau_d;
tau_r = param.tau_r; 

% static parameters
tref = 2;       % refractory time period (ms)
tau_m = 10;     % membrane time constant (ms)
vreset = -65;   % reset potential (V)
vthresh = -40;  % threshold potential (V) 
rng(1);         % every time the same random distribution 

%% Define weights
output_weights = weights.output;
static_weights = weights.static;
feedback_weights = weights.feedback;

%% Target function and time
Zx = target;

% set time and timesteps
T = length(target);
nt = T/dt;

%% Input
input = zeros(N, nt);

% adjust the time resolution of input to dt of the network
for t = 1:nt 
    
    if mod(t, 1/dt) == 0
        input(:, t) = thalamus_input(:, t*dt); 
    end
end

%% Storage parameters
% post synaptic current and sum synaptic input
Ipsc = zeros(N,1); 
Ispikes = 0*Ipsc; 

% initialize neuronal voltage with random distribtuions
v = vreset + rand(N,1)*(30-vreset); 

% refactory times, spiketimes and total number of spikes
tlast = zeros(N,1); 
tspikes = zeros(4*nt, 2);
%tspikes = [];
nspikes = 0;

% first, second and third storage variables for filtered rates
h = zeros(N,1); 
r = zeros(N,1); 
hr = zeros(N,1); 
 
% initialize z output and storage value
Z = 0; 
Z_out = zeros(T,1);

% initialize the correlation weight matrix for RLMS 
Pinv = eye(N)*alpha; 

%% MAIN NETWORK LOOP
for i = 1:1:nt
    
    % present a new datapoint every 1 ms 
    in = ceil(i * dt);
   
    % update the input current of the neurons
    I = Ipsc + feedback_weights*Z + Ibias + input(:,i);
    dv = (dt*i > tlast + tref).*(-v + I)/tau_m;
    v = v + dt*(dv);
    
    % find neurons that spiked
    spike_index = find(v >= vthresh);
    
    % get the increase in current due to spiking and store the spike times
    if ~isempty(spike_index)
       
        Ispikes = sum(static_weights(:, spike_index), 2);
        
        if ~FORCE
            tspikes(nspikes+1:nspikes+length(spike_index),:) =...
                [spike_index, 0*spike_index+dt*i];
            
            nspikes = nspikes + length(spike_index);
            
            % alternative spike storage
            %spikes = [spike_index,0*spike_index+dt*i];
            %tspikes = [tspikes; spikes];
        end
    end
    %% Implement RLMS with the FORCE method
    % calculate the network output and error
    Z = output_weights'*r;
    err = Z - Zx(:, in); 
    
    % store the output
    Z_out(in, 1) = Z;
    
    % RLMS, check if FORCE learning applies and if target input is not 0
    if FORCE && Zx(in) ~= 0
        
        % every step iterations the output is updated
        if mod(i, step) == 1
            cd = Pinv*r;
            output_weights = output_weights - (cd*err');
            Pinv = Pinv -((cd)*(cd'))/( 1 + (r')*(cd));
        end
    end
     
    % set the refractory period of the neurons
    tlast = tlast + (dt*i - tlast).*(v >= vthresh);
    
    % apply the double exponential filter for the postsynaptic current
    Ipsc = Ipsc*exp(-dt/tau_r) + h*dt;
    h = h*exp(-dt/tau_d) + Ispikes*(~isempty(spike_index))/(tau_r*tau_d);  
    
    % filter the spikes of the synaptic output
    r = r*exp(-dt/tau_r) + hr*dt; 
    hr = hr*exp(-dt/tau_d) + (v>=vthresh)/(tau_r*tau_d);
    
    % spike and reset the voltage of the neurons that fired
    v = v + (30 - v).*(v >= vthresh);
    v = v + (vreset - v).*(v >= vthresh);
end

% Mean Square Error between network output and target
%s_t = length(Zx) - 800 - 500; 
%error = immse(Z_out(s_t:end)', Zx(s_t:end));
error = immse(Z_out(1000:end)', Zx(1000:end));

% remove the zeros from the tspikes struct
tspikes = tspikes((tspikes(:, 1) ~= 0), :);

end

