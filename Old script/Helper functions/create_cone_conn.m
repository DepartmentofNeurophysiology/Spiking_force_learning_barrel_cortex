function [ Ein] = create_cone_conn( N_th, N, sigmaffwd, K_out, Win )
% N_th: This is the number of thalamus neurons
% N: Number of neurons in the reservoir.
% sigmaffd: Radius of the gaussian projection.
% K_out: Number of connections each thalamus neuron makes.
% Win: Scale the input weights.

[ indices ] = disk_connections( N_th, N, sigmaffwd, K_out); % Find the indices connecting the the input to the reservoir

%  Make the input weights matrix
                                          
Ein = zeros(N,N_th);                        
for i = 1:N_th
    Ein(indices(i,:), i) = Win * rand(K_out,1);     % Give all neurons to connect with a random weight.
end
end

