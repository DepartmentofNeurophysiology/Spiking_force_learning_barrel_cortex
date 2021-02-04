# Spiking FORCE learning barrel cortex

### Department of Neurophysiology, Donders Institute, Radboud University

## Description
This is MATLAB code which FORCE trains a spiking neural network with different whisker traces as input.

Most important scripts and functions of the code are:

1. **'run_sim.m'**
 This function prepares and runs the simulation. This function sets the fixed parameters, prepares the trials, and performs an optional parameter sweep.

2. **'LIF_training'** (*'Network functions' folder*)
This function is called in *'run_sim.m'* and prepares the weights and input of the spiking network, and runs the test and train trials of the network. During test and train trials the *'LIF_spiking_network'* function is called to simulate the network with the trial input and network weights. The network output is currently being saved in this function as well.

3. **'LIF_spiking_network'** (*'Network functions' folder*)
 This the function where the dynamics of a spiking neural network are computed and FORCE training is applied.

The thalamic spikes, ConvTraces and PSTH signal were made according to the 'cortical-representation-of-touch-in-silico' model and are saved per trial as a structure in the **'Spiking structures'** folder. The name of each structure is the trial number.

Additional functions and data can be found in the **'Helper data'**, **'Helper functions'**, **'Network functions'** and **'Thalamus functions'** folders.
A quick overview of different plots used to study the networks can be found in the **'Plot functions'** folder.

## Run code
To run an example of the network, run *'quick_example.m'*.

Important parameters set in this script are:
- The number of neurons in the reservoir and the number of thalamus neurons
- The number of train and test trials (always a 1:1 ratio of random proximal and distal trials)
- The number of epochs
- The scaling parameters for the weights (Win, G, Q). To initiate a parameter sweep simply enter multiple values for the parameters.
- The learning rate alpha
- A boolean whether FORCE learning is applied or not
- A boolean whsether the input spikes are made or load from the 'Spiking Structures' folder
- The input type, choose between 'ConvTrace', 'PSTH' or 'spikes'
- The percentage of excitatory neurons. When not 0, Dale's law is automatically applied.
