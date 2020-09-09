# Spiking FORCE learning barrel cortex

### Department of Neurophysiology, Donders Institute, Radboud University

## Description
This is MATLAB code which FORCE trains a spiking neural network with thalamic spikes as input.

Most important scripts and functions of the code are:

1. **'main_script_train_some'**
 This the main script from which (some) parameters can be defined and networks can be made and trained.
 This is also where the data of trained and tested networks can be saved.
 At the end of the script, the network output and network spikes can also be plotted.

2. **'LIF_train_some'** (*'Network loop functions' folder*)
 This is an important function which is called in the *'main_script_train_some'* script.
 Here a network with a chosen parameter combination is trained for a number of trials, after which it is tested.
 This is also where the important *'spiking_network'* function is called.

3. **'spiking_network'** (*'Network loop functions' folder*)
 This the function where the dynamics of a spiking neural network are computed. It is called in the *'LIF_train_some'* function.
 The dynamics of the network are computed for one trial each time.
 The parameter 'FORCE' is passed to this function, and tells the function to FORCE train the output weights or not.

The thalamic spikes were made according to the 'cortical-representation-of-touch-in-silico' model and are saved per trial as a structure in the **'Spiking structures'** folder. The name of each structure is the trial number.

Additional functions and data can be found in the **'Helper data'**, **'Helper functions'**, **'Network loop functions'** and **'Thalamus functions'** folders.

## Run code
Here follows a brief description of how the script can be run:

The *'main_script_train_some'* is the main script from which a network can be made and trained. Here, most important
parameters can be defined. Most importantly, you can set the parameter combinations with which you want to train networks.
These parameters can be defined as a 1D array, and need to be passed to the *'all_comb'* function. An example is
if you want to train networks with all combinations of G = [8, 9] and Q = [1, 2] then 'all_comb(G,Q)' will return all possible combinations of G and Q with which networks will be trained.

All parameters are passed as one single structure to the *'LIF_train_some'* function. This is a function that trains
a network and then tests it. All important data of the trained network is returned in the structure 'run'.
In this structure, all important data per network is stored, where for example, run(1).parameters_out will have all data of a network trained with G = 8 and Q = 1.

In the *'LIF_train_some'* function, the static weights and the feedback weights are initialized.
It is important that at the beginning of this function all parameters are defined correctly.
THIS MAY VARY DEPENDING ON THE PARAMETER SEARCH YOU ARE PERFORMING.

Then the network is trained. For each trial the network is trained, the spiking input is loaded from the 'Spiking structures' folder. This is done using the *'reservoir_input'* function for each trial.

The dynamics and the training of the network happens in the *'spiking_network'* function. Importantly, the FORCE parameter can be set to true/false. This either trains the network or not. The *'spiking_network'* returns important data such as the updated output-weights, which are needed for the next trial.

When the network is trained, it is tested. This happens in the loop, after the comment '% VALIDATION PART'
Here the *'spiking_network'* function is called again, but now the parameter 'FORCE' is set to 'false'.

When the network has been trained and tested, the relevant data is returned to the 'run' struct.

In the remaining part of the *'main_script_train_some'*, data from the 'run' structure can be saved.

At the end of the *'main_script_train_some'* script the network output and network spikes can be plotted for one test trial.
