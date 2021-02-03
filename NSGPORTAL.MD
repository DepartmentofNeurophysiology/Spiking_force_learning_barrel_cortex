# Tutorial Neuroscience Gateway Portal (NSG)

For more computation time you can consult the [Neuroscience Gateway Portal(NSG)](https://nsgprod.sdsc.edu:8443/portal2/login!input.action).

## Registration

To use NSG you need to be registered, which can be done through this [form](https://www.nsgportal.org/gest/reg.php).

## How to run a job on NGS

NSG has made a step-by-step guide how to use the portal interface (there have been made some small changes to the interface but the rest is still the same):
[How to use the NSG portal](http://www.nsgportal.org/HowtousetheNSGportal.pdf).

## Preparing the input file

To be able to run files on the NSG portal you need to properly create an input zip file.
NSG needs the input zip file to have a specific format. There should be only one top level directory that you compress and upload as the zip file.

Linux:
```
mkdir <model directory>
cd <model directory>
cp files/sub directories here
cd ..
zip -r input <model directory>
```

Windows:
- Make a new file (model directory)
- Go into that file and copy the files/sub directories you want to run there
- Go back up one directory
- Right click on the model directory, in the drop-down menu, choose "Send to" and then click "Compressed (zipped) folder.

## Useful tips
These are some useful tips when running this project on the NSG portal

- Downloading and unpacking the test results of a run can be time consuming. So if you want to know a specific output of the network (such as the accuracy), make sure this result is printed in the command window. Then you only have to download the command window output of the run (STDOUT.txt).
- If you want to run a network with different settings (which cannot be included in the parametersweep) it is easier to create an input file with multiple executable scripts instead of creating and uploading a new input file for each setting. The script that is to be executed can be chosen as one of the input parameters in the NSG portal.
- Use the *Run and Time* function in Matlab to get an indication of how long it will take to run a particular simulation.
- It is possible to see how far a run is, by going to the *status* of a task that is still running and by clicking the *list* of intermediate results. 
