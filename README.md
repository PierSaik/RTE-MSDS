# RTE - Machine Stable Domain Search

# Objectives 
The package provides the python code to create the simulation stable domain of a generator for the dynamic simulator Eurostag as well as several others script for analysing the computed domains and industrializing the search on multiple mahcines

# Requirements
Python 2.7 needs to be installed with the following packages  :
- scipy (version 0.12.0 or above for the ConvexHull function)
- numpy
- matplotlib
- argparse
- multiprocessing
- subprocess

Eurostag 5.1.1 with the EUROSTAG environment variable correctly set up. The python API of eurostag must be present as well. 

The scripts were tested on a Windows 7 OS and should work on UNIX if Eurostag is correctly installed and the python API is available.

*Optionnal* : a .off viewer might be useful to vizualise the domains in 3D such as [Mesh Viewer](http://mview.sourceforge.net/)

## Installing python 
[Anaconda](https://www.continuum.io/downloads) was used to developp and test the scripts and is thus recommened on windows. All the packages are already included. 

The other option is to manually install python 2.7 alongside [pip](https://pip.pypa.io/en/stable/installing/). You can then install each package manually with the following command : ">pip install 'name of package'"


# Installation & Usage
In this section all scripts will be describe along with a small use case for some. The scripts are divided in three categories, the core scripts, the analysers scripts, and the industrialization scripts. 

*Note* : the use cases need Eurostag file that are not provided in this repository. However the list of the required files is given for each case inside their folders. 

## Core scripts
The core scripts are where all the calculations are done to compute the stable simulation domain for a single generator. It is composed of two scripts :

- MachineStableDomainSearch.py : main script, calls the simulator.
- CreateMachineDomain.py       : subscript used to create the domain.

### Installation 
Copy MachineStableDomainSearch.py, CreateMachineDomain.py, dict.csv and the Example1 Folder in a directory. To get the list of inputs, go to the MachineStableDomainSearch.py directory in a console and type : ">python MachineStableDomainSearch.py -h". You must be in the directory of a generator to use this script : see Usage bellow.

### Usage 
Follow the sequence to execute the script on the Example1 generator :

1. Open a console and move to the Example1 directory : "cd ~\Example1"
2. Execute the script with the desired parameters, for example : ">python ..\MachineStableDomainSearch.py -v 0 -i 0".

*Note* : to suppress console outputs, type : ">python ..\MachineStableDomainSearch.py -v 0 -i 0 > outputs.txt". The outpurs will be stored in the outputs.txt file.

3. After it is done, you will have .off files for 3D vizualisation as well as the "ampl_generators_domain.txt" describing the domain as linear constraints.

## Analysers scripts
The analysers scripts are used to check whether the computed domains are correct by plotting UQ or PQ planes. Some analysers are also used to run analysis on power network situations. There are three analysers that will be describes separatly :

- CreatePlaneDomain.py : analyses the computed domains of one generators
- ProjectionDomain.py  : analyses the working points against the stable simulation domain of power network situations
- EurostagState.py and EurostagState_v2.py : analyses the status of a eurostag dynamic simulaiton

### CreatePlaneDomain.py

This script analyses the computed domains of one generators by comparing them betwen each other and also against a contractual or reference diagram provided by the producers or constructors if available. In the provided example, no reference diagram is available.

#### Installation 
Copy CreatePlaneDomain.py, dict.csv and Example1 in a directory. To get the list of inputs, go to the CreatePlaneDomain.py directory in a console and type : ">python CreatePlaneDomain.py -h". You must be in the directory of a generator to use this script and the stable simulation domain must have been computed before : see Usage bellow.

#### Usage 
Follow the sequence to execute the script on the Example1 generator :

1. Open a console and move to the Example1 directory : "cd ~\Example1"
2. Execute the script with the desired parameters, for example : ">python ..\MachineStableDomainSearch.py -p DiagUQ-Example1.csv -t "ZFN RPT RTE" -d 0 -s 1 -c 1 -a 0 -r 0"
3. Once done the results are saved in a .png file (Diag_ref.png if the reference is chosen, Diag_multiple.png if not)

*Note* : Even though a DiagUQ-Example1.csv is given, this file is empty, the parameters -r should stay at 0.


### ProjectionDomain.py

This script analyses all the power network situations that are present in the directory it is placed. The network situations must start by "itesla_projector_". The scripts checks if the working points of the generators in a network situation are inside the stable simulation domain. If not, it create a PQ and UQ planes which shows the distance to the domain. The results are all stored in "domain_to_network_dist_PQ.txt" and "domain_to_network_dist_UQ.txt" for all the network situations. 

*Note* : No example is provided on git hub as network sutation are classified informations.

### EurostagState.py 
This script analyses a specific eurostag simulation or all the eurostag simulation which folder starts with "itesla_eurostag_stabilization_" if no input is given. When executed, it tells wheter the initialisation of the eurostag simulation was stable or unstable. If it is unstable, a list of the unstable generators is given in outputs.

*Note* : The EurostagState_v2.py is another version that does not use argparse for cooperation issues. Here only EurostagState.py is used as it has more comprehensive inputs

#### Installation
Copy EurostagState.py and Example1 in a directory. To get the list of inputs, go to the EurostagState.py directory in a console and type : ">python EurostagState.py -h".

#### Usage
*Note* : Here, we will only use the script on one situation

1. Open a console and move to the directory where EurostagState.py was copied
2. Execute the script with the desired parameters in our case : "> python EurostagState.py -p Example1"
3. The outputs are 0 if the situation is stable, 1 if it isn't. If it's not stable a file "ecart_groupe.csv" is created with the information of the unstable generators.

## Industrialization scripts
As we have seen, the MachineStableDomainSearch.py looks for the stable domain of a single machine. These scripts aim to search the domains of multiple machines :

- IndusMachineDomain.py : main script of industrialization
- CreatedtaFile.py      : subscript that create a .dta file for IndusMachineDomain.py
- MachineDomain.bat     : bash script calling IndusMachineDomain.py for more comprehensive inputs

### IndusMachineDomain.py and MachineDomain.bat
These scripts launch the search of the stable simulation domain on several machine listed in "generators.txt" using a .dta file as base. The MachineDomain.bat is only used as a more comprehensive interface.

#### Installation
Copy MachineStableDomainSearch.py, CreateMachineDomaine.py, IndusMachineDoamin.py, MachineDomain.bat, dict.csv inside the IndusExample directory.

#### Usage
1. Double clic on MachineDomain.bat
2. Input the desired files
3. Once the execution is done the results for each machine is stored in their folders, but all the results are concatenated in the "ampl_generators_domain.txt" file in the working directory.

### CreatedtaFile.py 
This script is used to create a .dta and generators.txt file to be used in IndusMachineDomain.py from several network situation, keeping only unique occurence of generators. It requires to be placed in a folder containing .dta files and their associated dictionnaries. The results are stored in a new folder called "rez" containing a .dta file with all the information of generators, a new dictionnary and a "generators.txt" file ready to be used in IndusMachineDomain.py.

*Note* : No example is given in the repository

# License
This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Author
Pierre Saikaly (saikalypierre@gmail.com)

