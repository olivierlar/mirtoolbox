Congratulations! You have downloaded MIRtoolbox 1.2

The list of new features and bug fixes offered by this new version is shown in the ReleaseNotes text file located in the toolbox folder, and also available from our website. New versions of the toolbox will also be released at the same address:
http://www.jyu.fi/hum/laitokset/musiikki/en/research/coe/materials/mirtoolbox

Please inform me if you find any bug when using the toolbox.

**

MIRtoolbox requires Matlab version 7 and Mathworks' Signal Processing toolbox.

**

This distribution actually includes four different toolboxes:

- the MIRtoolbox itself, version 1.2,

- the Auditory toolbox, version 2, by Malcolm Slaney
(also available directly at http://cobweb.ecn.purdue.edu/~malcolm/interval/1998-010/),

- the Netlab toolbox, version 3.3, by Ian Nabney
(also available directly at http://www.ncrg.aston.ac.uk/netlab/)

- the SOM toolbox, version 2.0, by Esa Alhoniemi, Johan Himberg, Jukka Parviainen and Juha Vesanto
(also available directly at http://www.cis.hut.fi/projects/somtoolbox/)

**

MIRtoolbox license is based on GPL 2.0. As such, it can integrate codes from other GPL 2.0 projects, as long as their origins are explicitly stated.
- codes from the Music Analysis Toolbox by Elias Pampalk (2004), related to the computation of Terhardt outer ear modeling, Bark band decomposition and masking effects.

- openbdf and readbdf script by T.S. Lorig to read BDF files, based on openedf and readedf by Alois Schloegl. 

**

To install MIRtoolbox in your Matlab environment, copy all the toolboxes folders (or only those that are not installed yet in your computer) into your Matlab "toolbox" folder. Then add each folder in your Matlab path.

**

NOTICE: If you replace an older version of MIRtoolbox with a new one, please update your Matlab path using the following command:
rehash toolboxcache
Update also the class structure of the toolbox, either by restarting Matlab, or by typing the following command:
clear classes

**

To get an overview of the functions available in the toolbox, type:
help mirtoolbox

A complete documentation in PDF is provided in the main toolbox folder.

A short documentation for each function is available directly in Matlab using the help command. For instance: help miraudio

Examples of use of the toolbox are shown in the MIRToolboxDemos folder:
mirdemo
demo1basics
demo2timbre
demo3segmentation
demo4tempo
demo5export
demo6curves
demo7tonality
demo8classification

**

If you have any problem or request, please contact us.

Olivier Lartillot, Petri Toiviainen and Tuomas Eerola

Olivier.Lartillot@campus.jyu.fi