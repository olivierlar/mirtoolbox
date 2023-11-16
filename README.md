NOTE: The official MIRtoolbox webpage was until November 2023 hosted at the University of Jyväskylä. Since November 14th, 2023, the GitHub repository has become public.

# mirtoolbox

MIRtoolbox is a Matlab toolbox dedicated to the extraction of musical features from audio files, including routines for statistical analysis, segmentation and clustering. MIRtoolbox integrates a user-friendly syntax that enables to easily combine low and high-level operators into complex flowcharts. The modular design of MIRtoolbox is guided by a philosophy of expertise capitalization: techniques developed for certain domains of music analysis are turned into general operators that could be used for different analytical purposes.

Each feature extraction method can accept as argument an audio file, or any preliminary result from intermediary stages of the chain of operations. Also the same syntax can be used for analyses of single audio files, batches of files, series of audio segments, multi-channel signals, etc. For that purpose, the data and methods of the toolbox are organised in an object-oriented architecture.

Memory management mechanisms allow the analysis of large-scale corpus, through the integration of automated chunk decomposition mechanisms and of distinctive processes for flowchart design and evaluation. A set of meta-functions have been designed that enable the integration of user-defined algorithms with the help of simple templates.

# authors

Olivier Lartillot, Petri Toiviainen, Pasi Saari and Tuomas Eerola were members of the Finnish Centre of Excellence in Interdisciplinary Music Research, University of Jyväskylä, Finland.
