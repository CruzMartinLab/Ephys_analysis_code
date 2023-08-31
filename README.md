# Ephys_analysis_code




The description below details the contents of the ‘mEPSC_analysis’ folder you have downloaded. To ensure that the code, as presently constructed, works smoothly with your data, you MUST:
1.	Name and filter your raw traces in the same format as they exist in the sample_data. See details below for naming structure/filtering of data
2.	Your folder structure must also be identical to what is provided here in the sample data:
a.	____________/cell_type/Condition_x/cell_x
i.	_____whatever your pathing is to get to this point does not matter
ii.	cell_type: PV cell, pyramidal cell, etc.
iii.	Condition_x: KO, KI, KD, control, etc
iv.	cell_x: number cell, this folder contains all the recordings from your cell

Two folders:
1.	‘MATLAB’ contains:
a.	mPSC_detection.m: primary ‘master’ code you will need to open and execute
i.	line 73-77: change function depending on if you are running the code for pyramidal cells or PV cells
ii.	line 119: you must change if your recordings are anything other than 1 minute, gap-free recordings. If you take 30s recordings, for example, change 60 to 30
b.	mPSC_detection_functions: folder containing all necessary functions
i.	abfload_fcollman.m: reads in .abf files (raw traces from clampex)
1.	documentation available online, I did not write this
ii.	detrend_y_convert.m
1.	helps correct for drift in baseline, if there is any
2.	detrends your trace and removes any signal ~6pA above baseline (potential oscillation/noise in your rig)
iii.	find_EPSCs_PVIN.m
1.	primary function that you will run if you are analyzing mEPSCS from parvalbumin interneurons (PV-INs) (fast kinetics, larger amplitudes)
iv.	find_EPSCs_PYR.m
1.	primary function that you will run if you are analyzing mEPSCs from pyramidal neurons (PYR) (slower kinetics, smaller amplitudes)
2.	‘sample_data’ contains:
a.	PV-IN: data from PV-INs
b.	PYR: data from PYR
c.	Note: the folder structure of both PV-IN and PYR are identical, so I will describe the contents of both in general:
i.	Condition_X: folder containing data (multiple 1 minute gap-free raw .abf traces) from two cells each per ‘condition’
1.	cell_x: folder containing:
a.	raw_trace (.abf)
i.	all raw traces must be listed in the following format
1.	YYYY_MM_DD_xxxx_filtered.abf
a.	YYYY_MM_DD date of recording
b.	xxxx = number of recording that day (e.g. 0001 or 0019 for the first and nineteenth recording that day, respectively)
c.	filtered: all of my raw traces are filtered in clampfit. Analyze->Filter->Lowpass (boxcar, 7 smoothing points, filter the whole trace)
d.	.abf: the natural clampex file extension
ii.	AFTER YOU RUN THE CODE, EACH cell_x FOLDER WILL ALSO CONTAIN:
1.	YYYY_MM_DD_xxxx_filtered.abf_Cell_output.xlsx: excel that has all the data for that one recording to which the name corresponds
2.	YYYY_MM_DD_xxxx_filtered.abf_Figure: figure that can be opened in matlab that has the raw trace and all the events that were identified by the code
3.	cell_x_Final_Results.xlsx: excel that combines the average amplitude, frequency, rise, and decay from each of the recordings of a given cell. Also gives the overall average (and standard deviation) of the amplitude, frequency, rise, and decay of the cell overall by averaging (or taking the standard deviation) of the 4 metrics from each of that cell’s recordings
4.	cell_x_Final_Results_w_IEI.xlsx: gives the IEIs, amplitudes, rise, and decay of each individual event from all recordings from that given cell
ii.	Condition_X_all_IEIs.xlsx: excel that contains every inter-event interval and amplitude for each event from all recordings of given cell in a given condition_x folder, organized on a cell-by-cell basis (i.e. concatenates all the data from all the recordings of a given cell)
1.	This is used for cumulative frequency distributions
iii.	Condition_X_all_Kinetics.xlsx: excel that contains the Decay and Rise for each event from all recordings of given cell in a given condition_x folder, organized on a cell-by-cell basis (i.e. concatenates all the data from all the recordings of a given cell)
1.	This is used for cumulative frequency distributions
iv.	Condition_X_Final_Results.xlsx: excel that contains the average amplitude, frequency, rise, and decay for all cells within a condition_x folder
1.	This is the main data you would plot

