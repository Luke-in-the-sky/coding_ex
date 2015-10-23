-----
This is a ROOT macro. For more info, visit https://root.cern.ch
-----

This Macro analyzes the noise levels in (voltage) measurements performed under 
different conditions (temperature and electrical current).

We have 11 txt-input files. We want to monitor if some values of applied current/temperature 
have an impact on the noise level in the measured output current (2nd column), the measured 
temperature (4th column) or the measured voltage (5th column)
Each input file has data stored in a tabular fashion, 
but not all colums are of float type. This means we will have to read the files
line by line, select what we need and save it to a new, float-only TNtuple object
for analysis and plotting.
