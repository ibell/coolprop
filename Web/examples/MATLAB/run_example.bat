copy ..\..\..\wrappers\MATLAB\example.m
copy ..\..\..\wrappers\MATLAB\*rops.c
matlab -r MATLABBuilder
matlab -nojvm -nodesktop -r "Example; quit" -logfile Output.txt