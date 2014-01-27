copy ..\..\..\wrappers\MATLAB\example.m
copy ..\..\..\wrappers\MATLAB\*rops.c
matlab -nojvm -nodesktop -r MATLABBuilder
matlab -nojvm -nodesktop -r "Example; quit" -logfile Output.txt
erase *Props.c