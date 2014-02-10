copy ..\..\..\wrappers\MATLAB\example.m .
copy ..\..\..\wrappers\MATLAB\*rops.c .
matlab -wait -nojvm -nodesktop -r "MATLABBuilder; quit"
matlab -nojvm -nodesktop -r "Example; quit" -logfile Output.txt
erase *Props.c