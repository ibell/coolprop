copy ..\..\..\wrappers\MATLAB\Example.m .
copy ..\..\..\wrappers\MATLAB\*Props*.c .
matlab -wait -nojvm -nodesktop -r "MATLABBuilder; quit" -logfile build_log.txt
matlab -wait -nojvm -nodesktop -r "Example; quit" -logfile Output.txt
erase *Props*.c