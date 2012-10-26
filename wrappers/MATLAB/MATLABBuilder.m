
%The path to the main folder of the CoolProp source
path_to_src = '../../CoolProp/'

%All the include folders we need
include_string = [' -I',path_to_src,' -I',[path_to_src,'purefluids'],' -I',[path_to_src,'pseudopurefluids']]


if isempty(strfind(mexext(),'32'))
    mexopts_string = ' '
else
    mexopts_string = ' -f mexopts_w32.bat -DCONVENTION=__cdecl '
end
 
%List of files to be compiled to object files
pure_fluids = dir([path_to_src,'purefluids/*.cpp']);
pure_fluids = cellfun(@(x) fullfile(path_to_src, 'purefluids', x), {pure_fluids.name}, 'uniformoutput', false)';
ppure_fluids = dir([path_to_src,'pseudopurefluids/*.cpp']);
ppure_fluids = cellfun(@(x) fullfile(path_to_src, 'pseudopurefluids', x), {ppure_fluids.name}, 'uniformoutput', false)';
main_files = {'CoolProp.cpp','Brine.cpp','CoolPropTools.cpp','FluidClass.cpp','Helmholtz.cpp','PengRobinson.cpp','REFPROP.cpp','Solvers.cpp','Ice.cpp','HumidAirProp.cpp'}';

%Append path to source to the list of the CoolProp main files
for i=1:size(main_files,1)
    main_files{i,1} = [path_to_src,main_files{i,1}];
end
    
files=[pure_fluids; ppure_fluids; main_files];
o_files = '';
cpp_files = '';

for i=1:size(files,1)
	file = files{i,1};
	o_file = strrep(file,'.cpp','.obj');
	o_files = [o_files, ' ', o_file];
	cpp_files = [cpp_files, ' ',file];
    eval(['mex -v -c', include_string,mexopts_string,' -DCOOLPROP_LIB -outdir . ',file])
end

%Build the MEX files
eval(['mex -v ', include_string,mexopts_string,' Props.c *.obj'])
eval(['mex -v ', include_string,mexopts_string,' HAProps.c *.obj'])

%Clean up - delete the obj files
delete('*.obj')

%Quit MATLAB
quit
