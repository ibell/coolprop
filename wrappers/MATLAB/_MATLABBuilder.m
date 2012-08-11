
%The path to the main folder of the CoolProp source
path_to_src = '../../CoolProp/'

%Get the command line options passed to script
args = argv();

if sum(strcmp(args,'--octave'))>0
    _octave = true
else
    _octave = false
end

include_string = [' -I',path_to_src,' -I',[path_to_src,'purefluids'],' -I',[path_to_src,'pseudopurefluids']]

%List of files to be compiled to object files
pure_fluids = glob([path_to_src,'purefluids/*.cpp']);
ppure_fluids = glob([path_to_src,'pseudopurefluids/*.cpp']);
main_files = {'CoolProp.cpp','Brine.cpp','CoolPropTools.cpp','FluidClass.cpp','Helmholtz.cpp','PengRobinson.cpp','REFPROP.cpp','Solvers.cpp'}';

%Append path to source to the list of the CoolProp main files
for i=1:size(main_files)(1)
    main_files{i,1} = [path_to_src,main_files{i,1}];
end
    
files=[pure_fluids; ppure_fluids; main_files];
o_files = '';
cpp_files = '';

for i=1:size(files)(1)
	file = files{i,1};
	o_file = strrep(file,'.cpp','.o');
	o_files = [o_files, ' ', o_file];
	cpp_files = [cpp_files, ' ',file];
    if _octave
        system(['mkoctfile -v -c', include_string,' -o ',o_file,' ',file])
    else
        system(['mex -v -c', include_string,' -o ',o_file,' ',file])
    end
end

if _octave
    system(['mkoctfile -v --mex ', include_string,' Props.c',o_files])
else
    system(['mex -v ', include_string,' Props.c',o_files])
end

% Clean up - remove the object files
for i=1:size(files)(1)
	file = files{i,1};
	o_file = strrep(file,'.cpp','.o');
	unlink(o_file);
	disp(['deleting the file ',o_file]);
end
if !_octave
    unlink('CoolProp.exp')
    unlink('CoolProp.lib')
    unlink('CoolProp_wrap.o')
    unlink('CoolProp_wrap.cpp')
end