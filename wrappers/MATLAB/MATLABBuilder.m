
%The path to the main folder of the CoolProp source
path_to_src = '../../CoolProp/'

%Get the command line options passed to script
% args = varargin;

% if sum(strcmp(args,'--octave'))>0
%     is_octave = true
% else
%     is_octave = false
% end

is_octave = false

include_string = [' -I',path_to_src,' -I',[path_to_src,'purefluids'],' -I',[path_to_src,'pseudopurefluids']]

%List of files to be compiled to object files
pure_fluids = dir([path_to_src,'purefluids/*.cpp']);
pure_fluids = cellfun(@(x) fullfile(path_to_src, 'purefluids', x), {pure_fluids.name}, 'uniformoutput', false)';
ppure_fluids = dir([path_to_src,'pseudopurefluids/*.cpp']);
ppure_fluids = cellfun(@(x) fullfile(path_to_src, 'pseudopurefluids', x), {ppure_fluids.name}, 'uniformoutput', false)';
main_files = {'CoolProp.cpp','Brine.cpp','CoolPropTools.cpp','FluidClass.cpp','Helmholtz.cpp','PengRobinson.cpp','REFPROP.cpp','Solvers.cpp'}';

disp(main_files)
%Append path to source to the list of the CoolProp main files
for i=1:size(main_files,1)
    main_files{i,1} = [path_to_src,main_files{i,1}];
end
    
files=[pure_fluids; ppure_fluids; main_files];
o_files = '';
cpp_files = '';

for i=1:size(files,1)
	file = files{i,1};
	o_file = strrep(file,'.cpp','.o');
	o_files = [o_files, ' ', o_file];
	cpp_files = [cpp_files, ' ',file];
    if is_octave
        system(['mkoctfile -v -c', include_string,' -o ',o_file,' ',file])
    else
        eval(['mex -v -c', include_string,' -outdir . ',file])
    end
end

if is_octave
    system(['mkoctfile -v --mex ', include_string,' Props.c',o_files])
else
    eval(['mex -v ', include_string,' Props.c *.obj'])
end

% % Clean up - remove the object files
% for i=1:size(files,1)
% 	file = files{i,1};
% 	o_file = strrep(file,'.cpp','.obj');
% 	unlink(o_file);
% 	disp(['deleting the file ',o_file]);
% end
% if ~is_octave
%     unlink('CoolProp.exp')
%     unlink('CoolProp.lib')
%     unlink('CoolProp_wrap.o')
%     unlink('CoolProp_wrap.cpp')
% end