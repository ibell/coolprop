
%The path to the main folder of the CoolProp source
path_to_src = '../../CoolProp/'

%All the include folders we need
include_string = [' -I',path_to_src];

mexopts_string = ' '
 
%List of files to be compiled to object files
main_files = dir([path_to_src,'*.cpp']);
main_files = cellfun(@(x) x, {main_files.name}, 'uniformoutput', false)';

%Append path to source to the list of the CoolProp main files
for i=1:size(main_files,1)
    main_files{i,1} = [path_to_src,main_files{i,1}];
end
    
files=[main_files];
o_files = '';
cpp_files = '';

for i=1:size(files,1)
	file = files{i,1};
	o_file = strrep(file,'.cpp','.o');
	o_files = [o_files, ' ', o_file];
	cpp_files = [cpp_files, ' ',file];
    eval(['mex -v -c', include_string,mexopts_string,' -DCOOLPROP_LIB -outdir . ',file])
end

%Build the MEX files
eval(['mex -v ', include_string,mexopts_string,' Props.c *.o'])
eval(['mex -v ', include_string,mexopts_string,' HAProps.c *.o'])

%Clean up - delete the obj files
delete('*.o')

%Quit MATLAB
quit
