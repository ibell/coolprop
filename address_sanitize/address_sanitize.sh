/usr/bin/clang -fsanitize=address -c ../CoolProp/*.cpp -I../CoolProp
/usr/bin/clang -fsanitize=address -c ../main.cpp -I../CoolProp
/usr/bin/clang -fsanitize=address -lstdc++ *.o
ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.5 ./a.out

