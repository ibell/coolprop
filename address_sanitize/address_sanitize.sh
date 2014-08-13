/usr/bin/clang -fsanitize=memory -c ../CoolProp/*.cpp -I../CoolProp
/usr/bin/clang -fsanitize=memory -c ../main.cpp -I../CoolProp
/usr/bin/clang -fsanitize=memory -lstdc++ *.o
ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer-3.5 ./a.out

