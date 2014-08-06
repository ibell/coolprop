EXPORTS="['_main','_F2K','_Props1SI','_PropsSI','_get_global_param_string','_HAProps']"
em++ -O2 ../../CoolProp/*.cpp -I../../CoolProp -c  -s EXPORTED_FUNCTIONS=${EXPORTS} -s DISABLE_EXCEPTION_CATCHING=0
em++ -O2 -o coolprop.js ../../CoolProp/CoolPropDLL.cpp *.o -I../../CoolProp -DEXTERNC  -s EXPORTED_FUNCTIONS=${EXPORTS} -s DISABLE_EXCEPTION_CATCHING=0

# Can only use compression with HTML :(
# --compression ~/Code/emscripten/third_party/lzma.js/lzma-native,/home/xubuntu/Code/emscripten/third_party/lzma.js/lzma-decoder.js,LZMA.decompress

# Using closure compiler to compress javascript file
#java -jar compiler.jar --js precursor.js --js_output_file coolprop.js --compilation_level ADVANCED_OPTIMIZATIONS --language_in ECMASCRIPT5

rm *.o
