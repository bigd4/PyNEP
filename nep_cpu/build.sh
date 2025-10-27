c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) src/pynep.cpp src/nep.cpp -o nep$(python3-config --extension-suffix) -Wno-attributes -Wno-sign-compare -Wno-unused-variable
