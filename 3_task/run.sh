#!/bin/bash

echo Compiling...
gcc btw_openmp.cpp -o btw -lstdc++ -std=c++11 -fopenmp
echo Generation graph and solving it by python ...
python gen_graph.py > /dev/null
echo Solving BTW problem by C++ ...
./btw > /dev/null
echo "diff btwcheck_cpp btwcheck_py:"
diff btwcheck_cpp btwcheck_py
