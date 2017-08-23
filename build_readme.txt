If xsd-files are changed, create .cxx and .hxx files with the following command
# xsd cxx-tree --generate-serialization ***.xsd

need to install xerces-c (http://xerces.apache.org/) and xsd (http://www.codesynthesis.com/products/xsd/)

In order to build program:
cd ./build; rm CMakeCache.txt; cmake -D CMAKE_BUILD_TYPE=Release ..; make; cd ..;
