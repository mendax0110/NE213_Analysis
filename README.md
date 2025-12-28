# NE213 Analysis Tools

This repository contains a collection of analysis tools developed for NE213 liquid scintillator detectors.
These tools are designed to facilitate data processing and analysis of experimental results obtained using NE213 detectors.

## Programs Included
- The ne213_analysis program: The c++ implementation of Python Script DataVisualization.py, an attempt to load faster bigger file and also have a more responsive GUI.
- NOTE: Not all Graphs are properly implemented in the c++ version yet, still needs refinement, testing and optimizations...

## Build Instructions
To build the ne213_analysis program, follow these steps:

1. Run the install_deps.sh script to install necessary dependencies:
   ```bash
   ./scripts/install_deps.sh
   ```
2. Create a build directory and navigate into it:
   ```bash
    mkdir build
    cd build
    ```
3. Run CMake to configure the build system:
    ```bash
    cmake ..
    ```
4. Compile the program using cmake
    ```bash
    cmake --build .
    ```
5. After successful compilation, you can run the ne213_analysis program from the build directory:
    ```bash
    ./ne213_analysis /PATH/TO/YOUR/DATAFILE.txt
    ```
6. For comparision or better results for now, you can also run the original Python script:
    ```bash
    python3 python/DataVisualization.py # NOTE you must manually adapt the filepath in the python script
    ```
   
## More infos
For more detailed infos, please check the papers folder, there are two papers explaining the analysis methods, calculations and so on...