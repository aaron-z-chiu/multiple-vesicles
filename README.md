# multiple-vesicles
This project provides the source code for the paper, *Numerical simulation of multiple vesicles with a N-phase phase-field system*. The code is used to generate the results shown in Section 4.15 (figure 10 and 11) of the paper.

## Author
Yutong Wu, Zecheng Qiu, Junxiang Yang* \
\*Corresponding Author

## File Descriptions
* `*.cpp` & `*.h`: C++ source and header files for the calculation.
* `show.m`: MATLAB script to visualize the final state of the vesicles.
* `draw_area.m`:  MATLAB script to plot the change in vesicle area over time.
* `draw_mass.m`: MATLAB script to plot the change in mass over time.
* `makefile`:  For macOS and Linux systems, automates the compilation and cleaning process.

## Compilation and Execution
### MacOS/Linux

#### Compiling and running the calculation codes
Open your terminal, and run the following commands:
1. `cd /path/to/your/codes`
2. `make all`
3. `./bnsch.out`

#### Drawing figures
Run the different `.m` files in MATLAB.

#### Deleting all the results
Input `make delete` in your terminal.

##### Special Warning
Executing this command will **permanently remove** all generated data files (e.g., `.m` files containing calculated results) and figure files (e.g., `.eps` files) from the output directory. This operation is irreversible. Please ensure you have backed up any critical results before proceeding.

### Windows
A Makefile for Windows has not been provided at this moment. Please follow the steps below to compile and run the code manually:
1. Create folders named `data1`, `data2`, `data3` in the same path as codes.
2. Compile all the `*.cpp` & `*.h` files.
3. Run the generated `.exe` file.
4. Run the different `.m` files in MATLAB for drawing figures
