# multiple-vesicles
This project provides the source code for the paper, *Numerical simulation of multiple vesicles with a N-phase phase-field system*. The code is used to generate the results shown in Section 4.15 (figure 10 and 11) of the paper.

## Author
Yutong Wu, Zecheng Qiu, Juxiang Yang* \
\*Corresponding Author

## File Descriptions
* `*.cpp` & `*.h`: C++ source and header files for the calculation.
* `show.m`: MATLAB script to visualize the final state of the vesicles.
* `draw_area.m`:  MATLAB script to plot the change in vesicle area over time.
* `draw_mass.m`: MATLAB script to plot the change in mass over time.

## Compilation and Execution
### MacOS/Linux

#### Compiling and running the calculation codes
Open your terminal, and run the following commands:
1. `cd /path/to/your/codes`
2. `make all`
3. `./bnsch.out`

#### Drawing figures
Run the different `.m` files in MATLAB

#### Deleting all the results
Input `make clean` in your terminal

### Windows
We have not supported the makefile for Windows at present, please do the following in order to run the programs.
1. Create folders named `data1`, `data2`, `data3`
2. Compile all the `*.cpp` & `*.h` files
3. Run the generated `.exe` file.
4. Run the different `.m` files in MATLAB for drawing figures
