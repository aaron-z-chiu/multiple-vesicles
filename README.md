# multiple-vesicles
This project provides the source code for the paper, *Numerical simulation of multiple vesicles with a N-phase phase-field system*. The code is used to generate the results shown in Figure 10 (Section 4.15) of the paper.

## Author
Yutong Wu, Zecheng Qiu, Juxiang Yang* \
\*Corresponding Author

## Files
* `*.cpp` & `*.h`: Calculation codes
* `show.m`: demostrate the figure of models
* `draw_area.m`: Showing a line chart of the area with time
* `draw_mass.m`: Showing a line chart of the mass with time

## Compile and run
### MacOS/Linux

#### Compiling and running the calculation codes
Open your terminal, and run the following commands:
1. `cd [Your data path]`
2. `make all`
3. `./bnsch.out`

#### Drawing figures
Run the different `.m` files in MATLAB

#### Deleting all the results
Input `make clean` in your terminal

### Windows
We have not supported the makefile for Windows at present, please do the following in order to run the programs.
1. Make folders named `data1`, `data2`, `data3`
2. Compile all the `*.cpp` & `*.h` files
3. Run the generated `.exe` file.
4. Run the different `.m` files in MATLAB for drawing figures
