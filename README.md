# multiple-vesicles
This project provides the source code for the paper, *Numerical simulation of multiple vesicles with a N-phase phase-field system*. The code is used to generate the results shown in Section 4.3.5 (figure 11 and 12) of the paper.

## Author
Yutong Wu, Zecheng Qiu, Junxiang Yang* \
\*Corresponding Author

## File Descriptions
* `*.cpp` & `*.h`: C++ source and header files for the calculation.
* `show.m`: MATLAB script to visualize the final state of the vesicles.
* `draw_area.m`:  MATLAB script to plot the change in vesicle area over time.
* `draw_mass.m`: MATLAB script to plot the change in mass over time.
* `draw_energy.m`: MATLAB script to plot the change in energy over time.
* `makefile`:  For macOS, Linux, WSL (Windows) systems, automates the compilation and cleaning process.
* `CITATION.cff`: Standard citation format file for this program.

## Compilation and Execution
### MacOS/Linux/WSL (Windows)

#### Compiling and running the calculation codes
Open your terminal, and run the following commands:
1. `cd /path/to/your/codes`
2. `make all`
3. `./bnsch.out`

#### Drawing figures
Run the different `.m` files in MATLAB.

#### Deleting all the results
Input `make delete` in your terminal.

Warning: Executing this command will **permanently remove** all generated data files (e.g., `.m` files containing calculated results) and figure files (e.g., `.eps` files) from the output directory. This operation is irreversible. Please ensure you have backed up any critical results before proceeding.

#### Plot figures & graphs
* For vesicles visualization, run `show.m` in MATLAB.
* For line chart showing area over time, run `draw_area.m` in MATLAB.
* For line chart showing mass over time, run `draw_mass.m` in MATLAB.
* For line chart showing energy over time, run `draw_energy.m` in MATLAB.

## Citation
If you use this code for your research, please cite our paper:

Yutong Wu, Zecheng Qiu, Junxiang Yang. *Numerical simulation of multiple vesicles with a N-phase phase-field system*. Computer Physics Communications (2026).

**DOI**: [10.1016/j.cpc.2026.110053](https://doi.org/10.1016/j.cpc.2026.110053)

If you use BibTeX, please use the following entry:

```bibtex
@Article{Wu2026,
  author    = {Wu, Yutong and Qiu, Zecheng and Yang, Junxiang},
  journal   = {Comput. Phys. Commun.},
  title     = {A three-dimensional multi-phase-field vesicles model and its practical finite difference solver},
  year      = {2026},
  issn      = {0010-4655},
  pages     = {110053},
  doi       = {10.1016/j.cpc.2026.110053},
  fjournal  = {Computer Physics Communications},
  publisher = {Elsevier BV},
}
```
