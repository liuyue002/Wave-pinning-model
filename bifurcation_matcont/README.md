These Matlab code run bifurcation analysis on the wave pinning model with MATCONT (Dhooge, Govaerts & Kuznetsov, 2003). It more or less has the same capability of XPPAUT, but its command line API makes it much easier to make reproducible plots and integrate smoothly with Matlab.

There is a complete manual [*MATCONT: Continuation toolbox for ODEs in Matlab*](http://www.staff.science.uu.nl/~kouzn101/NBA/ManualJan2018.pdf) (Govaerts et al, 2018).

There are some very helpful tutorials for MATCONT by Holmes (2012)[here](http://www.math.ubc.ca/~keshet/MCB2012/wrholmes/MCB2012.html), which my scripts build upon. Also see [*Local Perturbation Analysis in MatCont: User’s Guide & Tutorials*](http://www.math.ubc.ca/~keshet/Papers/LPA_BJTools_Supp2.pdf) by Holmes, Mata & Edelstein-Keshet (2014).

The original authors' website is no longer working. The software can be downloaded from [SourceForge](https://sourceforge.net/projects/matcont/files/matcont/matcont6p11/). My codes used version 6.11.

Part of MATCONT is written in C, so make sure you already have a C compiler. During installation, run `mex -setup` in Matlab, then choose a compiler. Sometimes you get errors like "<Filename.c> not found". In these cases you need to manually compile the C codes by running `mex Filename`.

Each txt file here describe a model. To import a model to MATCONT, use the graphical interface by running `matcont`, then Select->System>New. The first three lines of the txt file corresponds to Name, Coordinates, Parameters, respectively. Leave Time as "t". For computing derivatives, use symbolic if you have the toolbox. In the large text box, input the rest of the txt file, then click Ok. This will create a new file in the Systems folder.

To run my scripts, make sure both MATCONT's installation folder and the Systems folder are in the path (you will need to modify the `addpath` commands at the beginning).

I also provide some functions that are helpful for creating AUTO-style plots.
