Bifurcation analysis of the wave pinning model with XPPAUT, a program for dynamical systems solver by [B. Ermentrout](http://www.math.pitt.edu/~bard/xpp/xpp.html), which has built-in support for AUTO, a numerical continuation software.

Instructions:

1. Download and install XPPAUT 8.0 from http://www.math.pitt.edu/~bard/bardware/binary/latest/
    * Older versions should work, but might crash frequently
2. Run XPPAUT, i.e. `xppaut wavepin_lpa.ode`
    * On Linux, you might get an error complaing about "libc.so.6" when trying to run XPPAUT. In that case install the correct version of C Standard Library for your system.
3. Tutorials for XPPAUT can be found [here](http://www.math.pitt.edu/~bard/xpp/help/xpphelp.html). Also see the book *Simulating, Analyzing, and Animating Dynamical Systems: A Guide to XPPAUT* (Ermentrout 2002).
