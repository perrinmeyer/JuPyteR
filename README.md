Test JuPyteR IJulia Notebook


https://mybinder.org/v2/gh/perrinmeyer/JuPyteR/29d130aeaf712d89d8926e15b78542b6dcef1cd4

click top file:  psmfirstnotebook.ipynb

click menu: Cell -> All Output -> clear

(all plots should go away)

Edit f2,H2,b,a,h,H= cp10hH(1000,0.5,6,96e3,2^12);

to -6dB

f2,H2,b,a,h,H= cp10hH(1000,0.5,-6,96e3,2^12);

click menu: Cell -> Run All

This will recalculate Julia file, you should see a negative 2nd order parametric
It might take up to 10-15 seconds to recalculate...

cp10hH(freq,Q,+/-dB,...)

if you change the center frequencey to 8000 and Cell->Clear All output ; Cell > Run All, you should notice that
the impulse response now decays much more quickly...  Math!



