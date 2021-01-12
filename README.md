# Fourmilab Floating Point Benchmarks

There are many disadvantages to being a balding geezer. In 
compensation, if you've managed to live through the second half of the 
twentieth century and been involved in computing, there's bearing 
personal witness to what happens when a technological transition goes 
into full-tilt exponential blow-off mode.  I'm talking about Moore's 
Law (actually, more of an observation than a law, since it's predicated 
on certain physical principles and can't go on forever)—computing power 
available at constant cost doubling every 18 months or so.  I've not 
only seen this happen,  I've—er—profited from it; had the 80286-based 
IBM PC/AT and its competitors not appeared when they did, Autodesk 
would have been stillborn as too early to market or drowned out by 
competitors as we arrived too late.

When Moore's Law (or whatever) is directly connected to your career and 
your bank account, it's nice to have a little thermometer you can use 
to see how it's going as the years roll by.  This repository contains 
two benchmarks I've used to evaluate computer performance ever
since 1980.  They focus on things which matter to me—floating point 
computation speed, evaluation of trigonometric functions, and matrix 
algebra.  If you're interested in text searching or database retrieval 
speed, you should run screaming from these benchmarks.  Hey, they work 
for me.

The original floating point benchmark, which is based upon an optical 
design ray tracing program I wrote in BASIC in December 1980, has been 
ported to many different programming languages, and may be used to 
evaluate the performance of these languages and different 
implementations of them.

Please see the following HTML documents for details and results.

* [fbench: Trigonometry Function Benchmark](Web_pages/fbench.html)
* [ffbench: Fast Fourier Transform Benchmark](Web_pages/ffbench.html)

All of this software is licensed under the Creative Commons 
Attribution-ShareAlike license.  Please see LICENSE.md in this 
repository for details.
