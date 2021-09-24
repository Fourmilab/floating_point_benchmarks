This directory contains an absolutely minimalist port of the Perl 5
version of Fbench to Raku.  It was created simply by fixing all of the
differences between the two languages.  The Perl 5 version used the
cot() function from the Math::Trig module.  To avoid the need to deal
with Raku modules, I replaced the one place cot() is used in the code
with (1 / tan(x)) to avoid calling a function.

It was developed and tested on Rakudo 2021.09 in September 2021.
