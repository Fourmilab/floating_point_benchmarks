
#   Reformat output from ALGOL 60 floating point benchmark and test

#   The "marst" ALGOL 60 to C translator used with the
#   floating point benchmark has only rudimentary
#   output formatting facililties for floating point
#   numbers.  They are simply printed with a C "%.12g"
#   format, which results in some values having too many
#   and other two few decimal places, and one being output
#   in scientific notation.  This makes if difficult to
#   compare the output from the program to that of other
#   implementations.

#   This little Perl program post-processes the output of
#   the benchmark into the standard format used by
#   implementations in other languages and checks the result
#   against a master copy of the correct results.

    use strict;
    use warnings;

    my @refarr = (  	    	    # Reference results.  These happen to
    	    	    	    	    # be derived from a run on Microsoft 
    	    	    	    	    # Quick BASIC on the IBM PC/AT.

            '   Marginal ray          47.09479120920   0.04178472683',
            '   Paraxial ray          47.08372160249   0.04177864821',
            'Longitudinal spherical aberration:        -0.01106960671',
            '    (Maximum permissible):                 0.05306749907',
            'Offense against sine condition (coma):     0.00008954761',
            '    (Maximum permissible):                 0.00250000000',
            'Axial chromatic aberration:                0.00448229032',
            '    (Maximum permissible):                 0.05306749907'
    );

    my ($l, $i, $errors);
    my @outarr;
    
    #	Reformat the object distance and slope angle for the
    #	marginal and paraxial rays.
    
    for ($i = 0; $i < 2; $i++) {
    	$l = <>;
	$l =~ m/(^[\sA-Za-z]+)([^\s]+)\s+([^\s]+)/ || die("Cannot parse $l");
	my ($label, $od, $sa) = ($1, $2, $3);
	$label =~ s/\s+$//;
    	push(@outarr, sprintf("%15s   %21.11f  %14.11f", $label, $od, $sa));
    }

    #	Reformat the design analysis.
    
    for ($i = 0; $i < 6; $i++) {
    	$l = <>;
	$l =~ m/(^[\s\(\):A-Za-z]+)([^\s]+)/ || die("Cannot parse $l");
	my ($label, $ev) = ($1, $2);
	$label =~ s/\s+$//;
	push(@outarr, sprintf("%-39s %16.11f", $label, $ev));
    }
    
    #	Print and compare the reformatted output with the
    #	expected output, reporting any errors.
    
    for ($i = $errors = 0; $i <= $#outarr; $i++) {
    	print("$outarr[$i]\n");
	if ($outarr[$i] ne $refarr[$i]) {
	    printf("Error on line %d.  Should be \"%s\"\n", $i + 1, $refarr[$i]);
	    $errors++;
	}
    }
    
    #	Print error summary.

    if ($errors > 0) {
       printf("\n%d error%s in results.  This is VERY SERIOUS.\n",
          $errors, $errors > 1 ? "s" : "");
    } else {
       printf("\nNo errors in results.\n");
    }
