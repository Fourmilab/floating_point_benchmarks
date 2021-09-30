#!/usr/bin/env perl

#        John Walker's Floating Point Benchmark, derived from...
#
#	Marinchip Interactive Lens Design System
#
#				     John Walker   December 1980
#
#	By John Walker
#	   http://www.fourmilab.ch/
#
#	This program may be used, distributed, and modified freely as
#	long as the origin information is preserved.
#
#	This is a complete optical design raytracing algorithm,
#	stripped of its user interface and recast into Perl.
#	It not only determines execution speed on an extremely
#	floating point (including trig function) intensive
#	real-world application, it checks accuracy on an algorithm
#	that is exquisitely sensitive to errors.  The performance of
#	this program is typically far more sensitive to changes in
#	the efficiency of the trigonometric library routines than the
#	average floating point program.

use strict;
use warnings;
use PDL::LiteF;
use PDL::Math;

my $ITERATIONS = 1000;

#   Local variables

my ($aberr_lspher, $aberr_osc,
    $aberr_lchrom, $max_lspher, $max_osc, $max_lchrom, $radius_of_curvature,
    $object_distance, $ray_height, $axis_slope_angle, $from_index, $to_index);

my @s;  	    	    	    # Design being traced
my @od_sa;	    	    	    # Object distance and slope angle

my @outarr;     	    	    # Computed output of program goes here

my $niter = $ITERATIONS;	    # Iteration counter

my $spectral_line = pdl(   	    # Spectral lines in which we trace
  # wavelengths of Fraunhofer A-F spectral lines
  undef,
  7621.0, 	    # A
  6869.955,	    # B
  6562.816,	    # C
  5895.944,	    # D
  5269.557,	    # E
  4861.344,	    # F
  4340.477,	    # G'
  3968.494	    # H
);
($radius_of_curvature, $object_distance, $ray_height, $axis_slope_angle, $from_index, $to_index) = map pdl(0), 1..6;

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

#	The test case used in this program is the design for a 4 inch f/12
#	achromatic telescope objective used as the example in Wyld's
#	classic work on ray tracing by hand, given in Amateur Telescope
#	Making, Volume 3, p597.

my $testcase = pdl(
        # radius of curvature, refractive index, Abbe number, axial distance
        [ 27.05, 1.5137, 63.6, 0.52 ], # crown
        [ -16.68, 1, 0, 0.138 ],
        [ -16.68, 1.6164, 36.7, 0.38 ], # flint
        [ -78.1, 1, 0, 0 ]
);

use Inline Pdlpp => <<'EOF';
pp_def('trace_line',
  Pars => '
    surface(elt=4,s); spectral_lines(l); int paraxial(); indx line(); ray_height();
    [o] object_distance(s); [o] slope_angle(s);',
  GenericTypes => ['D'],
  Code => <<'EOCODE',
  /*
    Inputs:
  radius_of_curvature	  Radius of curvature of surface
			  being crossed.  If 0, surface is plane.
  object_distance  	  Distance of object focus from
			  lens vertex.	If 0, incoming rays are parallel and
			  the following must be specified:
  ray_height		  Height of ray from axis.  Only
			  relevant if $object_distance == 0
  axis_slope_angle 	  Angle incoming ray makes with axis at intercept
  from_index		  Refractive index of medium being left
  to_index		  Refractive index of medium being entered.
    Outputs:
  object_distance	  Distance from vertex to object focus after refraction.
  axis_slope_angle	  Angle incoming ray makes with axis
			  at intercept after refraction.
  */
  $GENERIC() o_d = 0, s_a, r_h = $ray_height(), from_index = 1;
  loop(s) %{
    /* Perform ray trace in specific spectral line */
    $GENERIC() radius_of_curvature = $surface(elt=>0, s=>s),
      to_index = $surface(elt=>1, s=>s),
      abbe_number = $surface(elt=>2, s=>s),
      axial_distance = $surface(elt=>3, s=>s);
    PDL_Indx which_line = $line();
    if (to_index > 1) {
      to_index += (
         ($spectral_lines(l=>4) - $spectral_lines(l=>which_line)) /
         ($spectral_lines(l=>3) - $spectral_lines(l=>6))
        ) *
        ((to_index - 1) / abbe_number);
    }
    $GENERIC() iang_sin /* Incidence angle sin */,
       rang_sin /* Refraction angle sin */;
    if ($paraxial()) {
      if (radius_of_curvature != 0) {
        if (o_d == 0) {
          s_a = 0;
          iang_sin = r_h / radius_of_curvature;
        } else {
          iang_sin = ((o_d -
            radius_of_curvature) / radius_of_curvature) *
            s_a;
        }
        rang_sin = (from_index / to_index) *
          iang_sin;
        $GENERIC() old_axis_slope_angle = s_a;
        s_a += iang_sin - rang_sin;
        if (o_d != 0)
          r_h = o_d * old_axis_slope_angle;
        o_d = r_h / s_a;
      } else {
        o_d *= to_index / from_index;
        s_a *= from_index / to_index;
      }
    } else {
      if (radius_of_curvature != 0) {
        if (o_d == 0) {
          s_a = 0;
          iang_sin = r_h / radius_of_curvature;
        } else {
          iang_sin = ((o_d -
            radius_of_curvature) / radius_of_curvature) *
            sin(s_a);
        }
        $GENERIC() iang = asin(iang_sin); /* Incidence angle */
        rang_sin = (from_index / to_index) *
          iang_sin;
        $GENERIC() old_axis_slope_angle = s_a;
        s_a += iang - asin(rang_sin);
        $GENERIC() sagitta = sin((old_axis_slope_angle + iang) / 2);
        sagitta = 2 * radius_of_curvature * sagitta * sagitta;
        o_d = ((radius_of_curvature * sin(
          old_axis_slope_angle + iang)) /
          tan(s_a)) + sagitta;
      } else {
        $GENERIC() rang = -asin((from_index / to_index) *
          sin(s_a)); /* Refraction angle */
        o_d *= ((to_index *
          cos(-rang)) / (from_index *
          cos(s_a)));
        s_a = -rang;
      }
    }
    from_index = to_index;
    o_d -= axial_distance;
    $object_distance() = o_d;
    $slope_angle() = s_a;
  %}
EOCODE
);
EOF

# Process the number of iterations argument, if one is supplied.

if ($#ARGV >= 0) {
   $niter = shift(@ARGV);
   if (($niter !~ m/^\d+$/) || ($niter < 1)) {
      print << "EOD";
This is John Walker's floating point accuracy and
performance benchmark program.  You call it with

perl fbench.pl <itercount>

where <itercount> is the number of iterations
to be executed.  Archival timings should be made
with the iteration count set so that roughly five
minutes of execution is timed.
EOD
      exit;
   }
}

# Load test case into working array

my $clear_aperture = 4;

my $nik = $niter / 1000;
print << "EOD";
Ready to begin John Walker's floating point accuracy
and performance benchmark.  $niter iterations will be made.

Measured run time in seconds should be divided by $nik
to normalise for reporting results.  For archival results,
adjust iteration count so the benchmark runs about five minutes.

EOD
use Time::HiRes qw(gettimeofday tv_interval);
my @t=gettimeofday();

# Perform ray trace the specified number of times.

my ($od_fline, $od_cline);

for (my $itercount = 0; $itercount < $niter; $itercount++) {
  my @inputs = (
    $testcase,
    $spectral_line,
    pdl(0, 1, 0, 0), # paraxial
    pdl(4, 4, 3, 6), # spectral line - main trace in D light, marginal in C,F
    pdl($clear_aperture / 2), # ray height, threads so no need to repeat
  );
  my ($od, $sa) = PDL::trace_line(@inputs);
  my $pdl_od_sa = pdl($od, $sa)->slice("(3)")->transpose; # slice as only last col is of interest
  @od_sa = @{ $pdl_od_sa->slice(",0:1")->unpdl };
  ($od_cline, $od_fline) = @{ $pdl_od_sa->slice("(0),2:3")->unpdl };
  $aberr_lspher = $od_sa[1][0] - $od_sa[0][0];
  $aberr_osc = 1 - ($od_sa[1][0] * $od_sa[1][1]) /
     (sin($od_sa[0][1]) * $od_sa[0][0]);
  $aberr_lchrom = $od_fline - $od_cline;
  $max_lspher = sin($od_sa[0][1]);
  # D light
  $max_lspher = 0.0000926 / ($max_lspher * $max_lspher);
  $max_osc = 0.0025;
  $max_lchrom = $max_lspher;
}

my $interval = tv_interval(\@t);
print "Time taken: $interval\n";
print "Divided by $nik = ", $interval/$nik, "\n";

# Now evaluate the accuracy of the results from the last ray trace

$outarr[0] = sprintf("%15s   %21.11f  %14.11f",
   "Marginal ray", $od_sa[0][0], $od_sa[0][1]);
$outarr[1] = sprintf("%15s   %21.11f  %14.11f",
   "Paraxial ray", $od_sa[1][0], $od_sa[1][1]);
$outarr[2] = sprintf(
   "Longitudinal spherical aberration:      %16.11f",
   $aberr_lspher);
$outarr[3] = sprintf(
   "    (Maximum permissible):              %16.11f",
   $max_lspher);
$outarr[4] = sprintf(
   "Offense against sine condition (coma):  %16.11f",
   $aberr_osc);
$outarr[5] = sprintf(
   "    (Maximum permissible):              %16.11f",
   $max_osc);
$outarr[6] = sprintf(
   "Axial chromatic aberration:             %16.11f",
   $aberr_lchrom);
$outarr[7] = sprintf(
   "    (Maximum permissible):              %16.11f",
   $max_lchrom);

# Now compare the edited results with the master values from
# reference executions of this program.

my $errors = 0;
for (my $i = 0; $i < 8; $i++) {
   if ($outarr[$i] ne $refarr[$i]) {
      printf("\nError in results on line %d...\n", $i + 1);
      print("Expected:  $refarr[$i]\n");
      print("Received:  $outarr[$i]\n");
      print("(Errors)   ");
      my $k = length($refarr[$i]);
      for my $j (0..$k-1) {
         print(substr($refarr[$i], $j, 1) eq substr($outarr[$i], $j, 1) ? ' ' : '^');
         if (substr($refarr[$i], $j, 1) ne substr($outarr[$i], $j, 1)) {
            $errors++;
        }
      }
      print("\n");
   }
}
if ($errors > 0) {
   printf("\n%d error%s in results.  This is VERY SERIOUS.\n",
      $errors, $errors > 1 ? "s" : "");
} else {
   printf("\nNo errors in results.\n");
}
