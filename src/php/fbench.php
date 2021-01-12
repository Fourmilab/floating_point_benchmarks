<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<title>PHP Floating Point Benchmark</title>
<link rel="stylesheet" href="https://www.fourmilab.ch/documents/styles/standard_screen.css"
    type="text/css" />
<meta name="keywords" content="PHP, floating, point, benchmark, Fourmilab" />
<meta name="description" content="PHP Floating Point Benchmark" />
<meta name="author" content="John Walker"  />
<meta name="robots" content="index,nofollow" />
<style type="text/css">
    blockquote.rights {
        font-family: sans-serif;
        font-size: smaller;
    }

    h1, h2, h3 {
        margin-bottom: 0px;
    }

    p.cn {
        text-align: center;
        margin-top: 0px;
    }

    p {
        text-align: justify;
    }

    table.footer {
        width: 100%;
    }

    table.footer td.left {
        width: 50%;
        text-align: left;
        font-style: italic;
        vertical-align: top;
    }

    table.footer td.right {
        width: 50%;
        text-align: right;
        vertical-align: top;
    }

    table.footer table.buttons {
        margin-left: auto;
    }

    table.footer table.buttons td {
        text-align: center;
    }

    td.l {
        width: 20%;
        vertical-align: bottom;
    }

    td.m {
        width: 60%;
        text-align: center;
        vertical-align: middle;
    }
    td.r {
        width: 20%;
    }
</style>
</head>

<body class="standard">

<table width="100%">
<tr><td class="l">
&nbsp;
</td>
<td class="m">
<h1>PHP Floating Point Benchmark</h1>
</td>
<td class="r">
&nbsp;
</td>
</tr>
</table>
<hr />

<p>
This document provides an implementation of the floating
point benchmark in PHP.
</p>

<?php

/*      John Walker's Floating Point Benchmark, derived from...

        Marinchip Interactive Lens Design System

                                     John Walker   December 1980

        By John Walker
           http://www.fourmilab.ch/

        This program may be used, distributed, and modified freely as
        long as the origin information is preserved.

        This is a complete optical design raytracing algorithm,
        stripped of its user interface and recast into PHP. It
        not only determines execution speed on an extremely
        floating point (including trig function) intensive
        real-world application, it checks accuracy on an
        algorithm that is exquisitely sensitive to errors.  The
        performance of this program is typically far more
        sensitive to changes in the efficiency of the
        trigonometric library routines than the average floating
        point program.

        This program can be run from the PHP command line interpreter
        as:
            php fbench.php iteration_count >fbench.html
        or from a Web server with:
            http://www.example.com/path/to/fbench.php?iterations=iteration_count
        where iteration_count is the number of iterations to run.
*/

    $ITERATIONS = 100000;       // Default number of iterations

    /* Allow user to specify the number of iterations as the
       first command line argument when run from the command
       line or with a "?iterations=<n>" query when accessed
       via a Web server. */

    if (isset($argv[1])) {
        $ITERATIONS = $argv[1];
    }
    if (isset($_GET['iterations'])) {
        $ITERATIONS = $_GET['iterations'];
    }

    $niter = $ITERATIONS;           // Iteration counter

    $spectral_line = array(         // Spectral lines in which we trace
                                    // (Not all used in this program.)
                            'A'   => 7621.0,
                            'B'   => 6869.955,
                            'C'   => 6562.816,
                            'D'   => 5895.944,
                            'E'   => 5269.557,
                            'F'   => 4861.344,
                            'G\'' => 4340.477,
                            'H'   => 3968.494
                          );

    /*  The test case used in this program is the design for a 4
        inch achromatic telescope objective used as the example
        in Wyld's classic work on ray tracing by hand, given in
        Amateur Telescope Making, Volume 3 (Volume 2 in the 1996
        reprint edition).  */

    $wyld_clear_aperture = 4;
    $wyld_lens = array(
            array(  27.05, 1.5137, 63.6, 0.52  ),
            array( -16.68, 1,       0,   0.138 ),
            array( -16.68, 1.6164, 36.7, 0.38  ),
            array( -78.1,  1,       0,   0     )
    );

    $refarr = array(    /*  Reference results.  These happen to
                            be derived from a run on Microsoft
                            Quick BASIC on the IBM PC/AT.  */
            '   Marginal ray          47.09479120920   0.04178472683',
            '   Paraxial ray          47.08372160249   0.04177864821',
            'Longitudinal spherical aberration:        -0.01106960671',
            '    (Maximum permissible):                 0.05306749907',
            'Offense against sine condition (coma):     0.00008954761',
            '    (Maximum permissible):                 0.00250000000',
            'Axial chromatic aberration:                0.00448229032',
            '    (Maximum permissible):                 0.05306749907'
    );

    /*    Calculate passage through surface

          If the variable $paraxial is true, the trace through the
          surface will be done using the paraxial approximations.
          Otherwise, the normal trigonometric trace will be done.

          This subroutine takes the following global inputs:

          $radius_of_curvature    Radius of curvature of surface
                                  being crossed.  If 0, surface is
                                  plane.

          $object_distance        Distance of object focus from
                                  lens vertex.  If 0, incoming
                                  rays are parallel and
                                  the following must be specified:

          $ray_height             Height of ray from axis.  Only
                                  relevant if $object_distance == 0

          $axis_slope_angle       Angle incoming ray makes with axis
                                  at intercept

          $from_index             Refractive index of medium being left

          $to_index               Refractive index of medium being
                                  entered.

          The outputs are the following global variables:

          $object_distance        Distance from vertex to object focus
                                  after refraction.

          $axis_slope_angle       Angle incoming ray makes with axis
                                  at intercept after refraction.  */

    function transit_surface() {
        global $paraxial, $radius_of_curvature, $object_distance,
               $ray_height, $axis_slope_angle, $from_index,
               $to_index, $object_distance, $axis_slope_angle;

            if ($paraxial) {
               if ($radius_of_curvature != 0) {
                  if ($object_distance == 0) {
                     $axis_slope_angle = 0;
                     $iang_sin = $ray_height / $radius_of_curvature;
                  } else {
                     $iang_sin = (($object_distance -
                        $radius_of_curvature) / $radius_of_curvature) *
                        $axis_slope_angle;
                  }
                  $rang_sin = ($from_index / $to_index) *
                     $iang_sin;
                  $old_axis_slope_angle = $axis_slope_angle;
                  $axis_slope_angle = $axis_slope_angle +
                     $iang_sin - $rang_sin;
                  if ($object_distance != 0) {
                     $ray_height = $object_distance * $old_axis_slope_angle;
                  }
                  $object_distance = $ray_height / $axis_slope_angle;
                  return;
               }
               $object_distance = $object_distance * ($to_index / $from_index);
               $axis_slope_angle = $axis_slope_angle * ($from_index / $to_index);
               return;
            }

            if ($radius_of_curvature != 0) {
               if ($object_distance == 0) {
                  $axis_slope_angle = 0;
                  $iang_sin = $ray_height / $radius_of_curvature;
               } else {
                  $iang_sin = (($object_distance -
                     $radius_of_curvature) / $radius_of_curvature) *
                     sin($axis_slope_angle);
               }
               $iang = asin($iang_sin);
               $rang_sin = ($from_index / $to_index) *
                  $iang_sin;
               $old_axis_slope_angle = $axis_slope_angle;
               $axis_slope_angle = $axis_slope_angle +
                  $iang - asin($rang_sin);
               $sagitta = sin(($old_axis_slope_angle + $iang) / 2);
               $sagitta = 2 * $radius_of_curvature * $sagitta * $sagitta;
               $object_distance = (($radius_of_curvature * sin(
                  $old_axis_slope_angle + $iang)) *
                  (1.0 / tan($axis_slope_angle))) + $sagitta;
               return;
            }

            $rang = -asin(($from_index / $to_index) *
               sin($axis_slope_angle));
            $object_distance = $object_distance * (($to_index *
               cos(-$rang)) / ($from_index *
               cos($axis_slope_angle)));
            $axis_slope_angle = -$rang;
    }

    //  Perform ray trace in specific spectral line

    function trace_line($line, $ray_h) {
        global $paraxial, $radius_of_curvature, $object_distance,
               $ray_height, $axis_slope_angle, $from_index,
               $to_index, $object_distance, $axis_slope_angle;
        global $current_surfaces, $s, $spectral_line;

        $object_distance = 0;
        $ray_height = $ray_h;
        $from_index = 1;

        for ($i = 1; $i <= $current_surfaces; $i++) {
           $radius_of_curvature = $s[$i][1];
           $to_index = $s[$i][2];
           if ($to_index > 1) {
              $to_index = $to_index + (($spectral_line['D'] -
                 $spectral_line[$line]) /
                 ($spectral_line['C'] - $spectral_line['F'])) * (($s[$i][2] - 1) /
                 $s[$i][3]);
           }
           transit_surface();
           $from_index = $to_index;
           if ($i < $current_surfaces) {
              $object_distance = $object_distance - $s[$i][4];
           }
        }
    }

    //  Load test case into working array

    $clear_aperture = $wyld_clear_aperture;
    $current_surfaces = count($wyld_lens);
    for ($i = 0; $i < $current_surfaces; $i++) {
       for ($j = 0; $j < 4; $j++) {
          $s[$i + 1][$j + 1] = $wyld_lens[$i][$j];
       }
    }

    //  Edit the number of iterations with commas
    $tniter = strrev(preg_replace('/(\d\d\d)(?=\d)(?!\d*\.)/', "$1,", strrev($niter)));
    echo "<p>Running $tniter iterations.</p>\n";

    $start_time = microtime(true);      // Start of timed bechmark

    // Perform ray trace the specified number of times.

    for ($itercount = 0; $itercount < $niter; $itercount++) {

       for ($paraxial = 0; $paraxial <= 1; $paraxial++) {

          // Do main trace in D light

          trace_line('D', $clear_aperture / 2);
          $od_sa[$paraxial][0] = $object_distance;
          $od_sa[$paraxial][1] = $axis_slope_angle;
       }
       $paraxial = 0;

       // Trace marginal ray in C

       trace_line('C', $clear_aperture / 2);
       $od_cline = $object_distance;

       // Trace marginal ray in F

       trace_line('F', $clear_aperture / 2);
       $od_fline = $object_distance;

        /*  Compute aberrations of the design

            The longitudinal spherical aberration is just the
            difference between where the D line comes to focus
            for paraxial and marginal rays.  */

       $aberr_lspher = $od_sa[1][0] - $od_sa[0][0];

       /*   The offense against the sine condition is a measure
            of the degree of coma in the design.  We compute it
            as the lateral distance in the focal plane between
            where a paraxial ray and marginal ray in the D line
            come to focus.
             */
       $aberr_osc = 1 - ($od_sa[1][0] * $od_sa[1][1]) /
          (sin($od_sa[0][1]) * $od_sa[0][0]);

       /*   The axial chromatic aberration is the distance between
            where marginal rays in the C and F lines come to focus.  */

       $aberr_lchrom = $od_fline - $od_cline;
       $max_lspher = sin($od_sa[0][1]);

     /*     Compute maximum acceptable values for each aberration

            Maximum longitudinal spherical aberration, which is
            also the maximum for axial chromatic aberration.  This
            is computed for the D line.  */

       $max_lspher = 0.0000926 / ($max_lspher * $max_lspher);
       $max_osc = 0.0025;
       $max_lchrom = $max_lspher;
    }

    $end_time = microtime(true);        // End of timed benchmark

    // Now evaluate the accuracy of the results from the last ray trace

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

    //  Echo the design evaluation

    echo "<div style='margin-left: auto; margin-right: auto; padding: 8px; background-color: #F0F0F0; width: 60%;'>\n";
    echo "<pre>\n";
    for ($i = 0; $i < 8; $i++) {
        echo "$outarr[$i]\n";
    }

    /*  Now compare the edited results with the master values from
        reference executions of this program.  */

    $errors = 0;
    for ($i = 0; $i < 8; $i++) {
        if ($outarr[$i] !== $refarr[$i]) {
            echo "\n";
            printf("Error in results on line %d...\n", $i + 1);
            print("Expected:  $refarr[$i]\n");
            print("Received:  $outarr[$i]\n");
            echo "(Errors)   ";
            $k = strlen($refarr[$i]);
            for ($j = 0; $j < $k; $j++) {
               echo substr($refarr[$i], $j, 1) === substr($outarr[$i], $j, 1) ? ' ' : '^';
               if (substr($refarr[$i], $j, 1) !== substr($outarr[$i], $j, 1)) {
                  $errors++;
              }
            }
            echo "\n";
        }
    }
    if ($errors > 0) {
       printf("\n%d error%s in results.  This is VERY SERIOUS.\n",
          $errors, $errors > 1 ? "s" : "");
    } else {
       printf("\nNo errors in results.\n");
    }
    echo "</pre>\n";
    echo "</div>\n";

    //  Print the benchmark timing
    echo "<p>\nExecution time for $tniter iterations:\n";
    printf("%.3f seconds (%.4f microseconds/iteration).\n</p>\n",
        ($end_time - $start_time),
        (($end_time - $start_time) * 1e6) / $niter);
?>

<hr />

<blockquote class="rights">
<p>
This document is in the public domain.
</p>
</blockquote>

<table class="footer">
<tr>
<td class="left">
by <a href="/">John Walker</a><br />
January, 2018

<p />

<h3 class="nav"><a href="https://www.fourmilab.ch/fbench/fbench.html">Back to <code>fbench</code></a></h3>
<h3 class="nav"><a href="https://www.fourmilab.ch/">Fourmilab Home Page</a></h3>
</td>

<td class="right">
    <table class="buttons">
    <tr><td>
    <form id="feedback" method="post" action="https://www.fourmilab.ch/cgi-bin/FeedbackForm.pl">
    <div>
    <input type="hidden" name="pagetitle" value="PHP Floating Point Benchmark" />
    <input type="hidden" name="backlink"
        value="Back to &lt;cite&gt;PHP Floating Point Benchmark&lt;/cite&gt;" />
    <input type="submit" value=" Send Feedback " />
    </div>
    </form>
    </td></tr>
   </table>
</td>
</tr>
</table>

</body>
</html>
