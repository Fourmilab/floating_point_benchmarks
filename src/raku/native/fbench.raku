#! /usr/bin/raku
#`{
        John Walker's Floating Point Benchmark, derived from...

              Marinchip Interactive Lens Design System

                                    John Walker   December 1980

        By John Walker
          http://www.fourmilab.ch/

        This program may be used, distributed, and modified freely as
        long as the origin information is preserved.

        This is a complete optical design raytracing algorithm,
        stripped of its user interface and recast into Raku. It not
        only determines execution speed on an extremely floating point
        (including trig function) intensive real-world application, it
        checks accuracy on an algorithm that is exquisitely sensitive
        to errors.  The performance of this program is typically far
        more sensitive to changes in the efficiency of the
        trigonometric library routines than the average floating point
        program.
}

    use v6;
    use strict;

    #`{ Wavelengths of standard spectral lines in Angstroms
              (Not all are used in this program)  }

    class SpectralLine {
        our $A = 7621.0e0;
        our $B = 6869.955e0;
        our $C = 6562.816e0;
        our $D = 5895.944e0;
        our $E = 5269.557e0;
        our $F = 4861.344e0;
        our $Gprime = 4340.477e0;
        our $H = 3968.494e0;
    }

    #`{ The Surface class represents a lens surface within
        the design.  It has no methods. }

    class Surface {
        has Real $.curvature_Radius;
        has Real $.index_Of_Refraction;
        has Real $.dispersion;
        has Real $.edge_Thickness;
    }

    #`{ The Design class repesents a complete lens design,
        represented as a sequence of Surfaces and a clear
        aperture in inches. }

    class Design {
        has Real $.clearAperture;
        has Surface @.surf;
    }

    #`{ AxialIncidence specifies whether we are tracing a ray
        through the edge of the lens (Marginal_Ray) or one that
        passes through its centre (Paraxial_Ray). }

    enum AxialIncidence <Marginal_Ray Paraxial_Ray>;

    #`{ The TraceContext is where we get down to business and actually
        trace the ray through the lens Design.  It is initialised with
        the Design to be studied, and then its traceLine() method is
        called, specifying the wavelength and axial incidence of each
        line to be traced.  The return is a list of the object distance
        (focal point) and axis slope angle (angle from the axis of the
        ray at the object distance) where it comes to focus on the
        axis. }

    class TraceContext {
        has Design $.des;
        has AxialIncidence $.axial_incidence is rw;
        has Real $.line is rw;

        has Int $!cSurf = 0;
        has Real $!radius_of_curvature;
        has Real $!object_distance;
        has Real $!ray_height;
        has Real $!axis_slope_angle;
        has Real $!from_index;
        has Real $!to_index;

        #`{
        transitSurface propagates a ray across a Surface.

        axial_incidence         Axial incidence of ray

        line                    Wavelength of ray being traced

        radius_of_curvature     Radius of curvature of surface
                                being crossed.  If 0, surface is
                                plane.

        object_distance         Distance of object focus
                                from lens vertex.  If 0,
                                incoming rays are parallel
                                and the following must be
                                specified:

        ray_height              Height of ray from axis.  Only
                                relevant if object_distance == 0

        axis_slope_angle        Angle incoming ray makes with axis
                                at intercept

        from_index              Refractive index of medium being left


        to_index                Refractive index of medium being
                                entered.

        The result of a call on transitSurface updates the
        Trace_Context with the following components modified to
        reflect the geometry of the ray as it exits the surface.

        object_distance         Distance from vertex to object focus
                                after refraction.

        ray_height              Height of ray from axis.

        axis_slope_angle        Angle incoming ray makes with axis
                                at intercept after refraction.
        }

        method !transitSurface() {
            my Real ($iang,         # Incidence angle
                     $rang,         # Refraction angle
                     $iang_sin,     # Incidence angle sin
                     $rang_sin,     # Refraction angle sin
                     $old_axis_slope_angle, $sagitta);

            #   Check for paraxial ray and, if so, use faster
            #   approximations.

            if ($.axial_incidence == Paraxial_Ray) {
                if ($!radius_of_curvature != 0e0) {
                    if ($!object_distance == 0e0) {
                        $!axis_slope_angle = 0e0;
                        $iang_sin = $!ray_height / $!radius_of_curvature;
                    } else {
                        $iang_sin = (($!object_distance -
                        $!radius_of_curvature) / $!radius_of_curvature) *
                            $!axis_slope_angle;
                    }
                    $rang_sin = ($!from_index / $!to_index) * $iang_sin;
                    $old_axis_slope_angle = $!axis_slope_angle;
                    $!axis_slope_angle = $!axis_slope_angle +
                        $iang_sin - $rang_sin;
                    if ($!object_distance != 0e0) {
                        $!ray_height = $!object_distance * $old_axis_slope_angle;
                    }
                    $!object_distance = $!ray_height / $!axis_slope_angle;
                } else {
                    $!object_distance = $!object_distance *
                        ($!to_index / $!from_index);
                    $!axis_slope_angle = $!axis_slope_angle *
                        ($!from_index / $!to_index);
                }
            } else {

                #   Marginal ray.  Do full geometric trace

                if ($!radius_of_curvature != 0e0) {
                    #   Surface is not flat.  Perform general case trace
                    if ($!object_distance == 0e0) {
                        $!axis_slope_angle = 0e0;
                        $iang_sin = $!ray_height / $!radius_of_curvature;
                    } else {
                        $iang_sin = (($!object_distance -
                            $!radius_of_curvature) / $!radius_of_curvature) *
                            sin($!axis_slope_angle);
                    }
                    $iang = asin($iang_sin);
                    $rang_sin = ($!from_index / $!to_index) *
                        $iang_sin;
                    $old_axis_slope_angle = $!axis_slope_angle;
                    $!axis_slope_angle = $!axis_slope_angle +
                        $iang - asin($rang_sin);
                    $sagitta = sin(($old_axis_slope_angle + $iang) / 2e0);
                    $sagitta = 2e0 * $!radius_of_curvature * ($sagitta ** 2);
                    $!object_distance = (($!radius_of_curvature *
                        sin($old_axis_slope_angle + $iang)) *
                        (1e0 / tan($!axis_slope_angle))) + $sagitta;
                } else {

                    #   Surface is flat.  Use simpler calculation

                    $rang = -asin(($!from_index / $!to_index) *
                        sin($!axis_slope_angle));
                    $!object_distance = $!object_distance * (($!to_index *
                        cos(-$rang)) / ($!from_index *
                        cos($!axis_slope_angle)));
                    $!axis_slope_angle = -$rang;
                }
            }
        }

        # Trace a line through the entire design

        method traceLine(Real :$line, AxialIncidence :$incidence) {
            $.axial_incidence = $incidence;
            $.line = $line;
            $!cSurf = 0;
            $!radius_of_curvature = $!object_distance =
                $!axis_slope_angle = 0e0;
            $!to_index = 0e0;
            $!ray_height = $!des.clearAperture / 2;
            $!from_index = 1e0;

            repeat {
                $!radius_of_curvature = $!des.surf[$!cSurf].curvature_Radius;
                $!to_index = $!des.surf[$!cSurf].index_Of_Refraction;
                if ($!to_index > 1) {
                    $!to_index += (($SpectralLine::D - $!line) /
                        ($SpectralLine::C - $SpectralLine::F)) *
                        (($!des.surf[$!cSurf].index_Of_Refraction - 1e0) /
                         $!des.surf[$!cSurf].dispersion);
                }
                self!transitSurface();
                $!from_index = $!to_index;
                $!object_distance -= $!des.surf[$!cSurf].edge_Thickness;
                $!cSurf++;
            } while ($!cSurf < $!des.surf.elems);

            return ($!object_distance, $!axis_slope_angle);
        }
    }


    #`{ A DesignEvaluation provides tools to analyse designs.  It takes
        a design, traces rays through it in various wavelengths and
        axial incidences, and computes its aberrations compared to
        acceptable standards. }

    class DesignEvaluation {
        has Design $.des;

        has Real $!cMarginalOD;      # C marginal ray
        has Real $!fMarginalOD;      # F marginal ray

        has Real $!dMarginalOD;      # D Marginal ray
        has Real $!dMarginalSA;

        has Real $!dParaxialOD;      # D Paraxial ray
        has Real $!dParaxialSA;

        #   Computed aberrations of design
        has Real $!longitudinalSphericalAberration;
        has Real $!offenseAgainstSineCondition;
        has Real $!axialChromaticAberration;

        #   Acceptable maxima for aberrations
        has Real $!maxLongitudinalSphericalAberration;
        has Real $!maxOffenseAgainstSineCondition = 0.0025e0;
        has Real $!maxAxialChromaticAberration;

        has Str @.received;

        method evaluate() {

            my TraceContext $tc = TraceContext.new(des => $.des);

            #   D marginal ray
            ($!dMarginalOD, $!dMarginalSA) = $tc.traceLine(
                line => $SpectralLine::D, incidence => Marginal_Ray);

            #   D paraxial ray
            ($!dParaxialOD, $!dParaxialSA) = $tc.traceLine(
                line => $SpectralLine::D, incidence => Paraxial_Ray);

            #   C marginal ray
            ($!cMarginalOD, $_) = $tc.traceLine(
                line => $SpectralLine::C, incidence => Marginal_Ray);

            #   F marginal ray
            ($!fMarginalOD, $_) = $tc.traceLine(
                line => $SpectralLine::F, incidence => Marginal_Ray);

            #   Compute aberrations of the design

            #`{ The longitudinal spherical aberration is just the
                difference between where the D line comes to focus
                for paraxial and marginal rays. }
            $!longitudinalSphericalAberration = $!dParaxialOD - $!dMarginalOD;

            #`{ The offense against the sine condition is a measure
                of the degree of coma in the design.  We compute it
                as the lateral distance in the focal plane between
                where a paraxial ray and marginal ray in the D line
                come to focus. }
            $!offenseAgainstSineCondition = 1 - ($!dParaxialOD * $!dParaxialSA) /
                (sin($!dMarginalSA) * $!dMarginalOD);

            #`{ The axial chromatic aberration is the distance between
                where marginal rays in the C and F lines come to focus. }
            $!axialChromaticAberration = $!fMarginalOD - $!cMarginalOD;

            #   Compute maximum acceptable values for each aberration

            #`{ Maximum longitudinal spherical aberration, which is
                also the maximum for axial chromatic aberration.  This
                is computed for the D line. }
            $!maxLongitudinalSphericalAberration = 0.0000926 /
                (sin($!dMarginalSA) ** 2);
            # Same criteria for longitudinal spherical and chromatic aberration
            $!maxAxialChromaticAberration = $!maxLongitudinalSphericalAberration;
        }

        method report() {
            my Str $mp = "    (Maximum permissible):              %16.11f";
            my Str $ry = "%15s   %21.11f  %14.11f";

            @.received = ( );
            push(@.received,
                 sprintf($ry, "Marginal ray", $!dMarginalOD, $!dMarginalSA));
            push(@.received,
                 sprintf($ry, "Paraxial ray", $!dParaxialOD, $!dParaxialSA));
            push(@.received,
                 sprintf("Longitudinal spherical aberration:      %16.11f",
                    $!longitudinalSphericalAberration));
            push(@.received, sprintf($mp, $!maxLongitudinalSphericalAberration));
            push(@.received,
                 sprintf("Offense against sine condition (coma):  %16.11f",
                    $!offenseAgainstSineCondition));
            push(@.received, sprintf($mp, $!maxOffenseAgainstSineCondition));
            push(@.received,
                 sprintf("Axial chromatic aberration:             %16.11f",
                    $!axialChromaticAberration));
            push(@.received, sprintf($mp, $!maxAxialChromaticAberration));
        }

        method validate() {

            constant @expected = (  # Reference results.  These happen to
                                    # be derived from a run on Microsoft
                                    # Quick BASIC on the IBM PC/AT.

                "   Marginal ray          47.09479120920   0.04178472683",
                "   Paraxial ray          47.08372160249   0.04177864821",
                "Longitudinal spherical aberration:        -0.01106960671",
                "    (Maximum permissible):                 0.05306749907",
                "Offense against sine condition (coma):     0.00008954761",
                "    (Maximum permissible):                 0.00250000000",
                "Axial chromatic aberration:                0.00448229032",
                "    (Maximum permissible):                 0.05306749907"
            );
            my Int $errors = 0;
            loop (my Int $n = 0; $n < 8; $n++) {
               if (@.received[$n] ne @expected[$n]) {
                  printf("\nError in results on line %d...\n", $n + 1);
                  print("Expected:  @expected[$n]\n");
                  print("Received:  @.received[$n]\n");
                  print("(Errors)   ");
                  my Int $k = @expected[$n].chars;
                  loop (my Int $j = 0; $j < $k; $j++) {
                     print(substr(@expected[$n], $j, 1) eq substr(@.received[$n], $j, 1) ?? ' ' !! '^');
                     if (substr(@expected[$n], $j, 1) ne substr(@.received[$n], $j, 1)) {
                        $errors++;
                    }
                  }
                  print("\n");
               }
            }
            return $errors;
        }
    }

    #`{
        The test case used in this program is the design
        for a 4 inch f/12 achromatic telescope objective
        used as the example in Wyld's classic work on ray
        tracing by hand, given in Amateur Telescope Making,
        Volume 3 (Volume 2 in the 1996 reprint edition).
    }

    constant $WyldLens = Design.new(
                            clearAperture => 4.0e0,
                            surf => (
                                Surface.new(
                                    curvature_Radius => 27.05e0,
                                    index_Of_Refraction => 1.5137e0,
                                    dispersion => 63.6e0,
                                    edge_Thickness => 0.52e0
                                ),
                                Surface.new(
                                    curvature_Radius => -16.68e0,
                                    index_Of_Refraction => 1.0e0,
                                    dispersion => 0.0e0,
                                    edge_Thickness => 0.138e0
                                ),
                                Surface.new(
                                    curvature_Radius => -16.68e0,
                                    index_Of_Refraction => 1.6164e0,
                                    dispersion => 36.7e0,
                                    edge_Thickness => 0.38e0
                                ),
                                Surface.new(
                                    curvature_Radius => -78.1e0,
                                    index_Of_Refraction => 1.0e0,
                                    dispersion => 0.0e0,
                                    edge_Thickness => 0.0e0
                                )
                            )
             );


    # Process the number of iterations argument, if one is supplied.
    my Int $niter = 1000;
    my Bool $batch = False;

    if (@*ARGS.end >= 0) {
        my Str $siter = shift(@*ARGS);
        if ($siter !~~ m/^\-?\d+$/) {
            print q:to/EOD/;
                This is John Walker's floating point accuracy and
                performance benchmark program.  You call it with

                    raku fbench.pl <itercount>

                where <itercount> is the number of iterations to be
                executed.  Archival timings should be made with the
                iteration count set so that roughly five minutes of
                execution is timed.  For non-interactive mode, use a
                negative <itercount>.
                EOD
            exit(0);
        }
        $niter = Int($siter);
        if ($niter < 0) {
            $niter = abs($niter);
            $batch = True;
        }
    }

    my DesignEvaluation $de = DesignEvaluation.new(des => $WyldLens);
    if (!$batch) {
        my Real $nik = $niter / 1000;
        print qq:to/EOD/;
            Ready to begin John Walker's floating point accuracy and
            performance benchmark.  $niter iterations will be made.

            Measured run time in seconds should be divided by $nik to
            normalise for reporting results.  For archival results,
            adjust iteration count so the benchmark runs about five
            minutes.

            EOD
        prompt("Press return to begin benchmark: ");
    }
    my Instant $stime = now;
    for (1 .. $niter) {
        $de.evaluate();
    }
    my Instant $etime = now;
    prompt("Stop the timer: ") unless $batch;

    say();
    $de.report();
    say(join("\n", $de.received()));
    my Int $errs = $de.validate();
    if ($errs == 0) {
        printf("\nNo errors in results.\n");
    } else {
        printf("\n%d error%s in results.  This is VERY SERIOUS.\n",
          $errs, $errs > 1 ?? "s" !! "");
    }

    my Duration $rtime = ($etime - $stime);
    printf("Run time for %d iterations: %.7f seconds.\n", $niter, $rtime);
    printf("Run time per 1000 iterations: %.7f seconds.\n", $rtime / ($niter / 1000));
    printf("Microseconds per iteration: %.4f\n", ($rtime * 1_000_000) / $niter);
