/*  John Walker's Floating Point Benchmark, derived from...

             Marinchip Interactive Lens Design System
                    John Walker   December 1980

                        by John Walker
                    http://www.fourmilab.ch/

    This program may be used, distributed, and modified freely
    as long as the origin information is preserved.

    This is a complete optical design raytracing algorithm,
    stripped of its user interface and recast into Chapel.  It
    not only determines execution speed on an extremely floating
    point (including trig function) intensive real-world
    application, it checks accuracy on an algorithm that is
    exquisitely sensitive to errors.  The performance of this
    program is typically far more sensitive to changes in the
    efficiency of the trigonometric library routines than the
    average floating point program.

    Implemented in October 2017 by John Walker.  */

module fbench {

//  use IO.FormattedIO;     // As of chpl 1.16.0, used automatically

    //  You can set iterations with ./fbench --iterations=n
    config const iterations: uint = 1000000;

    //  Compile with "--set partrace=true" to enable parallel trace
    config param partrace: bool = false;    // Parallelise trace ?

    /*  Compile with "--set pariter=n" to enable parallel iteration
        of the main iterations with n tasks.  The default of 1
        suppresses parallel execution of iterations.  */
    config param pariter: uint = 1;         // Paralleise iterations

    /*  Wavelengths of standard spectral lines in Angstroms
                (Not all are used in this program)  */

    type Wavelength = real;

    param a_line:  Wavelength = 7621.0,
          b_line:  Wavelength = 6869.955,
          c_line:  Wavelength = 6562.816,
          d_line:  Wavelength = 5895.944,
          e_line:  Wavelength = 5269.557,
          f_line:  Wavelength = 4861.344,
          gp_line: Wavelength = 4340.477,
          h_line:  Wavelength = 3968.494;

    record Surface_Property {
        var curvature_Radius,
              index_Of_Refraction,
              dispersion,
              edge_Thickness : real;
    }

    /*  The test case used in this program is the design for a
        4 inch f/12 achromatic telescope objective used as the
        example in Wyld's classic work on ray tracing by hand,
        given in Amateur Telescope Making, Volume 3 (Volume 2 in
        the 1996 reprint edition).  */

    param wyldClearAperture = 4.0;

    const wyldLens = [ //     CurRad  Index   Disp  Edge
        new Surface_Property(  27.05, 1.5137, 63.6, 0.52  ),
        new Surface_Property( -16.68, 1.0,     0.0, 0.138 ),
        new Surface_Property( -16.68, 1.6164, 36.7, 0.38  ),
        new Surface_Property( -78.1,  1.0,     0.0, 0.0   )
    ];

    /*  The Axial_Incidence type specifies whether we're tracing a
        marginal or paraxial ray in transit_surface.  */

    enum Axial_Incidence {
        Marginal_Ray,
        Paraxial_Ray
    }

    //  The cot function is defined in terms of tan

    proc cot(x: real) {
        return 1 / tan(x);
    }

    /*
                       Trace_Context

        axial_incidence         Axial incidence of ray

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

        The result of a call on transitSurface is a
        Trace_Context with the following components modified
        to reflect the geometry of the ray as it exits the
        surface.

        object_distance         Distance from vertex to object focus
                                after refraction.

        ray_height              Height of ray from axis.

        axis_slope_angle        Angle incoming ray makes with axis
                                at intercept after refraction.
    */

    record Trace_Context {
        var axial_incidence: Axial_Incidence;
        var radius_of_curvature,
            object_distance,
            ray_height,
            axis_slope_angle,
            from_index,
            to_index: real;
    }

    proc transitSurface(ref tc: Trace_Context) {
        select tc.axial_incidence {

            when Axial_Incidence.Paraxial_Ray {

                //  Paraxial ray

                if (tc.radius_of_curvature != 0) {

                    //  Curved surface

                    const odz = tc.object_distance == 0;
                    var asaprime, iangsin, rangsin, asadoubleprime,
                        rayheightprime, objectdistanceprime: real;
                    if (odz) {
                        asaprime = 0;
                        iangsin = tc.ray_height / tc.radius_of_curvature;
                    } else {
                        asaprime = tc.axis_slope_angle;
                        iangsin = ((tc.object_distance - tc.radius_of_curvature) /
                                     tc.radius_of_curvature) * tc.axis_slope_angle;
                    }
                    rangsin = (tc.from_index / tc.to_index) * iangsin;
                    asadoubleprime = asaprime + iangsin - rangsin;
                    if (odz) {
                        rayheightprime = tc.ray_height;
                    } else {
                        rayheightprime = tc.object_distance * asaprime;
                    }
                    objectdistanceprime = rayheightprime / asadoubleprime;

                    tc.object_distance = objectdistanceprime;
                    tc.ray_height = rayheightprime;
                    tc.axis_slope_angle = asadoubleprime;
                } else {

                    //  Flat surface

                    tc.object_distance = tc.object_distance * (tc.to_index / tc.from_index);
                    tc.axis_slope_angle = tc.axis_slope_angle * (tc.from_index / tc.to_index);
                }
            }

            when Axial_Incidence.Marginal_Ray {

                //  Marginal ray

                if (tc.radius_of_curvature != 0) {

                    //  Curved surface

                    const odz = tc.object_distance == 0;
                    var asaprime, iangsin, iang, rangsin, asadoubleprime,
                        sinasaiang, sagitta, rayheightprime, objectdistanceprime: real;
                    if (odz) {
                        asaprime = 0;
                        iangsin = tc.ray_height / tc.radius_of_curvature;
                    } else {
                        asaprime = tc.axis_slope_angle;
                        iangsin = ((tc.object_distance - tc.radius_of_curvature) /
                                     tc.radius_of_curvature) * sin(tc.axis_slope_angle);
                    }
                    iang = asin(iangsin);
                    rangsin = (tc.from_index / tc.to_index) * iangsin;
                    asadoubleprime = asaprime + iang - asin(rangsin);
                    sinasaiang = sin((asaprime + iang) / 2);
                    sagitta = 2 * tc.radius_of_curvature * sinasaiang * sinasaiang;
                    if (odz) {
                        rayheightprime = tc.ray_height;
                    } else {
                        rayheightprime = tc.object_distance * asaprime;
                    }
                    objectdistanceprime = ((tc.radius_of_curvature *
                                            sin(asaprime + iang)) *
                                            cot(asadoubleprime)) + sagitta;

                    tc.object_distance = objectdistanceprime;
                    tc.ray_height = rayheightprime;
                    tc.axis_slope_angle = asadoubleprime;
                } else {

                    //  Flat surface

                    const rang = -(asin((tc.from_index / tc.to_index))) *
                                   sin(tc.axis_slope_angle);
                    tc.object_distance = tc.object_distance * ((tc.to_index *
                                cos(-rang)) / (tc.from_index *
                                cos(tc.axis_slope_angle)));
                    tc.axis_slope_angle = -rang;
                }
            }
        }
    }


    /*  traceLens propagates a ray of the specified spectral
        line and axial incidence through a design and returns a
        tuple of the final object distance and axis slope angle,
        from which we can compute the aberrations of the
        design.  */

    proc traceLens(design: [], aperture: real,
                   line: Wavelength, incidence: Axial_Incidence) : (real, real) {
        var ctx = new Trace_Context(incidence,
                               0,
                               0,
                               aperture / 2,
                               0,
                               1,
                               0);

        for surf in design {
            ctx.radius_of_curvature = surf.curvature_Radius;
            ctx.to_index = surf.index_Of_Refraction;
            if (ctx.to_index > 1) {
                ctx.to_index += ((d_line -
                    line) / (c_line - f_line)) *
                    ((surf.index_Of_Refraction - 1) / surf.dispersion);
            }
            transitSurface(ctx);
            ctx.from_index = ctx.to_index;
            ctx.object_distance -= surf.edge_Thickness;
        }

        return (ctx.object_distance, ctx.axis_slope_angle);
    }

    /*  The evaluateDesign function performs a ray trace on a given design
        with a specified clear aperture and returns a DesignEvaluation
        which includes the results for the D line and calculation of
        spherical aberration, coma, and chromatic aberration, along with
        the conventional acceptable upper bounds for these quantities.  */

    record DesignEvaluation {
        var dMarginalOD,            // D Marginal ray
            dMarginalSA,

            dParaxialOD,            // D Paraxial ray
            dParaxialSA,

        //  Computed aberrations of design
            longitudinalSphericalAberration,
            offenseAgainstSineCondition,
            axialChromaticAberration,

        //  Acceptable maxima for aberrations
            maxLongitudinalSphericalAberration,
            maxOffenseAgainstSineCondition = 0.0025,
            maxAxialChromaticAberration : real;
    }

    proc evaluateDesign(design: [], aperture: real) : DesignEvaluation {
        var de: DesignEvaluation;
        var cMarginalOD, fMarginalOD: real;

        if (!partrace) {
           //  D marginal ray
           (de.dMarginalOD, de.dMarginalSA) = traceLens(design, aperture, d_line,
               Axial_Incidence.Marginal_Ray);

           //  D paraxial ray
           (de.dParaxialOD, de.dParaxialSA) = traceLens(design, aperture, d_line,
               Axial_Incidence.Paraxial_Ray);

           //  C marginal ray
           (cMarginalOD, _) = traceLens(design, aperture,
               c_line, Axial_Incidence.Marginal_Ray);

           //  F marginal ray
           (fMarginalOD, _) = traceLens(design, aperture,
               f_line, Axial_Incidence.Marginal_Ray);

        } else {
            var dmOD, dmSA, dpOD, dpSA: real;

            /*  Let me say a few words about the parallel computation
                done below.  We perform the ray trace of the four
                spectral lines as separate tasks, specifying a task
                intent of "ref" so that each task can communicate
                its results (object distance and axis slope angle)
                to the enclosing context when it's done.  This
                doesn't require any synchronisation since each task
                sets a variable or variable which isn't referenced
                by any other task or the parent until all of the tasks
                are done.  We could actually store directly into
                fields of the DesignEvaluation if we passed it by
                reference, but people might be scared by that, so
                we use variables private to this context.  */

            cobegin with (ref dmOD, ref dmSA, ref dpOD, ref dpSA,
                    ref cMarginalOD, ref fMarginalOD) {
                {
                    //  D marginal ray
                    (dmOD, dmSA) = traceLens(design, aperture, d_line,
                        Axial_Incidence.Marginal_Ray);
                }

                {
                    //  D paraxial ray
                    (dpOD, dpSA) = traceLens(design, aperture, d_line,
                        Axial_Incidence.Paraxial_Ray);
                }

                {
                    //  C marginal ray
                    (cMarginalOD, _) = traceLens(design, aperture,
                        c_line, Axial_Incidence.Marginal_Ray);
                }

                {
                    //  F marginal ray
                    (fMarginalOD, _) = traceLens(design, aperture,
                        f_line, Axial_Incidence.Marginal_Ray);
                }
            }
            (de.dMarginalOD, de.dMarginalSA) = (dmOD, dmSA);
            (de.dParaxialOD, de.dParaxialSA) = (dpOD, dpSA);
        }

        //  Compute aberrations of the design

        /*  The longitudinal spherical aberration is just the
            difference between where the D line comes to focus
            for paraxial and marginal rays.  */
        de.longitudinalSphericalAberration = de.dParaxialOD - de.dMarginalOD;

        /*  The offense against the sine condition is a measure
            of the degree of coma in the design.  We compute it
            as the lateral distance in the focal plane between
            where a paraxial ray and marginal ray in the D line
            come to focus.  */
        de.offenseAgainstSineCondition = 1 - (de.dParaxialOD * de.dParaxialSA) /
            (sin(de.dMarginalSA) * de.dMarginalOD);

        /*  The axial chromatic aberration is the distance between
            where marginal rays in the C and F lines come to focus.  */
        de.axialChromaticAberration = fMarginalOD - cMarginalOD;

        //  Compute maximum acceptable values for each aberration

        /*  Maximum longitudinal spherical aberration, which is
            also the maximum for axial chromatic aberration.  This
            is computed for the D line.  */
        var sin_dm_sa = sin(de.dMarginalSA);
        de.maxLongitudinalSphericalAberration = 0.0000926 / (sin_dm_sa * sin_dm_sa);
        de.maxAxialChromaticAberration =
            de.maxLongitudinalSphericalAberration; // Same criterion

        return de;
    }

    /*  The evaluationReport function takes a DesignEvaluation
        and returns an Array of strings containing a
        primate-readable version of the evaluation.  */

    proc evaluationReport(de: DesignEvaluation) : [1..8] string {
        const mp = "    (Maximum permissible):              %16.11dr";
        var report: [1..8] string;

        try! {
            report = [
                "   Marginal ray        %16.11dr  %14.11dr".format(
                    de.dMarginalOD, de.dMarginalSA),
                "   Paraxial ray        %16.11dr  %14.11dr".format(
                    de.dParaxialOD, de.dParaxialSA),
                "Longitudinal spherical aberration:      %16.11dr".format(
                    de.longitudinalSphericalAberration),
                mp.format(de.maxLongitudinalSphericalAberration),
                "Offense against sine condition (coma):  %16.11dr".format(
                    de.offenseAgainstSineCondition),
                mp.format(de.maxOffenseAgainstSineCondition),
                "Axial chromatic aberration:             %16.11dr".format(
                    de.axialChromaticAberration),
                mp.format(de.maxAxialChromaticAberration)
            ];
        }

        return report;
    }

    /*  The validateResults function compares a primate-readable report
        from evaluationReport with the archival results from the reference
        implementation (which all language implementations must reproduce
        character-by-character [apart from trivia such as end of line
        conventions and trailing white space]).  It returns a Boolean
        indicating whether the results matched.  */

    proc validateResults(er: [1..8] string) : bool {
        /* Reference results.  These happen to be derived from a run
           on Microsoft Quick BASIC on the IBM PC/AT.  */
        const expectedResults = [
            "   Marginal ray          47.09479120920   0.04178472683",
            "   Paraxial ray          47.08372160249   0.04177864821",
            "Longitudinal spherical aberration:        -0.01106960671",
            "    (Maximum permissible):                 0.05306749907",
            "Offense against sine condition (coma):     0.00008954761",
            "    (Maximum permissible):                 0.00250000000",
            "Axial chromatic aberration:                0.00448229032",
            "    (Maximum permissible):                 0.05306749907"
        ];

        var j = 1;
        var errs = 0;
        for s in er {
            if (s != expectedResults[j]) {
                try! {
                writeln('Validation failed on line %u'.format(j));
                writeln('  Expected: "%s"'.format(expectedResults[j]));
                writeln('  Received: "%s"'.format(s));
                }
                errs += 1;
            }
            j += 1;
        }

        return errs == 0;
    }

    //  Run the benchmark for the specified number of iterations.

    proc runBenchmark(iterations: uint) {

        if pariter == 1 {
            for n in 1..iterations {
                const de = evaluateDesign(wyldLens, wyldClearAperture);
                if n == iterations {
                    const dr = evaluationReport(de);
                    if !validateResults(dr) {
                        writeln("Error(s) detected in results.  This is VERY SERIOUS.");
                    }
                }
            }
        } else {
            const extras = iterations % pariter;
            coforall i in 1..pariter {
                for n in 1..((iterations / pariter) +
                        (if i == 0 then extras else 0)) {
                    const de = evaluateDesign(wyldLens, wyldClearAperture);
                    if i == pariter && n == iterations {
                        const dr = evaluationReport(de);
                        if !validateResults(dr) {
                            writeln("Error(s) detected in results.  This is VERY SERIOUS.");
                        }
                    }
                }
            }
        }
    }


    runBenchmark(iterations);
}
