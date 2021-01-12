    #=  John Walker's Floating Point Benchmark, derived from...

            Marinchip Interactive Lens Design System
                   John Walker   December 1980

                       by John Walker
                   http://www.fourmilab.ch/

        This program may be used, distributed, and modified
        freely as long as the origin information is preserved.

        This is a complete optical design raytracing algorithm,
        stripped of its user interface and recast into Julia. It
        not only determines execution speed on an extremely
        floating point (including trig function) intensive
        real-world application, it checks accuracy on an
        algorithm that is exquisitely sensitive to errors. The
        performance of this program is typically far more
        sensitive to changes in the efficiency of the
        trigonometric library routines than the average floating
        point program.

        Implemented in November 2017 by John Walker.  =#


    #=  Wavelengths of standard spectral lines in Angstroms
              (Not all are used in this program)  =#

    const A_line      = 7621.0;
    const B_line      = 6869.955;
    const C_line      = 6562.816;
    const D_line      = 5895.944;
    const E_line      = 5269.557;
    const F_line      = 4861.344;
    const Gprime_line = 4340.477;
    const H_line      = 3968.494;

    #=  A Surface describes the boundary between two components
        in the Design.  =#

    struct Surface
        curvature_Radius :: Float64
        index_Of_Refraction :: Float64
        dispersion :: Float64
        edge_Thickness :: Float64
    end

    #=  A Design is the specification of the optical assembly
        to be evaluated.  =#

    struct Design
        clearAperture :: Float64
        surf :: Array{Surface}
    end

    #=  The test case used in this program is the design
        for a 4 inch f/12 achromatic telescope objective
        used as the example in Wyld's classic work on ray
        tracing by hand, given in Amateur Telescope Making,
        Volume 3 (Volume 2 in the 1996 reprint edition).  =#

    const WyldLens =
                Design(4.0, #   CurRad  Index   Disp  Edge
                    [ Surface(  27.05, 1.5137, 63.6, 0.52  ),
                      Surface( -16.68, 1.0,     0.0, 0.138 ),
                      Surface( -16.68, 1.6164, 36.7, 0.38  ),
                      Surface( -78.1,  1.0,     0.0, 0.0   )
                     ]);

    @enum AxialIncidence Marginal_Ray Paraxial_Ray

    #   The following global static variables serve as the trace context

    cSurf = 0;
    object_distance = 0.0;
    ray_height = 0.0;
    axis_slope_angle = 0.0;
    from_index = 0.0;

    #=  transitSurface propagates a ray through a Design.

        design                  Design

        axial_incidence         Axial incidence of ray

        line                    Wavelength of ray being traced

        The trace context provides input and output as follows.

        cSurf                   Current surface number

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

        cSurf                   Next surface number

        object_distance         Distance from vertex to object focus
                                after refraction.

        ray_height              Height of ray from axis.

        axis_slope_angle        Angle incoming ray makes with axis
                                at intercept after refraction.

        from_index              The to_index of the input context.
    =#

    function transitSurface(
        axial_incidence :: AxialIncidence,
        d :: Design,
        line :: Float64)

        #   Set context variables from current surface

        radius_of_curvature = d.surf[cSurf].curvature_Radius;
        #   Note that to_index is our local copy
        to_index = d.surf[cSurf].index_Of_Refraction;

        #=  Adjust effective to_index according to the dispersion
            and the wavelength of the line being traced.  =#

        if to_index > 1
            to_index += ((D_line - line) /
                (C_line - F_line)) *
                ((d.surf[cSurf].index_Of_Refraction - 1) /
                d.surf[cSurf].dispersion)
        end

        if axial_incidence == Paraxial_Ray

            #   Paraxial ray

            if radius_of_curvature != 0

                #   Curved surface

                odz = object_distance == 0;
                asaprime = odz ? 0 : axis_slope_angle;
                iangsin = odz ? ray_height / radius_of_curvature :
                                ((object_distance - radius_of_curvature) /
                                 radius_of_curvature) * axis_slope_angle;
                rangsin = (from_index / to_index) * iangsin;
                asadoubleprime = asaprime + iangsin - rangsin;
                rayheightprime = odz ? ray_height :
                                        object_distance * asaprime;
                objectdistanceprime = rayheightprime / asadoubleprime;

                global object_distance = objectdistanceprime -
                    d.surf[cSurf].edge_Thickness;
                global ray_height = objectdistanceprime;
                global axis_slope_angle = asadoubleprime;
                global from_index = to_index;
                global cSurf = cSurf + 1;
            else

                #   Flat surface

                global object_distance = (object_distance *
                    (to_index / from_index)) -
                    d.surf[cSurf].edge_Thickness;
                global axis_slope_angle = axis_slope_angle *
                    (from_index / to_index);
                global from_index = to_index;
                global cSurf = cSurf + 1;
            end

         else

            #   Marginal ray

            if radius_of_curvature != 0

                #   Curved surface

                odz = object_distance == 0;
                asaprime = odz ? 0 : axis_slope_angle;
                iangsin = odz ? ray_height / radius_of_curvature :
                                ((object_distance - radius_of_curvature) /
                                 radius_of_curvature) * sin(axis_slope_angle);
                iang = asin(iangsin);
                rangsin = (from_index / to_index) * iangsin;
                asadoubleprime = asaprime + iang - asin(rangsin);
                sinasaiang = sin((asaprime + iang) / 2);
                sagitta = 2 * radius_of_curvature * sinasaiang * sinasaiang;
                rayheightprime = odz ? ray_height :
                                        object_distance * asaprime;
                objectdistanceprime = ((radius_of_curvature *
                                          sin(asaprime + iang)) *
                                          cot(asadoubleprime)) + sagitta;

                global object_distance = objectdistanceprime -
                    d.surf[cSurf].edge_Thickness;
                global ray_height = rayheightprime;
                global axis_slope_angle = asadoubleprime;
                global from_index = to_index;
                global cSurf = cSurf + 1;
              else

                #   Flat surface

                rang = -(asin((from_index / to_index))) *
                              sin(axis_slope_angle);

                global object_distance = (object_distance *
                    ((to_index * cos(-rang)) / (from_index *
                    cos(axis_slope_angle)))) -
                    d.surf[cSurf].edge_Thickness;
                global axis_slope_angle = -rang;
                global from_index = to_index;
                global cSurf = cSurf + 1;
            end
        end

    end

    #=  traceLine traces a spectral line through a design with a
        specified axial incidence and returns the object
        distance and slope angle where the ray comes to focus.  =#

    function traceLine(
        axial_incidence :: AxialIncidence,
        d :: Design,
        line :: Float64)

        #   Initial trace context

        global cSurf = 1;
        global object_distance = 0;
        global ray_height = d.clearAperture / 2;
        global axis_slope_angle = 0;
        global from_index = 1;

        for i = 1 : length(d.surf)
            transitSurface(axial_incidence, d, line);
        end

        return (object_distance, axis_slope_angle);
    end

    #=  A DesignEvaluation describes the aberrations of a
        Design computed from tracing a variety of rays through
        it.  =#

    struct DesignEvaluation
        dMarginalOD :: Float64;     # D Marginal ray
        dMarginalSA :: Float64;

        dParaxialOD :: Float64;     # D Paraxial ray
        dParaxialSA :: Float64;

        #   Computed aberrations of design
        longitudinalSphericalAberration :: Float64;
        offenseAgainstSineCondition :: Float64;
        axialChromaticAberration :: Float64;

        #   Acceptable maxima for aberrations
        maxLongitudinalSphericalAberration :: Float64;
        maxOffenseAgainstSineCondition :: Float64;
        maxAxialChromaticAberration :: Float64;

    end

    function evaluateDesign(d :: Design) :: DesignEvaluation

        #   D marginal ray
        dMarginalOD, dMarginalSA =
            traceLine(Marginal_Ray, WyldLens, D_line);

        #   D paraxial ray
        dParaxialOD, dParaxialSA =
            traceLine(Paraxial_Ray, WyldLens, D_line);

        #   C marginal ray
        cMarginalOD, _ =
            traceLine(Marginal_Ray, WyldLens, C_line);

        #   F marginal ray
        fMarginalOD, _ =
            traceLine(Marginal_Ray, WyldLens, F_line);

        #   Compute aberrations of the design

        #=  The longitudinal spherical aberration is just the
            difference between where the D line comes to focus
            for paraxial and marginal rays.  =#
        longitudinalSphericalAberration = dParaxialOD - dMarginalOD;

        #=  The offense against the sine condition is a measure
            of the degree of coma in the design.  We compute it
            as the lateral distance in the focal plane between
            where a paraxial ray and marginal ray in the D line
            come to focus.  =#
        offenseAgainstSineCondition = 1 - (dParaxialOD * dParaxialSA) /
            (sin(dMarginalSA) * dMarginalOD);

        #=  The axial chromatic aberration is the distance between
            where marginal rays in the C and F lines come to focus.  =#
        axialChromaticAberration = fMarginalOD - cMarginalOD;

        #   Compute maximum acceptable values for each aberration

        #=  Maximum longitudinal spherical aberration, which is
            also the maximum for axial chromatic aberration.  This
            is computed for the D line.  =#
        sin_dm_sa = sin(dMarginalSA);
        maxLongitudinalSphericalAberration = 0.0000926 / (sin_dm_sa * sin_dm_sa);
        maxOffenseAgainstSineCondition = 0.0025;
        maxAxialChromaticAberration = maxLongitudinalSphericalAberration; #  Same criterion

        return DesignEvaluation(dMarginalOD, dMarginalSA,
                                dParaxialOD, dParaxialSA,
                                longitudinalSphericalAberration,
                                offenseAgainstSineCondition,
                                axialChromaticAberration,
                                maxLongitudinalSphericalAberration,
                                maxOffenseAgainstSineCondition,
                                maxAxialChromaticAberration);
    end

    function evaluationReport(e :: DesignEvaluation) :: Array{String}
        return [
            @sprintf("%15s   %21.11f  %14.11f", "Marginal ray",
                e.dMarginalOD, e.dMarginalSA),
            @sprintf("%15s   %21.11f  %14.11f", "Paraxial ray",
                e.dParaxialOD, e.dParaxialSA),
            @sprintf("Longitudinal spherical aberration:      %16.11f",
                e.longitudinalSphericalAberration),
            @sprintf("    (Maximum permissible):              %16.11f",
                e.maxLongitudinalSphericalAberration),
            @sprintf("Offense against sine condition (coma):  %16.11f",
                e.offenseAgainstSineCondition),
            @sprintf("    (Maximum permissible):              %16.11f",
                e.maxOffenseAgainstSineCondition),
            @sprintf("Axial chromatic aberration:             %16.11f",
                e.axialChromaticAberration),
            @sprintf("    (Maximum permissible):              %16.11f",
                e.maxAxialChromaticAberration)
        ];
    end

    function validateReport(received :: Array{String}) :: UInt

        #=  The test case used in this program is the design
            for a 4 inch f/12 achromatic telescope objective
            used as the example in Wyld's classic work on ray
            tracing by hand, given in Amateur Telescope Making,
            Volume 3 (Volume 2 in the 1996 reprint edition).  =#

        const expected = [
            "   Marginal ray          47.09479120920   0.04178472683",
            "   Paraxial ray          47.08372160249   0.04177864821",
            "Longitudinal spherical aberration:        -0.01106960671",
            "    (Maximum permissible):                 0.05306749907",
            "Offense against sine condition (coma):     0.00008954761",
            "    (Maximum permissible):                 0.00250000000",
            "Axial chromatic aberration:                0.00448229032",
            "    (Maximum permissible):                 0.05306749907"
        ];

        errors = 0;
        if received != expected
            for i = 1 : length(expected)
                if received[i] != expected[i]
                    errors += 1;
                    println("Error in results on line $(i)");
                    println("Expected:  \" $(expected[i])\"");
                    println("Received:  \" $(received[i])\"");
                    print("(Errors)     ");
                    for j = 1 : length(expected[i])
                        print(received[i][j] == expected[i][j] ?
                            " " : "^");
                    end
                    println();
                end
            end
        end

        return errors;
    end

    #=  runBenchmark runs evaluateDesign for the number of
        iterations given by the argument and leaves the result of
        the last run in the global eDesign for analysis and the
        accuracy check.  =#

    function runBenchmark(niter :: UInt)
        for n = 1 : niter
            global eDesign = evaluateDesign(WyldLens);
        end
    end

    iterations = UInt(1000);

    #=  If we're called with a command line argument, use as the
        iteration count.  =#

    if length(ARGS) > 0
        iterations = parse(UInt, ARGS[1]);
    end

    println("Iterations: ", iterations);

    @time runBenchmark(iterations);

    const received = evaluationReport(eDesign);
    const errors = validateReport(received);

    if errors > 0
        println("$(errors) errors in results.  This is VERY SERIOUS.");
    else
       println("No errors in results.");
    end
