/*  John Walker's Floating Point Benchmark, derived from...

             Marinchip Interactive Lens Design System
                    John Walker   December 1980

                        by John Walker
                    http://www.fourmilab.ch/

    This program may be used, distributed, and modified freely
    as long as the origin information is preserved.

    This is a complete optical design raytracing algorithm,
    stripped of its user interface and recast into Scala.  It not
    only determines execution speed on an extremely floating
    point (including trig function) intensive real-world
    application, it checks accuracy on an algorithm that is
    exquisitely sensitive to errors.  The performance of this
    program is typically far more sensitive to changes in the
    efficiency of the trigonometric library routines than the
    average floating point program.

    Implemented in September 2017 by John Walker.  */

    import scala.math

    object fbench {

        // The cot function defined in terms of math.tan

        def cot(x: Double): Double = 1 / math.tan(x)

        //  Wavelengths of standard spectral lines in Angstroms
        //          (Not all are used in this program)

        val (a_line, b_line, c_line, d_line,
             e_line, f_line, gp_line, h_line) =
            (7621.0, 6869.955, 6562.816, 5895.944,
             5269.557, 4861.344, 4340.477, 3968.494)

        /*  The  test case used in this program is the design
            for a 4 inch f/12 achromatic telescope objective
            used as the example in Wyld's classic work on ray
            tracing by hand, given in Amateur Telescope Making,
            Volume 3 (Volume 2 in the 1996 reprint edition).  */

        class Surface_Property(cr: Double,
                               ir: Double,
                               di: Double,
                               et: Double) {
            var curvature_Radius = cr
            var index_Of_Refraction = ir
            var dispersion = di
            var edge_Thickness = et
        }

        val wyldClearAperture = 4.0 // Clear aperture of the following lens design
        val wyldLens = Array(
                        //        CurRad  Index   Disp  Edge
            new Surface_Property(  27.05, 1.5137, 63.6, 0.52  ),
            new Surface_Property( -16.68, 1.0,     0.0, 0.138 ),
            new Surface_Property( -16.68, 1.6164, 36.7, 0.38  ),
            new Surface_Property( -78.1,  1.0,     0.0, 0.0   )
        )

        //  Nomenclature for axial incidence of traced rays

        object Axial_Incidence extends Enumeration {
            type Axial_Incidence = Value
            val Marginal_Ray, Paraxial_Ray = Value
        }

        import Axial_Incidence._

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

        class Trace_Context(ai: Axial_Incidence,
                            rc: Double,
                            od: Double,
                            rh: Double,
                            as: Double,
                            fi: Double,
                            ti: Double) {
            var axial_incidence = ai
            var radius_of_curvature = rc
            var object_distance = od
            var ray_height = rh
            var axis_slope_angle = as
            var from_index = fi
            var to_index = ti
        }

        def transitSurface(tc: Trace_Context): Trace_Context = {
            val tcp = tc

            //  Paraxial ray

            tc.axial_incidence match {

            case Paraxial_Ray =>

                if (tc.radius_of_curvature != 0) {

                    //  Curved surface

                    val odz = tc.object_distance == 0
                    val asaprime = if (odz) { 0 } else { tc.axis_slope_angle }
                    val iangsin = if (odz) { tc.ray_height / tc.radius_of_curvature
                                           } else {
                                    ((tc.object_distance - tc.radius_of_curvature) /
                                     tc.radius_of_curvature) * tc.axis_slope_angle
                                           }
                    val rangsin = (tc.from_index / tc.to_index) * iangsin
                    val asadoubleprime = asaprime + iangsin - rangsin
                    val rayheightprime = if (odz) { tc.ray_height } else {
                                            tc.object_distance * asaprime }
                    val objectdistanceprime = rayheightprime / asadoubleprime

                    tcp.object_distance = objectdistanceprime
                    tcp.ray_height = rayheightprime
                    tcp.axis_slope_angle = asadoubleprime
                } else {

                    //  Flat surface

                    tcp.object_distance = tc.object_distance * (tc.to_index / tc.from_index)
                    tcp.axis_slope_angle = tc.axis_slope_angle * (tc.from_index / tc.to_index)
                }

            case Marginal_Ray =>

                //  Marginal ray

                if (tc.radius_of_curvature != 0) {

                    //  Curved surface

                    val odz = tc.object_distance == 0
                    val asaprime = if (odz) { 0 } else { tc.axis_slope_angle }
                    val iangsin = if (odz) { tc.ray_height / tc.radius_of_curvature
                                           } else {
                                    ((tc.object_distance - tc.radius_of_curvature) /
                                     tc.radius_of_curvature) * math.sin(tc.axis_slope_angle)
                                           }
                    val iang = math.asin(iangsin)
                    val rangsin = (tc.from_index / tc.to_index) * iangsin
                    val asadoubleprime = asaprime + iang - math.asin(rangsin)
                    val sinasaiang = math.sin((asaprime + iang) / 2)
                    val sagitta = 2 * tc.radius_of_curvature * sinasaiang * sinasaiang
                    val rayheightprime = if (odz) { tc.ray_height } else {
                                            tc.object_distance * asaprime }
                    val objectdistanceprime = ((tc.radius_of_curvature *
                                              math.sin(asaprime + iang)) *
                                              cot(asadoubleprime)) + sagitta

                    tcp.object_distance = objectdistanceprime
                    tcp.ray_height = rayheightprime
                    tcp.axis_slope_angle = asadoubleprime
                } else {

                    //  Flat surface

                    val rang = -(math.asin((tc.from_index / tc.to_index))) *
                               math.sin(tc.axis_slope_angle)

                    tcp.object_distance = tc.object_distance * ((tc.to_index *
                                math.cos(-rang)) / (tc.from_index *
                                math.cos(tc.axis_slope_angle)))
                    tcp.axis_slope_angle = -rang
                }
            }

            return tcp
        }

        /*  traceLens propagates a ray of the specified spectral line
            and axial incidence through a design and returns the final
            object distance and axis slope angle, from which we can
            compute the aberrations of the design.  */

        def traceLens(design: Array[Surface_Property], aperture: Double,
                       line: Double, incidence: Axial_Incidence):
                (Double, Double) = {
            var ctx = new Trace_Context(incidence,
                                   0,
                                   0,
                                   aperture / 2,
                                   0,
                                   1,
                                   0)

            for (surf <- design) {
                ctx.radius_of_curvature = surf.curvature_Radius
                ctx.to_index = surf.index_Of_Refraction
                if (ctx.to_index > 1) {
                    ctx.to_index += ((d_line -
                        line) / (c_line - f_line)) *
                        ((surf.index_Of_Refraction - 1) / surf.dispersion)
                }
                ctx = transitSurface(ctx)
                ctx.from_index = ctx.to_index
                ctx.object_distance -= surf.edge_Thickness
            }

            return (ctx.object_distance, ctx.axis_slope_angle)
        }

        /*  The evaluateDesign function performs a ray trace on a given design
            with a specified clear aperture and returns a DesignEvaluation
            which includes the results for the D line and calculation of
            spherical aberration, coma, and chromatic aberration, along with
            the conventional acceptable upper bounds for these quantities.  */

        class DesignEvaluation {
            var dMarginalOD = 0.0       // D Marginal ray
            var dMarginalSA = 0.0

            var dParaxialOD = 0.0       // D Paraxial ray
            var dParaxialSA = 0.0

            //  Computed aberrations of design
            var longitudinalSphericalAberration = 0.0
            var offenseAgainstSineCondition = 0.0
            var axialChromaticAberration = 0.0

            //  Acceptable maxima for aberrations
            var maxLongitudinalSphericalAberration = 0.0
            var maxOffenseAgainstSineCondition = 0.0025
            var maxAxialChromaticAberration = 0.0

            override def toString: String = {
                "Design evaluation:\n" +
                "  D-marginal OD: " + f"$dMarginalOD%16.11f\n" +
                "  D_marginal SA: " + f"$dMarginalSA%16.11f\n" +
                "  D-paraxial OD: " + f"$dParaxialOD%16.11f\n" +
                "  D_paraxial SA: " + f"$dParaxialSA%16.11f\n" +
                "  Long sph aber: " + f"$longitudinalSphericalAberration%16.11f\n" +
                "  Off sine cond: " + f"$offenseAgainstSineCondition%16.11f\n" +
                "  Chromat aberr: " + f"$axialChromaticAberration%16.11f\n" +
                "  Max long aber: " + f"$maxLongitudinalSphericalAberration%16.11f\n" +
                "  Max off sine:  " + f"$maxOffenseAgainstSineCondition%16.11f\n" +
                "  Max chrom abb: " + f"$maxAxialChromaticAberration%16.11f"
            }
        }

        def evaluateDesign(design: Array[Surface_Property], aperture: Double):
                DesignEvaluation = {
            val de = new DesignEvaluation

            //  D marginal ray
            val (od_mD, sa_mD) = traceLens(design, aperture, d_line, Marginal_Ray)
            de.dMarginalOD = od_mD
            de.dMarginalSA = sa_mD

            //  D paraxial ray
            val (od_pD, sa_pD) = traceLens(design, aperture, d_line, Paraxial_Ray)
            de.dParaxialOD = od_pD
            de.dParaxialSA = sa_pD

            //  C marginal ray
            val (cMarginalOD, _) = traceLens(design, aperture,
                c_line, Marginal_Ray)

            //  F marginal ray
            val (fMarginalOD, _) = traceLens(design, aperture,
                f_line, Marginal_Ray)

            //  Compute aberrations of the design

            /*  The longitudinal spherical aberration is just the
                difference between where the D line comes to focus
                for paraxial and marginal rays.  */
            de.longitudinalSphericalAberration = de.dParaxialOD - de.dMarginalOD

            /*  The offense against the sine condition is a measure
                of the degree of coma in the design.  We compute it
                as the lateral distance in the focal plane between
                where a paraxial ray and marginal ray in the D line
                come to focus.  */
            de.offenseAgainstSineCondition = 1 - (de.dParaxialOD * de.dParaxialSA) /
                (math.sin(de.dMarginalSA) * de.dMarginalOD)

            /*  The axial chromatic aberration is the distance between
                where marginal rays in the C and F lines come to focus.  */
            de.axialChromaticAberration = fMarginalOD - cMarginalOD

            //  Compute maximum acceptable values for each aberration

            /*  Maximum longitudinal spherical aberration, which is
                also the maximum for axial chromatic aberration.  This
                is computed for the D line.  */
            val sin_dm_sa = math.sin(de.dMarginalSA)
            de.maxLongitudinalSphericalAberration = 0.0000926 / (sin_dm_sa * sin_dm_sa)
            de.maxAxialChromaticAberration = de.maxLongitudinalSphericalAberration // Same criterion

            return de
        }

        /*  The evaluationReport function takes a
            DesignEvaluation and returns an Array of strings
            containing a primate-readable version of the
            evaluation.  */

        def evaluationReport(de: DesignEvaluation): Array[String] = {
            val mp = "    (Maximum permissible):              "

            Array(
                "   Marginal ray        " +
                    f"${de.dMarginalOD}%16.11f" + "  " +
                    f"${de.dMarginalSA}%14.11f",
                "   Paraxial ray        " +
                    f"${de.dParaxialOD}%16.11f" + "  " +
                    f"${de.dParaxialSA}%14.11f",
                "Longitudinal spherical aberration:      " +
                    f"${de.longitudinalSphericalAberration}%16.11f",
                mp +
                    f"${de.maxLongitudinalSphericalAberration}%16.11f",
                "Offense against sine condition (coma):  " +
                    f"${de.offenseAgainstSineCondition}%16.11f",
                mp +
                    f"${de.maxOffenseAgainstSineCondition}%16.11f",
                "Axial chromatic aberration:             " +
                    f"${de.axialChromaticAberration}%16.11f",
                mp +
                    f"${de.maxAxialChromaticAberration}%16.11f"
            )
        }

        /*  The validateResults function compares a primate-readable report
            from evaluationReport with the archival results from the reference
            implementation (which all language implementations must reproduce
            character-by-character [apart from trivia such as end of line
            conventions and trailing white space]).  It returns a Boolean
            indicating whether the results matched.  */

        def validateResults(er: Array[String]): Boolean = {
            /* Reference results.  These happen to be derived from a run
               on Microsoft Quick BASIC on the IBM PC/AT.  */
            val expectedResults = Array(
                "   Marginal ray          47.09479120920   0.04178472683",
                "   Paraxial ray          47.08372160249   0.04177864821",
                "Longitudinal spherical aberration:        -0.01106960671",
                "    (Maximum permissible):                 0.05306749907",
                "Offense against sine condition (coma):     0.00008954761",
                "    (Maximum permissible):                 0.00250000000",
                "Axial chromatic aberration:                0.00448229032",
                "    (Maximum permissible):                 0.05306749907"
            )

            var j = 0
            var errs = 0
            for (s <- er) {
                if (s != expectedResults(j)) {
                    println("Validation failed on line " + f"${j + 1}%d")
                    println("  Expected: \"(" + expectedResults(j) + ")\"")
                    println("  Received: \"(" + s + ")\"")
                    errs += 1
                }
                j += 1
            }

            errs == 0
        }

        //  Run the benchmark for the specified number of iterations.

        def runBenchmark(iterations: Int) {

            for (n <- 1 to iterations) {
                val de = evaluateDesign(wyldLens, wyldClearAperture)
                if (n == iterations) {
                    val dr = evaluationReport(de)
                    val ok = validateResults(dr)
                    if (!ok) {
                        println("Error(s) detected in results.  This is VERY SERIOUS.")
                    }
                }
            }
        }

        def main(args: Array[String]): Unit = {

            /*  The argument specifies the number for iterations to
                run.  For archival purposes, set the iteration count
                to achieve a run time of around five minutes.  For
                the runs on my machine, I used an iteration count of
                126612394.  */
		
	    val niter = if (args.length > 0) { args(0).toInt } else { 1000 }

            runBenchmark(niter)
        }
    }
