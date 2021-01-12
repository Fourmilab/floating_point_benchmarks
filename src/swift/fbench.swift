/*  John Walker's Floating Point Benchmark, derived from...

        Marinchip Interactive Lens Design System
               John Walker   December 1980

                   by John Walker
               http://www.fourmilab.ch/

    This program may be used, distributed, and modified freely as
    long as the origin information is preserved.

    This is a complete optical design raytracing algorithm,
    stripped of its user interface and recast into Swift.
    It not only determines execution speed on an extremely
    floating point (including trig function) intensive
    real-world application, it checks accuracy on an algorithm
    that is exquisitely sensitive to errors. The performance of
    this program is typically far more sensitive to changes in
    the efficiency of the trigonometric library routines than the
    average floating point program.

    Implemented in December 2016 by John Walker.  */

    import Foundation

    /* Wavelengths of standard spectral lines in Angstroms
              (Not all are used in this program) */

    enum spectralLine: Double {
    	case A = 7621.0
	case B = 6869.955
	case C = 6562.816
	case D = 5895.944
	case E = 5269.557
	case F = 4861.344
	case Gprime = 4340.477
	case H = 3968.494
    }

    /* The  test case used in this program is the design for a 4 inch
       f/12 achromatic telescope objective used as the example in Wyld's
       classic work on ray tracing by hand, given in Amateur Telescope
       Making, Volume 3 (Volume 2 in the 1996 reprint edition). */

    struct SurfaceProperty {
    	var curvatureRadius: Double
	var indexOfRefraction: Double
	var dispersion: Double
	var edgeThickness: Double
    }

    let wyldClearAperture = 4.0
    let wyldLens: [SurfaceProperty] = [
    	    	      SurfaceProperty(	    	    	// Crown glass front element
			curvatureRadius:    27.05,
			indexOfRefraction:   1.5137,
			dispersion: 	    63.6,
			edgeThickness:	     0.52
		      ),
    	    	      SurfaceProperty(
			curvatureRadius:   -16.68,
			indexOfRefraction:   1.0,
			dispersion: 	     0.0,
			edgeThickness:	     0.138
		      ),
    	    	      SurfaceProperty(	    	    	// Flint glass rear element
			curvatureRadius:   -16.68,
			indexOfRefraction:   1.6164,
			dispersion: 	    36.7,
			edgeThickness:	     0.38
		      ),
    	    	      SurfaceProperty(
			curvatureRadius:   -78.1,
			indexOfRefraction:   1.0,
			dispersion: 	     0.0,
			edgeThickness:	     0.0
		      )
    	       ]

    /*
    	transitSurface is the heart of the ray tracing algorithm.  It
	takes a TraceContext as input, computes the refraction of the
	ray as the result of traversing the surface, and returns a
	TraceContext containing the state of the ray after it emerges
	from the transition.  The fields in the TraceContext are as
	follows.

    	axialIncidence	    	Is ray marginal or paraxial?

        radiusOfCurvature       Radius of curvature of surface
                                being crossed.  If 0, surface is
                                plane.

        objectDistance          Distance of object focus from
                                lens vertex.  If 0, incoming
                                rays are parallel and
                                the following must be specified:

        rayHeight               Height of ray from axis.  Only
                                relevant if objectDistance == 0

        axisSlopeAngle          Angle incoming ray makes with axis
                                at intercept

        fromIndex               Refractive index of medium being left


        toIndex                 Refractive index of medium being
                                entered

        The result of a call on transitSurface is a TraceContext
        with the following components modified to reflect the geometry
        of the ray as it exits the surface.

        objectDistance          Distance from vertex to object focus
                                after refraction

        rayHeight               Height of ray from axis

        axisSlopeAngle          Angle incoming ray makes with axis
                                at intercept after refraction
    */

    enum AxialIncidence {   	// Specify whether the incident ray is
    	case marginalRay    	// marginal or paraxial
	case paraxialRay
    }

    struct TraceContext {
    	var axialIncidence: AxialIncidence
	var radiusOfCurvature: Double
	var objectDistance: Double
	var rayHeight: Double
	var axisSlopeAngle: Double
	var fromIndex: Double
	var toIndex: Double
    }

    func transitSurface(_ tc: TraceContext) -> TraceContext {
    	var tcp = tc

	/*  There are four cases: the ray can be either paraxial or
	    marginal and, for each, either traversing a flat or a
	    curved surface.  */

	if (tc.axialIncidence == AxialIncidence.paraxialRay) {

	    //  Paraxial ray

	    if tc.radiusOfCurvature != 0 {

	    	//  Curved surface

		let odz = tc.objectDistance == 0
		let asaprime = odz ? 0 : tc.axisSlopeAngle
		let iangsin = odz ? tc.rayHeight / tc.radiusOfCurvature :
		    	    	    ((tc.objectDistance - tc.radiusOfCurvature) /
				     tc.radiusOfCurvature) * tc.axisSlopeAngle
		let rangsin = (tc.fromIndex / tc.toIndex) * iangsin
		let asadoubleprime = asaprime + iangsin - rangsin
		let rayheightprime = odz ? tc.rayHeight : tc.objectDistance * asaprime
		let objectdistanceprime = rayheightprime / asadoubleprime

		tcp.objectDistance = objectdistanceprime
		tcp.rayHeight = rayheightprime
		tcp.axisSlopeAngle = asadoubleprime

	    } else {

	    	//  Flat surface

		tcp.objectDistance = tc.objectDistance * (tc.toIndex / tc.fromIndex)
		tcp.axisSlopeAngle = tc.axisSlopeAngle * (tc.fromIndex / tc.toIndex)
	    }

	} else { // tc.axialIncidence == AxialIncidence.marginalRay

	    //	Marginal ray
	    if tc.radiusOfCurvature != 0 {
	    	func cot(_ x: Double) -> Double {
		    return 1 / tan(x)
		}

	    	//  Curved surface

		let odz = tc.objectDistance == 0
		let asaprime = odz ? 0 : tc.axisSlopeAngle
    	    	let iangsin = odz ? tc.rayHeight / tc.radiusOfCurvature :
		    	((tc.objectDistance - tc.radiusOfCurvature) /
			tc.radiusOfCurvature) * sin(tc.axisSlopeAngle)
		let iang = asin(iangsin)
		let rangsin = (tc.fromIndex / tc.toIndex) * iangsin
		let asadoubleprime = asaprime + iang - asin(rangsin)
		let sinasaiang = sin((asaprime + iang) / 2)
            	let sagitta = 2 * tc.radiusOfCurvature * sinasaiang * sinasaiang
        	let rayheightprime = odz ? tc.rayHeight : tc.objectDistance * asaprime
        	let objectdistanceprime = ((tc.radiusOfCurvature *
                                	    sin(asaprime + iang)) *
                                    	    cot(asadoubleprime)) + sagitta

		tcp.objectDistance = objectdistanceprime
		tcp.rayHeight = rayheightprime
		tcp.axisSlopeAngle = asadoubleprime

	    } else {

	    	//  Flat surface

		let rang = -(asin((tc.fromIndex / tc.toIndex))) * sin(tc.axisSlopeAngle)

		tcp.objectDistance = tc.objectDistance * ((tc.toIndex *
		    	    cos(-rang)) / (tc.fromIndex * cos(tc.axisSlopeAngle)))
		tcp.axisSlopeAngle = -rang
	    }
	}

    	return tcp
    }

    /*  traceLens propagates a ray of the specified spectral line
    	and axial incidence through a design and returns the final
	object distance and axis slope angle, from which we can
	compute the aberrations of the design.  */

    func traceLens(_ design: [SurfaceProperty], _ aperture: Double,
    	    	   _ line: spectralLine, _ incidence: AxialIncidence) ->
		    	(odist: Double, asa: Double) {
    	var ctx = TraceContext(axialIncidence: incidence,
	    	    	       radiusOfCurvature: 0,
			       objectDistance: 0,
			       rayHeight: aperture / 2,
			       axisSlopeAngle: 0,
			       fromIndex: 1,
			       toIndex: 0)

	for surf in design {
	    ctx.radiusOfCurvature = surf.curvatureRadius
	    ctx.toIndex = surf.indexOfRefraction
	    if (ctx.toIndex > 1) {
	    	ctx.toIndex += ((spectralLine.D.rawValue -
		    line.rawValue) / (spectralLine.C.rawValue - spectralLine.F.rawValue)) *
		    ((surf.indexOfRefraction - 1) / surf.dispersion)
	    }
	    ctx = transitSurface(ctx)
	    ctx.fromIndex = ctx.toIndex
	    ctx.objectDistance -= surf.edgeThickness 	// edgeThickness == 0 for last surface
	}

	return (ctx.objectDistance, ctx.axisSlopeAngle)
    }

    /*  The evaluateDesign function performs a ray trace on a given design
	with a specified clear aperture and returns a DesignEvaluation
	which includes the results for the D line and calculation of
	spherical aberration, coma, and chromatic aberration, along with
	the conventional acceptable upper bounds for these quantities.  */

    struct DesignEvaluation {
    	var dMarginalOD = 0.0	    // D Marginal ray
	var dMarginalSA = 0.0

    	var dParaxialOD = 0.0	    // D Paraxial ray
	var dParaxialSA = 0.0

	//  Computed aberrations of design
        var longitudinalSphericalAberration = 0.0
        var offenseAgainstSineCondition = 0.0
        var axialChromaticAberration = 0.0

        //  Acceptable maxima for aberrations
        var maxLongitudinalSphericalAberration = 0.0
        var maxOffenseAgainstSineCondition = 0.0025
        var maxAxialChromaticAberration = 0.0
    }

    func evaluateDesign(_ design: [SurfaceProperty], _ aperture: Double) -> DesignEvaluation {
    	var de = DesignEvaluation()

	//  D marginal ray
	(de.dMarginalOD, de.dMarginalSA) = traceLens(design, aperture,
	    spectralLine.D, AxialIncidence.marginalRay)

    	//  D paraxial ray
	(de.dParaxialOD, de.dParaxialSA) = traceLens(design, aperture,
	    spectralLine.D, AxialIncidence.paraxialRay)

	//  C marginal ray
	let (cMarginalOD, _) = traceLens(design, aperture,
	    spectralLine.C, AxialIncidence.marginalRay)

	//  F marginal ray
	let (fMarginalOD, _) = traceLens(design, aperture,
	    spectralLine.F, AxialIncidence.marginalRay)

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
	    (sin(de.dMarginalSA) * de.dMarginalOD)

        /*  The axial chromatic aberration is the distance between
            where marginal rays in the C and F lines come to focus.  */
        de.axialChromaticAberration = fMarginalOD - cMarginalOD

        //  Compute maximum acceptable values for each aberration

        /*  Maximum longitudinal spherical aberration, which is
            also the maximum for axial chromatic aberration.  This
            is computed for the D line.  */
        let sin_dm_sa = sin(de.dMarginalSA)
        de.maxLongitudinalSphericalAberration = 0.0000926 / (sin_dm_sa * sin_dm_sa)
	de.maxAxialChromaticAberration = de.maxLongitudinalSphericalAberration // Same criterion

	return de
    }

    /*  The evaluationReport function takes a DesignEvaluation and returns a
	list of strings containing a primate-readable version of the evaluation.
	For the purposes of the benchmark, this function is not timed; it serves
	only to create textual results which can be compared against those from
	the reference implementation.  */

    func evaluationReport(_ de: DesignEvaluation) -> [String] {
    	let mp = "    (Maximum permissible):              %16.11f"

    	return [
	    String(format: "   Marginal ray        %16.11f  %14.11f",
	    	de.dMarginalOD, de.dMarginalSA),
	    String(format: "   Paraxial ray        %16.11f  %14.11f",
	    	de.dParaxialOD, de.dParaxialSA),
	    String(format: "Longitudinal spherical aberration:      %16.11f",
	    	de.longitudinalSphericalAberration),
	    String(format: mp,
	    	de.maxLongitudinalSphericalAberration),
	    String(format: "Offense against sine condition (coma):  %16.11f",
	    	de.offenseAgainstSineCondition),
	    String(format: mp,
	    	de.maxOffenseAgainstSineCondition),
	    String(format: "Axial chromatic aberration:             %16.11f",
	    	de.axialChromaticAberration),
	    String(format: mp,
	    	de.maxAxialChromaticAberration)
	       ]
    	}

    /*  The validateResults function compares a primate-readable report
	from evaluationReport with the archival results from the reference
	implementation (which all language implementations must reproduce
	character-by-character [apart from trivia such as end of line
	conventions and trailing white space]).  It returns a Bool
	indicating whether the results. This function is not timed in the
	benchmark.  */

    func validateResults(_ er: [String]) -> Bool {
	/* Reference results.  These happen to be derived from a run
	   on Microsoft Quick BASIC on the IBM PC/AT.  */
    	let expectedResults = [
            "   Marginal ray          47.09479120920   0.04178472683",
            "   Paraxial ray          47.08372160249   0.04177864821",
            "Longitudinal spherical aberration:        -0.01106960671",
            "    (Maximum permissible):                 0.05306749907",
            "Offense against sine condition (coma):     0.00008954761",
            "    (Maximum permissible):                 0.00250000000",
            "Axial chromatic aberration:                0.00448229032",
            "    (Maximum permissible):                 0.05306749907"
	    	    	      ]
	var j = 0
	var errs = 0
	for s in er {
	    if s != expectedResults[j] {
		print("Validation failed on line \(j + 1):")
		print("  Expected: \"\(s)\"")
		print("  Received: \"\(expectedResults[j])\"")
	    	errs += 1
	    }
	    j += 1
	}
	return errs == 0
    }

    func runBenchmark(_ iterations: Int) {
    	var de: DesignEvaluation
    	for _ in 1...iterations - 1 {
	    de = evaluateDesign(wyldLens, wyldClearAperture)
	}
	//  I do it this way to avoid uninitialised structure error
	de = evaluateDesign(wyldLens, wyldClearAperture)
	let dr = evaluationReport(de)
	let ok = validateResults(dr)
    	if !ok {
    	    print("Error(s) detected in results.  This is VERY SERIOUS.")
    	}
    }

    /*	The argument specifies the number for iterations to
    	run.  For archival purposes, set the iteration count
	to achieve a run time of around five minutes.  For the
	runs on my machine, I used an iteration count of 66561695.  */

    runBenchmark(100000)

