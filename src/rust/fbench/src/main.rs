/* John Walker's Floating Point Benchmark, derived from...

        Marinchip Interactive Lens Design System
               John Walker   December 1980

                   by John Walker
               http://www.fourmilab.ch/

    This program may be used, distributed, and modified freely as
    long as the origin information is preserved.

    This is a complete optical design raytracing algorithm, stripped of
    its user interface and recast into Rust. It not only determines
    execution speed on an extremely floating point (including trig
    function) intensive real-world application, it checks accuracy on
    an algorithm that is exquisitely sensitive to errors. The
    performance of this program is typically far more sensitive to
    changes in the efficiency of the trigonometric library routines
    than the average floating point program.

    Implemented in October 2014 by John Walker.
    Extensively revised in January 2021 to cope with changes in the
    Rust language as of version 1.43.0.
*/

extern crate time;

fn main() {

    /*  The test case used in this program is the design for a 4 inch
        f/12 achromatic telescope objective used as the example in Wyld's
        classic work on ray tracing by hand, given in Amateur Telescope
        Making, Volume 3 (Volume 2 in the 1996 reprint edition).  */

    struct SurfaceProperty {
        curvature_radius: f64,
        index_of_refraction: f64,
        dispersion: f64,
        edge_thickness: f64
    }

    const WYLD_CLEAR_APERTURE : f64 = 4.0;

    const WYLD_LENS : [ SurfaceProperty;4 ] = [
        SurfaceProperty {
            curvature_radius: 27.05,
            index_of_refraction: 1.5137,
            dispersion: 63.6,
            edge_thickness: 0.52
        },

        SurfaceProperty {
            curvature_radius: -16.68,
            index_of_refraction: 1.0,
            dispersion: 0.0,
            edge_thickness: 0.138
        },

        SurfaceProperty {
            curvature_radius: -16.68,
            index_of_refraction: 1.6164,
            dispersion: 36.7,
            edge_thickness: 0.38
        },

        SurfaceProperty {
            curvature_radius: -78.1,
            index_of_refraction: 1.0,
            dispersion: 0.0,
            edge_thickness: 0.0
        }
    ];

    //  Wavelengths of standard spectral lines in Angstroms

    const C_LINE: f64 = 6562.816;
    const D_LINE: f64 = 5895.944;
    const F_LINE: f64 = 4861.344;

    /*  To avoid the ugly syntax of math functions being written
        as methods of numbers, we define syntactic sugar functions
        to wrap them.  */

    fn asin(x: f64) -> f64 { x.asin() }
    fn cos(x: f64) -> f64 { x.cos() }
    fn sin(x: f64) -> f64 { x.sin() }
    fn cot(x: f64) -> f64 { 1.0 / x.tan() }

    /*  The heavy lifting of the ray trace is performed by the
        transit_surface function.  It is called with a
        TraceContext consisting of the following items.  The
        first is an integer indicating whether the ray is
        marginal or paraxial and the balance are all real
        numbers.

        ray                     This is either MARGINAL_RAY or
                                PARAXIAL_RAY depending on the
                                offset of the initial ray from
                                the axis.

        radius                  Radius of curvature of surface
                                being crossed.  If 0, surface is
                                plane.

        odist                   Distance of object focus from
                                lens vertex.  If 0, incoming
                                rays are parallel and
                                the following must be specified:

        rheight                 Height of ray from axis.  Only
                                relevant if odist == 0

        asa                     Angle incoming ray makes with axis
                                at intercept

        fromindex               Refractive index of medium being left


        toindex                 Refractive index of medium being
                                entered.

        The result of a call on transit_surface is a TraceContext
        with the following components modified to reflect the geometry
        of the ray as it exits the surface.

        odist                   Distance from vertex to object focus
                                after refraction.

        rheight                 Height of ray from axis.

        asa                     Angle incoming ray makes with axis
                                at intercept after refraction.  */

    const MARGINAL_RAY: u8 = 0;
    const PARAXIAL_RAY: u8 = 1;

    struct TraceContext {
        ray: u8,
        radius: f64,
        odist: f64,
        rheight: f64,
        asa: f64,
        fromindex: f64,
        toindex: f64
    }

    /*  We break the function down into four cases depending
        upon whether the ray is marginal or paraxial and
        whether the surface is flat or curved.  */

    //  Paraxial ray with flat radius of curvature

    fn transit_surface_pf(c: TraceContext) -> TraceContext
    {   TraceContext {
            ray: c.ray,
            radius: 0.0,
            odist: (c.odist * (c.toindex / c.fromindex)),
            rheight: c.rheight,
            asa: (c.asa * (c.fromindex / c.toindex)),
            fromindex: c.fromindex,
            toindex: c.toindex
        }
    }

    //  Paraxial ray with curved surface

    fn transit_surface_pc(c: TraceContext) -> TraceContext
    {
        let odz = c.odist == 0.0;

        let asaprime = if odz { 0.0 } else { c.asa };

        let iangsin = if odz { c.rheight / c.radius} else
                             { ((c.odist - c.radius) / c.radius) * c.asa };

        let rangsin = (c.fromindex / c.toindex) * iangsin;

        let asadoubleprime = asaprime + iangsin - rangsin;

        let rayheightprime = if odz { c.rheight } else
                                    { c.odist * asaprime };

        let objectdistanceprime = rayheightprime / asadoubleprime;

        TraceContext {
            ray: c.ray,
            radius: c.radius,
            odist: objectdistanceprime,
            rheight: rayheightprime,
            asa: asadoubleprime,
            fromindex: c.fromindex,
            toindex: c.toindex
        }
    }

    //  Marginal ray with flat radius of curvature

    fn transit_surface_mf(c: TraceContext) -> TraceContext
    {
        let rang = -(asin(c.fromindex / c.toindex)) *
                    sin(c.asa);

        TraceContext {
            ray: c.ray,
            radius: 0.0,
            odist: (c.odist * ((c.toindex * cos(-rang)) /
                   (c.fromindex * cos(c.asa)))),
            rheight: c.rheight,
            asa: -rang,
            fromindex: c.fromindex,
            toindex: c.toindex
        }
    }

    //  Marginal ray with curved surface

    fn transit_surface_mc(c: TraceContext) -> TraceContext
    {
        let odz = c.odist == 0.0;

        let asaprime = if odz { 0.0 } else { c.asa };

        let iangsin = if odz { c.rheight / c.radius } else
                             { ((c.odist - c.radius) /
                                c.radius) * sin(c.asa) };
        let iang = asin(iangsin);
        let rangsin = (c.fromindex / c.toindex) * iangsin;
        let asadoubleprime = asaprime + iang - asin(rangsin);
        let sinasaiang = sin((asaprime + iang) / 2.0);
        let sagitta = 2.0 * c.radius * sinasaiang * sinasaiang;
        let ray_heightprime = if odz { c.rheight } else
                                     { c.odist * asaprime };
        let object_distanceprime = ((c.radius * sin(asaprime + iang)) *
                                    cot(asadoubleprime)) + sagitta;

        TraceContext {
            ray: c.ray,
            radius: c.radius,
            odist: object_distanceprime,
            rheight: ray_heightprime,
            asa: asadoubleprime,
            fromindex: c.fromindex,
            toindex: c.toindex
        }
    }

    /*  transit_surface is a wrapper which calls one of the four
        functions for the various cases.  We could write the
        whole thing as a doubly-nested conditional, but this makes
        it easier to read.  */

    fn transit_surface(c: TraceContext) -> TraceContext
    {
        if c.ray == MARGINAL_RAY {
            if c.radius == 0.0 {
                transit_surface_mf(c)
            } else {
                transit_surface_mc(c)
           }
        } else { // c.ray == PARAXIAL_RAY
            if c.radius == 0.0 {
                transit_surface_pf(c)
            } else {
                transit_surface_pc(c)
           }
        }
    }

    /*  The trace_line function walks a given ray through the
        surfaces of a design, invoking transit_surface for each
        surface.  It is called with the lens design, the
        spectral line to be traced, and a TraceContext
        as described above for transit_surface.  It walks through
        the design list and uses transit_surface to follow the ray
        through the surfaces.  */

    fn trace_line(surf: &[SurfaceProperty], n: usize, sline: f64, c: TraceContext) ->
        (&[SurfaceProperty], usize, f64, TraceContext) {
        if n < surf.len() {
            let to_index_prime = if surf[n].index_of_refraction > 1.0 {
                                        surf[n].index_of_refraction +
                                        ((D_LINE - sline) /
                                        (C_LINE - F_LINE)) *
                                        ((surf[n].index_of_refraction - 1.0) /
                                        surf[n].dispersion)
                                    } else {
                                        surf[n].index_of_refraction
                                    };
            let cp = transit_surface(TraceContext{ray: c.ray,
                        radius: surf[n].curvature_radius,
                        odist: c.odist, rheight: c.rheight, asa: c.asa,
                        fromindex: c.fromindex, toindex: to_index_prime});

            let (_, _, _, cprime) = trace_line(surf, n + 1, sline,
                TraceContext{ray: cp.ray, radius: cp.radius,
                             odist: cp.odist - surf[n].edge_thickness,
                             rheight: cp.rheight, asa: cp.asa,
                             fromindex: cp.toindex, toindex: 0.0});

            (surf, n + 1, sline, cprime)

        } else {
            (surf, n + 1, sline, c)
        }
    }

    /*  The trace_lens function is a little bit of syntactic sugar
        to simplify invoking trace_line.  It creates an initial
        TraceContext for trace_line which specifies the axial_incidence
        and ray_height (computed from the clear aperture of the
        design), invokes trace_line with the supplied design and
        spectral line, then returns a tuple containing the
        object_distance and axis_slope_angle resulting from the ray
        trace.  */

    fn trace_lens(surf: &[SurfaceProperty], sline: f64, aperture: f64, axial_incidence: u8) ->
        (f64, f64)
    {
        let (_, _, _, ctx) = trace_line(surf, 0, sline,
            TraceContext{ray: axial_incidence, radius: 0.0, odist: 0.0,
                         rheight: aperture / 2.0, asa: 0.0,
                         fromindex: 1.0, toindex: 0.0});
        (ctx.odist, ctx.asa)
    }

    /*  The evaluate_design function performs a ray trace on a given design
        with a specified clear aperture and returns a DesignEvaluation
        which includes the results for the D line and calculation of
        spherical aberration, coma, and chromatic aberration, along with
        the conventional acceptable upper bounds for these quantities.  */

    struct DesignEvaluation {
        //  Results from ray trace
        marginal_ray_od: f64, marginal_ray_sa: f64,
        paraxial_ray_od: f64, paraxial_ray_sa: f64,
        longitudinal_spherical_aberration: f64,
        offense_against_sine_condition: f64,
        axial_chromatic_aberration: f64,

        //  Acceptable maxima for aberrations
        max_longitudinal_spherical_aberration: f64,
        max_offense_against_sine_condition: f64,
        max_axial_chromatic_aberration: f64
    }

    fn evaluate_design(surf: &[SurfaceProperty], aperture: f64) -> DesignEvaluation
    {

        let d_marginal = trace_lens(surf, D_LINE, aperture, MARGINAL_RAY);
        let (d_m_od, d_m_sa) = d_marginal;

        let d_paraxial = trace_lens(surf, D_LINE, aperture, PARAXIAL_RAY);
        let (d_p_od, d_p_sa) = d_paraxial;

        let c_marginal = trace_lens(surf, C_LINE, aperture, MARGINAL_RAY);
        let (c_m_od, _) = c_marginal;

        let f_marginal = trace_lens(surf, F_LINE, aperture, MARGINAL_RAY);
        let (f_m_od, _) = f_marginal;

        //  Compute aberrations of the design

        /*  The longitudinal spherical aberration is just the
            difference between where the D line comes to focus
            for paraxial and marginal rays.  */
        let aberr_lspher = d_p_od - d_m_od;

        /*  The offense against the sine condition is a measure
            of the degree of coma in the design.  We compute it
            as the lateral distance in the focal plane between
            where a paraxial ray and marginal ray in the D line
            come to focus.  */
        let aberr_osc = 1.0 - (d_p_od * d_p_sa) / (sin(d_m_sa) * d_m_od);

        /*  The axial chromatic aberration is the distance between
            where marginal rays in the C and F lines come to focus.  */
        let aberr_lchrom = f_m_od - c_m_od;

        //  Compute maximum acceptable values for each aberration

        /*  Maximum longitudinal spherical aberration, which is
            also the maximum for axial chromatic aberration.  This
            is computed for the D line.  */
        let sin_dm_sa = sin(d_m_sa);
        let max_lspher = 0.0000926 / (sin_dm_sa * sin_dm_sa);

        DesignEvaluation {
            //  Results from ray trace
            marginal_ray_od: d_m_od, marginal_ray_sa: d_m_sa,
            paraxial_ray_od: d_p_od, paraxial_ray_sa: d_p_sa,
            longitudinal_spherical_aberration: aberr_lspher,
            offense_against_sine_condition: aberr_osc,
            axial_chromatic_aberration: aberr_lchrom,

            //  Acceptable maxima for aberrations
            max_longitudinal_spherical_aberration: max_lspher,
            max_offense_against_sine_condition: 0.0025,
            max_axial_chromatic_aberration: max_lspher  // (NOT an error: see above)
        }
    }

    /*  The evaluation_report function takes a DesignEvaluation and prints
        a primate-readable report of its contents.  */

    fn evaluation_report(de: DesignEvaluation) -> ()
    {
        println!("   Marginal ray        {:16.11}  {:14.11}",
            de.marginal_ray_od, de.marginal_ray_sa);
        println!("   Paraxial ray        {:16.11}  {:14.11}",
            de.paraxial_ray_od, de.paraxial_ray_sa);
        println!("Longitudinal spherical aberration:      {:16.11}",
            de.longitudinal_spherical_aberration);
        println!("    (Maximum permissible):              {:16.11}",
            de.max_longitudinal_spherical_aberration);
        println!("Offense against sine condition (coma):  {:16.11}",
            de.offense_against_sine_condition);
        println!("    (Maximum permissible):              {:16.11}",
            de.max_offense_against_sine_condition);
        println!("Axial chromatic aberration:             {:16.11}",
            de.axial_chromatic_aberration);
        println!("    (Maximum permissible):              {:16.11}",
            de.max_axial_chromatic_aberration);
    }

    /*  The run_benchmark function runs the raytrace and analysis the
        specified number of times.  The DesignEvaluation from the last
        run is returned.  */

    fn run_benchmark(iterations: u32, surf: &[SurfaceProperty], aperture: f64) -> DesignEvaluation
    {
        for _ in 1 .. iterations {
            evaluate_design(surf, aperture);
        }
        evaluate_design(surf, aperture)
    }

    //  Actually run the benchmark, show results, and display time

    //  For official runs, adjust ITERATIONS so the benchmark
    //  runs about five minutes.
//    const ITERATIONS: u32 = 340_023_965;
    const ITERATIONS: u32 = 1_000_000;

    let start = time::precise_time_ns();
    let de = run_benchmark(ITERATIONS, &WYLD_LENS, WYLD_CLEAR_APERTURE);
    let fin = time::precise_time_ns();
    let lapse: f64 = ((fin - start) as f64) / 1_000_000_000.0;
    evaluation_report(de);

    println!("Time for {} iterations: {} seconds, {} \u{00B5}sec/iteration.",
        ITERATIONS, lapse, (lapse * 1_000_000.0) / (ITERATIONS as f64));
}
