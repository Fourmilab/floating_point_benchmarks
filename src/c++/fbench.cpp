    /*  John Walker's Floating Point Benchmark, derived from...

            Marinchip Interactive Lens Design System
                   John Walker   December 1980

                       by John Walker
                   http://www.fourmilab.ch/

        This program may be used, distributed, and modified
        freely as long as the origin information is preserved.

        This is a complete optical design raytracing algorithm,
        stripped of its user interface and recast into C++. It
        not only determines execution speed on an extremely
        floating point (including trig function) intensive
        real-world application, it checks accuracy on an
        algorithm that is exquisitely sensitive to errors. The
        performance of this program is typically far more
        sensitive to changes in the efficiency of the
        trigonometric library routines than the average floating
        point program.

        Implemented in November 2017 by John Walker.  */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

    using namespace std;

#   define Throw(exception, message) { \
        ostringstream em; \
        em << __FILE__ << "(" << __LINE__ << "): " << message; \
        throw(exception(em.str())); \
    }

    /*  You can configure this program to use floating point
        numbers of different precisions by defining the following
        symbols to be 1.
            (none)          C++ "double"
            LONG_DOUBLE     C++ "long double"
            FLOAT128        GCC's 128 bit floating point type
            FLOAT_MPFR      MPFR multiple-precision package:
                            MPFR_PRECISION sets mantissa precision in bits
        Note that the meaning of these symbols is dependent upon
        the compiler and the architecture on which it is running.
        With GCC 4.10.0 on an X86_64 machine, all precisions
        yielded identical results in the computed output.  */

#if FLOAT128
#include <quadmath.h>
    typedef __float128 Real;
#define sin(x)  sinq(x)
#define asin(x) asinq(x)
#define cos(x)  cosq(x)
#define tan(x)  tanq(x)
    //  Implement printing __float128 for debug output
    ostream& operator<<(ostream &os, Real t) {
        char s[80];
        quadmath_snprintf(s, sizeof s, "%16.11Qf", t);
        os << s;
        return os;
    }
#elif FLOAT_MPFR
/*  This uses our local copy of the MPFR C++ binding of the MPFR
    library.  If you have a later version of this software on your
    system, you'll probably want to replace the follwing statement
    with:
        #include <mpreal.h>
*/
#include "mpfrc++-3.6.3/mpreal.h"
    using namespace mpfr;
    typedef mpreal Real;
#   define Provides_cot
#   ifndef MPFR_PRECISION
#       define MPFR_PRECISION 128
#   endif
#elif LONG_DOUBLE
    typedef long double Real;
#   define RealFormat "Lf"
#else
    typedef double Real;
#   define RealFormat "f"
#endif

#ifndef Provides_cot
    //  Define cot() in terms of tan()
#   define cot(x) (1.0 / tan(x))
#endif

    /*  Wavelengths of standard spectral lines in Angstroms
              (Not all are used in this program)  */

    typedef Real Wavelength;

    class SpectralLine {
    public:
#if !FLOAT_MPFR
        const static Wavelength A = 7621.0,
                                B = 6869.955,
                                C = 6562.816,
                                D = 5895.944,
                                E = 5269.557,
                                F = 4861.344,
                                Gprime = 4340.477,
                                H = 3968.494;
#else
        //  MPFR doesn't allow us to initialise these in the class
        static Wavelength A, B, C, D, E, F, Gprime, H;
#endif
    };
#if FLOAT_MPFR
    Wavelength SpectralLine::A = 7621.0,
               SpectralLine::B = 6869.955,
               SpectralLine::C = 6562.816,
               SpectralLine::D = 5895.944,
               SpectralLine::E = 5269.557,
               SpectralLine::F = 4861.344,
               SpectralLine::Gprime = 4340.477,
               SpectralLine::H = 3968.494;
#endif

    /*  A surface describes the boundary between two components
        in the Design.  */

    class Surface {
    public:
        Real curvature_Radius,
             index_Of_Refraction,
             dispersion,
             edge_Thickness;

        //  Constructor
        Surface(Real r, Real i, Real d, Real e) {
            curvature_Radius = r;
            index_Of_Refraction = i;
            dispersion = d;
            edge_Thickness = e;
        }

        //  Dump a surface for debugging
        void show(ostream &os) {
            os << setw(8) << curvature_Radius <<
                  setw(8) << index_Of_Refraction <<
                  setw(8) << dispersion <<
                  setw(8) << edge_Thickness << endl;
        }
    };

    /*  A Design is the specification of the optical assembly
        to be evaluated.  */

    class Design {
    private:
        const static int MaxSurfaces = 10;

    public:
        Real clearAperture;
        unsigned int nSurfaces;
        Surface *surf[MaxSurfaces];

        Design(void) {
            nSurfaces = 0;
        }

        Design(Real ca, unsigned int ns) {
            clearAperture = ca;
            if (ns >= MaxSurfaces) {
                Throw(out_of_range, "number of surfaces exceeds MaxSurfaces");
            }
            nSurfaces = ns;
        }

        ~Design() {
            //  Delete surfaces, which were dynamically allocated
            for (unsigned int i = 0; i < nSurfaces; i++) {
                if (surf[i] != NULL) {
                    delete surf[i];
                }
            }
        }

        void setSurf(int n, Surface *s) {
            if (n >= MaxSurfaces) {
                Throw(out_of_range, "surface index exceeds MaxSurfaces");
            }
            surf[n] = s;
        }

        //  Dump a design for debugging
        void show(ostream &os) {
            os << "Clear aperture: " << clearAperture << endl;
            os << " Surf      Radius           Index        " <<
                  "Dispersion       Edge Thick" << endl;
            for (unsigned int i = 0; i < nSurfaces; i++) {
                os << "  " << i << ": ";
                surf[i]->show(os);
            }
        }
    };

    enum AxialIncidence { Marginal_Ray, Paraxial_Ray };

    class TraceContext {
    private:
        Design *d;
        AxialIncidence axial_incidence;
        Wavelength line;
        unsigned int cSurf;
        Real radius_of_curvature,
             object_distance,
             ray_height,
             axis_slope_angle,
             from_index,
             to_index;

        /*  Transit a surface.  Returns true if this was the last
            surface in the design.  */
        bool transitSurface(void);

    public:

        //  Constructor
        TraceContext(Design &des, Wavelength w, AxialIncidence ai) {
            set(des, w, ai);
        }

        //  Reset a TraceContext to new values as by the constructor
        void set(Design &des, Wavelength w, AxialIncidence ai) {
            d = &des;
            axial_incidence = ai;
            line = w;
            cSurf = 0;
            radius_of_curvature = object_distance =
                axis_slope_angle = to_index = 0;
            ray_height = d->clearAperture / 2;
            from_index = 1;
        }

        //  Trace a spectral line through the design
        void traceLine(Real &od, Real &sa);
        void traceLine(Real &od);

        //  Dump a TraceContext for debugging
        void show(ostream &os) {
            streamsize op = os.precision(13);
            os << "OD: " << object_distance <<
                  "  SA: " << axis_slope_angle <<
                  "  ROC: " << radius_of_curvature <<
                  "  RH: " << ray_height <<
                  "  FI: " << from_index <<
                  "  TI: " << to_index <<
                  endl;
            os.precision(op);
        }
    };

    /*  TraceContext::transitSurface propagates a ray through a
        Design.

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

        transitSurface returns false when the last surface has been
        traversed and true if more surfaces remain to be traced.
    */


    bool TraceContext::transitSurface(void) {

        //  Set context variables from current surface

        radius_of_curvature = d->surf[cSurf]->curvature_Radius;
        to_index = d->surf[cSurf]->index_Of_Refraction;
        if (to_index > 1) {
            to_index += ((SpectralLine::D - line) /
                (SpectralLine::C - SpectralLine::F)) *
                ((d->surf[cSurf]->index_Of_Refraction - 1) / d->surf[cSurf]->dispersion);
        }

        if (axial_incidence == Paraxial_Ray) {

            //  Paraxial ray

            if (radius_of_curvature != 0) {

                //  Curved surface

                const bool odz = object_distance == 0;
                const Real asaprime = odz ? 0 : axis_slope_angle;
                const Real iangsin = odz ? ray_height / radius_of_curvature :
                                ((object_distance - radius_of_curvature) /
                                 radius_of_curvature) * axis_slope_angle;
                const Real rangsin = (from_index / to_index) * iangsin;
                const Real asadoubleprime = asaprime + iangsin - rangsin;
                const Real rayheightprime = odz ? ray_height :
                                        object_distance * asaprime;
                const Real objectdistanceprime = rayheightprime / asadoubleprime;

                object_distance = objectdistanceprime;
                ray_height = rayheightprime;
                axis_slope_angle = asadoubleprime;

            } else {

                //  Flat surface

                object_distance = object_distance * (to_index / from_index);
                axis_slope_angle = axis_slope_angle * (from_index / to_index);
            }
         } else {

            //  Marginal ray

            if (radius_of_curvature != 0) {

                //  Curved surface

                const bool odz = object_distance == 0;
                const Real asaprime = odz ? 0 : axis_slope_angle;
                const Real iangsin = odz ? ray_height / radius_of_curvature :
                                ((object_distance - radius_of_curvature) /
                                 radius_of_curvature) * sin(axis_slope_angle);
                const Real iang = asin(iangsin);
                const Real rangsin = (from_index / to_index) * iangsin;
                const Real asadoubleprime = asaprime + iang - asin(rangsin);
                const Real sinasaiang = sin((asaprime + iang) / 2);
                const Real sagitta = 2 * radius_of_curvature * sinasaiang * sinasaiang;
                const Real rayheightprime = odz ? ray_height :
                                        object_distance * asaprime;
                const Real objectdistanceprime = ((radius_of_curvature *
                                          sin(asaprime + iang)) *
                                          cot(asadoubleprime)) + sagitta;

                object_distance = objectdistanceprime;
                ray_height = rayheightprime;
                axis_slope_angle = asadoubleprime;

            } else {

                //  Flat surface

                const Real rang = -(asin((from_index / to_index))) *
                              sin(axis_slope_angle);

                object_distance = object_distance * ((to_index *
                            cos(-rang)) / (from_index *
                            cos(axis_slope_angle)));
                axis_slope_angle = -rang;
            }
        }

        from_index = to_index;
        object_distance -= d->surf[cSurf]->edge_Thickness;
        cSurf++;

        return cSurf >= d->nSurfaces;
    }

    void TraceContext::traceLine(Real &od, Real &sa) {
        do {
        } while (!transitSurface());
        od = object_distance;
        sa = axis_slope_angle;
    }

    void TraceContext::traceLine(Real &od) {
        do {
        } while (!transitSurface());
        od = object_distance;
    }

    /*  A DesignEvaluation provides tools to analyse designs.  It
        takes a design, traces rays through it in various wavelengths
        and axial incidences, and computes its aberrations compared to
        acceptable standards.  */

    class DesignEvaluation {
    private:
        Design *d;

        Real cMarginalOD;           // C marginal ray
        Real fMarginalOD;           // F marginal ray

        char received[8][80];       // Edited results of evaluation

    public:
        Real dMarginalOD;           // D Marginal ray
        Real dMarginalSA;

        Real dParaxialOD;           // D Paraxial ray
        Real dParaxialSA;

        //  Computed aberrations of design
        Real longitudinalSphericalAberration;
        Real offenseAgainstSineCondition;
        Real axialChromaticAberration;

        //  Acceptable maxima for aberrations
        Real maxLongitudinalSphericalAberration;
        Real maxOffenseAgainstSineCondition;
        Real maxAxialChromaticAberration;

        //  Construct a DesignEvaluation
        DesignEvaluation(Design &des) {
            d = &des;
            maxOffenseAgainstSineCondition = 0.0025;
        }

        //  Evaluate the design
        void evaluate(void);

        //  Edit the evaluation into a primate-readable report
        void report(void);

        //  Print the report
        void print(ostream &os) {
            for (int i = 0; i < 8; i++) {
                os << received[i] << endl;
            }
        }

        //  Validate the report
        unsigned int validate(ostream &os);
    };

    void DesignEvaluation::evaluate(void) {

        //  D marginal ray
        TraceContext tc(*d, SpectralLine::D, Marginal_Ray);
        tc.traceLine(dMarginalOD, dMarginalSA);

        //  D paraxial ray
        tc.set(*d, SpectralLine::D, Paraxial_Ray);
        tc.traceLine(dParaxialOD, dParaxialSA);

        //  C marginal ray
        tc.set(*d, SpectralLine::C, Marginal_Ray);
        tc.traceLine(cMarginalOD);

        //  F marginal ray
        tc.set(*d, SpectralLine::F, Marginal_Ray);
        tc.traceLine(fMarginalOD);

        //  Compute aberrations of the design

        /*  The longitudinal spherical aberration is just the
            difference between where the D line comes to focus
            for paraxial and marginal rays.  */
        longitudinalSphericalAberration = dParaxialOD - dMarginalOD;

        /*  The offense against the sine condition is a measure
            of the degree of coma in the design.  We compute it
            as the lateral distance in the focal plane between
            where a paraxial ray and marginal ray in the D line
            come to focus.  */
        offenseAgainstSineCondition = 1 - (dParaxialOD * dParaxialSA) /
            (sin(dMarginalSA) * dMarginalOD);

        /*  The axial chromatic aberration is the distance between
            where marginal rays in the C and F lines come to focus.  */
        axialChromaticAberration = fMarginalOD - cMarginalOD;

        //  Compute maximum acceptable values for each aberration

        /*  Maximum longitudinal spherical aberration, which is
            also the maximum for axial chromatic aberration.  This
            is computed for the D line.  */
        const Real sin_dm_sa = sin(dMarginalSA);
        maxLongitudinalSphericalAberration = 0.0000926 / (sin_dm_sa * sin_dm_sa);
        maxAxialChromaticAberration = maxLongitudinalSphericalAberration; // Same criterion
    }

    void DesignEvaluation::report(void) {
#if FLOAT128
        /*  The libquadmath installation on the system with which I
            developed this code does not provide XXprintf variants
            that support formatting such numbers, but only its own
            quadmath_snprintf() which accepts a single format code.
            The following code formats the numbers and then edits
            them into the output as strings.  */
        const static char mp[] =
            "    (Maximum permissible):              %s",
                          ry[] = "%15s   %s  %s",
                          qf[] =
            "%16.11Qf";
        char b1[80], b2[80];
        quadmath_snprintf(b1, sizeof b1, "%21.11Qf", dMarginalOD);
        quadmath_snprintf(b2, sizeof b2, "%14.11Qf", dMarginalSA);
        snprintf(received[0], 80, ry, "Marginal ray", b1, b2);
        quadmath_snprintf(b1, sizeof b1, "%21.11Qf", dParaxialOD);
        quadmath_snprintf(b2, sizeof b2, "%14.11Qf", dParaxialSA);
        snprintf(received[1], 80, ry, "Paraxial ray", b1, b2);
#       define Qe(x) ((void) quadmath_snprintf(b1, sizeof b1, qf, x), b1)
        snprintf(received[2], 80,
           "Longitudinal spherical aberration:      %s",
           Qe(longitudinalSphericalAberration));
        snprintf(received[3], 80, mp, Qe(maxLongitudinalSphericalAberration));
        snprintf(received[4], 80,
           "Offense against sine condition (coma):  %s",
           Qe(offenseAgainstSineCondition));
        snprintf(received[5], 80, mp, Qe(maxOffenseAgainstSineCondition));
        snprintf(received[6], 80,
           "Axial chromatic aberration:             %s",
           Qe(axialChromaticAberration));
        snprintf(received[7], 80, mp, Qe(maxAxialChromaticAberration));
#       undef  Qe
#elif FLOAT_MPFR
        /*  The MPFR C++ package does not allow its mpreal values to be
            edited by XXprintf functions, but provides a toString method
            which accepts XXprintf-like format codes.  The following code
            uses this method to format the output of the evaluation into
            the received[] array.  */
        const static char mp[] =
            "    (Maximum permissible):              %s",
                          ry[] = "%15s   %s  %s";
        snprintf(received[0], 80, ry, "Marginal ray",
            dMarginalOD.toString("%21.11Rf").c_str(),
            dMarginalSA.toString("%14.11Rf").c_str());
        snprintf(received[1], 80, ry, "Paraxial ray",
            dParaxialOD.toString("%21.11Rf").c_str(),
            dParaxialSA.toString("%14.11Rf").c_str());
#       define Qe(x) (x.toString("%16.11Rf").c_str())
        snprintf(received[2], 80,
           "Longitudinal spherical aberration:      %s",
           Qe(longitudinalSphericalAberration));
        snprintf(received[3], 80, mp, Qe(maxLongitudinalSphericalAberration));
        snprintf(received[4], 80,
           "Offense against sine condition (coma):  %s",
           Qe(offenseAgainstSineCondition));
        snprintf(received[5], 80, mp, Qe(maxOffenseAgainstSineCondition));
        snprintf(received[6], 80,
           "Axial chromatic aberration:             %s",
           Qe(axialChromaticAberration));
        snprintf(received[7], 80, mp, Qe(maxAxialChromaticAberration));
#       undef  Qe
#else
        const static char mp[] =
            "    (Maximum permissible):              %16.11" RealFormat,
                          ry[] =
            "%15s   %21.11" RealFormat "  %14.11" RealFormat;

        snprintf(received[0], 80, ry, "Marginal ray", dMarginalOD, dMarginalSA);
        snprintf(received[1], 80, ry, "Paraxial ray", dParaxialOD, dParaxialSA);
        snprintf(received[2], 80,
           "Longitudinal spherical aberration:      %16.11" RealFormat,
           longitudinalSphericalAberration);
        snprintf(received[3], 80, mp, maxLongitudinalSphericalAberration);
        snprintf(received[4], 80,
           "Offense against sine condition (coma):  %16.11" RealFormat,
           offenseAgainstSineCondition);
        snprintf(received[5], 80, mp, maxOffenseAgainstSineCondition);
        snprintf(received[6], 80,
           "Axial chromatic aberration:             %16.11" RealFormat,
           axialChromaticAberration);
        snprintf(received[7], 80, mp, maxAxialChromaticAberration);
#endif
    }

    unsigned int DesignEvaluation::validate(ostream &os) {
        unsigned int errors = 0;
        /*  Reference results.  These happen to be derived from
            a run on Microsoft Quick BASIC on the IBM PC/AT.  */
        const static char * const expected[] = {
            "   Marginal ray          47.09479120920   0.04178472683",
            "   Paraxial ray          47.08372160249   0.04177864821",
            "Longitudinal spherical aberration:        -0.01106960671",
            "    (Maximum permissible):                 0.05306749907",
            "Offense against sine condition (coma):     0.00008954761",
            "    (Maximum permissible):                 0.00250000000",
            "Axial chromatic aberration:                0.00448229032",
            "    (Maximum permissible):                 0.05306749907"
        };

        for (int i = 0; i < 8; i++) {
           if (strcmp(received[i], expected[i]) != 0) {
              os << "Error in results on line " << i + 1 << "..." << endl;
              os << "Expected:  \"" << expected[i] << "\"" << endl;
              os << "Received:  \"" << received[i] << "\"" << endl;
              os << "(Errors)    ";
              int k = strlen(expected[i]);
              for (int j = 0; j < k; j++) {
                 os << ((expected[i][j] == received[i][j]) ? ' ' : '^');
                 if (expected[i][j] != received[i][j]) {
                    errors++;
                 }
              }
              os << endl;
           }
        }

        return errors;
    }

    int main(int argc, char *argv[]) {

        long iterations = 1000000;

        if (argc > 1) {
            iterations = atol(argv[1]);
        }

#if FLOAT_MPFR
        mpreal::set_default_prec(MPFR_PRECISION);
#endif

        /*  The test case used in this program is the design
            for a 4 inch f/12 achromatic telescope objective
            used as the example in Wyld's classic work on ray
            tracing by hand, given in Amateur Telescope Making,
            Volume 3 (Volume 2 in the 1996 reprint edition).  */

        Design WyldLens(4.0, 4);
                               //        CurRad  Index   Disp  Edge
        WyldLens.setSurf(0, new Surface(  27.05, 1.5137, 63.6, 0.52 ));
        WyldLens.setSurf(1, new Surface( -16.68, 1.0,     0.0, 0.138));
        WyldLens.setSurf(2, new Surface( -16.68, 1.6164, 36.7, 0.38 ));
        WyldLens.setSurf(3, new Surface( -78.1,  1.0,     0.0, 0.0  ));
//WyldLens.show(cout);

        DesignEvaluation de(WyldLens);
        for (long l = 0; l < iterations; l++) {
            de.evaluate();
        }
        de.report();
//de.print(cout);
        unsigned int errors;
        if ((errors = de.validate(cout)) > 0) {
            cout << errors << " error" << (errors > 1 ? "s" : "") <<
                " in results.  This is VERY SERIOUS." << endl;
        } else {
           cout << "No errors in results." << endl;
        }

        return 0;
    }
