/*

        John Walker's Floating Point Benchmark, derived from...

        Marinchip Interactive Lens Design System

                                     John Walker   December 1980

        By John Walker
           http://www.fourmilab.ch/

        This  program may be used, distributed, and modified freely as
        long as the origin information is preserved.

        This is a special hacked version of the main C benchmark
        which prints intermediate values at each surface transition
        of the trace of each spectral line.  It makes just one
        iteration, and is not intended for timing or accuracy
        checks.  Its sole purpose is to serve as a model when
        implementing the benchmark in other languages, as it
        allows comparing output at each step in the ray trace
        process to determine where the results diverged.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define cot(x) (1.0 / tan(x))

/*  Dump a double  */
#define DumpD(d) printf(#d " = %.11f\n", d)

#define TRUE  1
#define FALSE 0

#define max_surfaces 10

/*  Local variables  */

static short current_surfaces;
static short paraxial;

static double clear_aperture;

static double aberr_lspher;
static double aberr_osc;
static double aberr_lchrom;

static double max_lspher;
static double max_osc;
static double max_lchrom;

static double radius_of_curvature;
static double object_distance;
static double ray_height;
static double axis_slope_angle;
static double from_index;
static double to_index;

static double spectral_line[9];
static double s[max_surfaces][5];
static double od_sa[2][2];

static char outarr[8][80];         /* Computed output of program goes here */

int itercount;                     /* The iteration counter for the main loop
                                      in the program is made global so that
                                      the compiler should not be allowed to
                                      optimise out the loop over the ray
                                      tracing code. */

#define ITERATIONS 1

int niter = ITERATIONS;            /* Iteration counter */

static char *refarr[] = {          /* Reference results.  These happen to
                                      be derived from a run on Microsoft
                                      Quick BASIC on the IBM PC/AT. */

        "   Marginal ray          47.09479120920   0.04178472683",
        "   Paraxial ray          47.08372160249   0.04177864821",
        "Longitudinal spherical aberration:        -0.01106960671",
        "    (Maximum permissible):                 0.05306749907",
        "Offense against sine condition (coma):     0.00008954761",
        "    (Maximum permissible):                 0.00250000000",
        "Axial chromatic aberration:                0.00448229032",
        "    (Maximum permissible):                 0.05306749907"
};

/* The  test  case  used  in  this program is the  design for a 4 inch
   achromatic telescope  objective  used  as  the  example  in  Wyld's
   classic  work  on  ray  tracing by hand, given in Amateur Telescope
   Making, Volume 3.  */

static double testcase[4][4] = {
        {27.05, 1.5137, 63.6, 0.52},
        {-16.68, 1, 0, 0.138},
        {-16.68, 1.6164, 36.7, 0.38},
        {-78.1, 1, 0, 0}
};

/*            Calculate passage through surface

              If  the variable PARAXIAL is true, the trace through the
              surface will be done using the paraxial  approximations.
              Otherwise,  the normal trigonometric trace will be done.

              This routine takes the following inputs:

              RADIUS_OF_CURVATURE         Radius of curvature of surface
                                          being crossed.  If 0, surface is
                                          plane.

              OBJECT_DISTANCE             Distance of object focus from
                                          lens vertex.  If 0, incoming
                                          rays are parallel and
                                          the following must be specified:

              RAY_HEIGHT                  Height of ray from axis.  Only
                                          relevant if OBJECT.DISTANCE == 0

              AXIS_SLOPE_ANGLE            Angle incoming ray makes with axis
                                          at intercept

              FROM_INDEX                  Refractive index of medium being left

              TO_INDEX                    Refractive index of medium being
                                          entered.

              The outputs are the following variables:

              OBJECT_DISTANCE             Distance from vertex to object focus
                                          after refraction.

              AXIS_SLOPE_ANGLE            Angle incoming ray makes with axis
                                          at intercept after refraction.

*/

static void transit_surface() {
        double iang,               /* Incidence angle */
               rang,               /* Refraction angle */
               iang_sin,           /* Incidence angle sin */
               rang_sin,           /* Refraction angle sin */
               old_axis_slope_angle, sagitta;

        if (paraxial) {
           if (radius_of_curvature != 0.0) {
              if (object_distance == 0.0) {
                 axis_slope_angle = 0.0;
                 iang_sin = ray_height / radius_of_curvature;
              } else
                 iang_sin = ((object_distance -
                    radius_of_curvature) / radius_of_curvature) *
                    axis_slope_angle;

              rang_sin = (from_index / to_index) *
                 iang_sin;
              old_axis_slope_angle = axis_slope_angle;
              axis_slope_angle = axis_slope_angle +
                 iang_sin - rang_sin;
              if (object_distance != 0.0)
                 ray_height = object_distance * old_axis_slope_angle;
              object_distance = ray_height / axis_slope_angle;
              return;
           }
           object_distance = object_distance * (to_index / from_index);
           axis_slope_angle = axis_slope_angle * (from_index / to_index);
           return;
        }

        if (radius_of_curvature != 0.0) {
           if (object_distance == 0.0) {
              axis_slope_angle = 0.0;
              iang_sin = ray_height / radius_of_curvature;
           } else {
              iang_sin = ((object_distance -
                 radius_of_curvature) / radius_of_curvature) *
                 sin(axis_slope_angle);
           }
/* DumpD(ray_height); DumpD(iang_sin); DumpD(axis_slope_angle); */
           iang = asin(iang_sin);
           rang_sin = (from_index / to_index) *
              iang_sin;
           old_axis_slope_angle = axis_slope_angle;
           axis_slope_angle = axis_slope_angle +
              iang - asin(rang_sin);
           sagitta = sin((old_axis_slope_angle + iang) / 2.0);
           sagitta = 2.0 * radius_of_curvature*sagitta*sagitta;
           object_distance = ((radius_of_curvature * sin(
              old_axis_slope_angle + iang)) *
              cot(axis_slope_angle)) + sagitta;
           return;
        }

        rang = -asin((from_index / to_index) *
           sin(axis_slope_angle));
        object_distance = object_distance * ((to_index *
           cos(-rang)) / (from_index *
           cos(axis_slope_angle)));
        axis_slope_angle = -rang;
}

/*  Perform ray trace in specific spectral line  */

static void trace_line(line, ray_h)
int line;
double ray_h;
{
        int i;

printf("Tracing %c (%d) %s line at height %g:\n",
    'A' + (line - 1), line, paraxial ? "paraxial" : "marginal", ray_h);
        object_distance = 0.0;
        ray_height = ray_h;
        from_index = 1.0;

        for (i = 1; i <= current_surfaces; i++) {
           radius_of_curvature = s[i][1];
/* DumpD(radius_of_curvature); */
           to_index = s[i][2];
           if (to_index > 1.0)
              to_index = to_index + ((spectral_line[4] -
                 spectral_line[line]) /
                 (spectral_line[3] - spectral_line[6])) * ((s[i][2] - 1.0) /
                 s[i][3]);
/* DumpD(to_index); */
           transit_surface();
printf("  Surface %d  object_distance = %16.11f  axis_slope_angle = %16.11f\n",
    i, object_distance, axis_slope_angle);
           from_index = to_index;
           if (i < current_surfaces)
              object_distance = object_distance - s[i][4];
        }
printf("  Result     object_distance = %16.11f  axis_slope_angle = %16.11f\n",
    object_distance, axis_slope_angle);
}

/*  Initialise when called the first time  */

int main(argc, argv)
int argc;
char *argv[];
{
        int i, j, k, errors;
        double od_fline, od_cline;

        spectral_line[1] = 7621.0;       /* A */
        spectral_line[2] = 6869.955;     /* B */
        spectral_line[3] = 6562.816;     /* C */
        spectral_line[4] = 5895.944;     /* D */
        spectral_line[5] = 5269.557;     /* E */
        spectral_line[6] = 4861.344;     /* F */
        spectral_line[7] = 4340.477;     /* G'*/
        spectral_line[8] = 3968.494;     /* H */

        /* Load test case into working array */

        clear_aperture = 4.0;
        current_surfaces = 4;
        for (i = 0; i < current_surfaces; i++)
           for (j = 0; j < 4; j++)
              s[i + 1][j + 1] = testcase[i][j];

        /* Perform ray trace the specified number of times. */

        for (itercount = 0; itercount < niter; itercount++) {

           for (paraxial = 0; paraxial <= 1; paraxial++) {

              /* Do main trace in D light */

              trace_line(4, clear_aperture / 2.0);
              od_sa[paraxial][0] = object_distance;
              od_sa[paraxial][1] = axis_slope_angle;
           }
           paraxial = FALSE;

           /* Trace marginal ray in C */

           trace_line(3, clear_aperture / 2.0);
           od_cline = object_distance;

           /* Trace marginal ray in F */

           trace_line(6, clear_aperture / 2.0);
           od_fline = object_distance;

           aberr_lspher = od_sa[1][0] - od_sa[0][0];
           aberr_osc = 1.0 - (od_sa[1][0] * od_sa[1][1]) /
              (sin(od_sa[0][1]) * od_sa[0][0]);
           aberr_lchrom = od_fline - od_cline;
           max_lspher = sin(od_sa[0][1]);

           /* D light */

           max_lspher = 0.0000926 / (max_lspher * max_lspher);
           max_osc = 0.0025;
           max_lchrom = max_lspher;
        }

        /* Now evaluate the accuracy of the results from the last ray trace */

        sprintf(outarr[0], "%15s   %21.11f  %14.11f",
           "Marginal ray", od_sa[0][0], od_sa[0][1]);
        sprintf(outarr[1], "%15s   %21.11f  %14.11f",
           "Paraxial ray", od_sa[1][0], od_sa[1][1]);
        sprintf(outarr[2],
           "Longitudinal spherical aberration:      %16.11f",
           aberr_lspher);
        sprintf(outarr[3],
           "    (Maximum permissible):              %16.11f",
           max_lspher);
        sprintf(outarr[4],
           "Offense against sine condition (coma):  %16.11f",
           aberr_osc);
        sprintf(outarr[5],
           "    (Maximum permissible):              %16.11f",
           max_osc);
        sprintf(outarr[6],
           "Axial chromatic aberration:             %16.11f",
           aberr_lchrom);
        sprintf(outarr[7],
           "    (Maximum permissible):              %16.11f",
           max_lchrom);

        /* Now compare the edited results with the master values from
           reference executions of this program. */

        errors = 0;
        for (i = 0; i < 8; i++) {
           if (strcmp(outarr[i], refarr[i]) != 0) {
              printf("\nError in results on line %d...\n", i + 1);
              printf("Expected:  \"%s\"\n", refarr[i]);
              printf("Received:  \"%s\"\n", outarr[i]);
              printf("(Errors)    ");
              k = strlen(refarr[i]);
              for (j = 0; j < k; j++) {
                 printf("%c", refarr[i][j] == outarr[i][j] ? ' ' : '^');
                 if (refarr[i][j] != outarr[i][j])
                    errors++;
              }
              printf("\n");
           }
        }
        if (errors > 0) {
           printf("\n%d error%s in results.  This is VERY SERIOUS.\n",
              errors, errors > 1 ? "s" : "");
        } else {
           printf("\nNo errors in results.\n");
        }
        
        return 0;
}
