    %   John Walker's Floating Point Benchmark, derived from...
    %
    %       Marinchip Interactive Lens Design System
    %              John Walker   December 1980
    %
    %                  by John Walker
    %              http://www.fourmilab.ch/
    %
    %   This program may be used, distributed, and modified freely as
    %   long as the origin information is preserved.
    %
    %   This is a complete optical design raytracing algorithm,
    %   stripped of its user interface and recast into Prolog.
    %   It not only determines execution speed on an extremely
    %   floating point (including trig function) intensive
    %   real-world application, it checks accuracy on an algorithm
    %   that is exquisitely sensitive to errors. The performance of
    %   this program is typically far more sensitive to changes in
    %   the efficiency of the trigonometric library routines than the
    %   average floating point program.
    %
    %   Implemented in September 2017 by John Walker.

    %   Wavelengths of standard spectral lines in Angstroms
    %           (Not all are used in this program)

    spectral_line(a, 7621.0  ).         % A
    spectral_line(b, 6869.955).         % B
    spectral_line(c, 6562.816).         % C
    spectral_line(d, 5895.944).         % D
    spectral_line(e, 5269.557).         % E
    spectral_line(f, 4861.344).         % F
    spectral_line(gPrime, 4340.477).    % G'
    spectral_line(h, 3968.494).         % H

    %   The  test case used in this program is the design for a 4 inch
    %   f/12 achromatic telescope objective used as the example in Wyld's
    %   classic work on ray tracing by hand, given in Amateur Telescope
    %   Making, Volume 3 (Volume 2 in the 1996 reprint edition).

    wyldClearAperture(4.0).
    wyldLens( [ [   27.05,  1.5137, 63.6, 0.52  ],
                [ (-16.68), 1.0,     0.0, 0.138 ],
                [ (-16.68), 1.6164, 36.7, 0.38  ],
                [ (-78.1),  1.0,     0.0, 0.0   ]
              ]
            ).

    %   The transit_surface rule traces a ray through a surface
    %   boundary.  Arguments specify the the condition of the
    %   ray entering the boundary.
    %
    %       Radius_of_curvature     Radius of curvature of surface
    %                               being crossed.  If 0, surface is
    %                               plane.
    %
    %       Object_distance         Distance of object focus from
    %                               lens vertex.  If 0, incoming
    %                               rays are parallel and
    %                               the following must be specified:
    %
    %       Ray_height              Height of ray from axis.  Only
    %                               relevant if object_distance =:= 0
    %
    %       Axis_slope_angle        Angle incoming ray makes with axis
    %                               at intercept
    %
    %       From_index              Refractive index of medium being left
    %
    %       To_index                Refractive index of medium being
    %                               entered.
    %
    %   The state of the ray as it exits the boundary is returned in
    %   the following arguments.
    %
    %       OD                          Distance from vertex to object focus
    %                                   after refraction.
    %
    %       RH                          Height of ray from axis.
    %
    %       SA                          Angle incoming ray makes with axis
    %                                   at intercept after refraction.

    %   We define the four cases for paraxial and marginal rays
    %   for flat and curved surfaces.

    %   Paraxial ray with flat radius of curvature

    transit_surface(paraxial,
                    0.0,
                    Object_distance,
                    Ray_height,
                    Axis_slope_angle,
                    From_index,
                    To_index,
                    OD,
                    RH,
                    SA) :-
        OD is Object_distance * (To_index / From_index),
        RH is Ray_height,
        SA is Axis_slope_angle * (From_index / To_index)
    .

    %   Paraxial ray with curved surface

    transit_surface(paraxial,
                    Radius_of_curvature,
                    Object_distance,
                    Ray_height,
                    Axis_slope_angle,
                    From_index,
                    To_index,
                    OD,
                    RH,
                    SA) :-

        ( Object_distance =:= 0 ->
            Axis_slope_angleP is 0,
            Iang_sin is Ray_height / Radius_of_curvature
        ;
            Axis_slope_angleP is Axis_slope_angle,
            Iang_sin is ((Object_distance - Radius_of_curvature) /
                         Radius_of_curvature) * Axis_slope_angle
        ),

        Rang_sin is (From_index / To_index) * Iang_sin,
        Axis_slope_anglePP is (Axis_slope_angleP + Iang_sin) - Rang_sin,

        ( Object_distance =\= 0 ->
            Ray_heightP = Object_distance * Axis_slope_angleP
        ;
            Ray_heightP = Ray_height
        ),

        OD is Ray_heightP / Axis_slope_anglePP,
        RH is Ray_heightP,
        SA is Axis_slope_anglePP
    .

    %   Marginal ray with flat radius of curvature

    transit_surface(marginal,
                    0.0,
                    Object_distance,
                    Ray_height,
                    Axis_slope_angle,
                    From_index,
                    To_index,
                    OD,
                    RH,
                    SA) :-

        Rang is (-(asin((From_index / To_index))) *
                    sin(Axis_slope_angle)),
        OD is Object_distance * ((To_index * cos(-Rang)) /
            (From_index * cos(Axis_slope_angle))),
        RH is Ray_height,
        SA is -Rang
    .

    %   Marginal ray with curved surface

    transit_surface(marginal,
                    Radius_of_curvature,
                    Object_distance,
                    Ray_height,
                    Axis_slope_angle,
                    From_index,
                    To_index,
                    OD,
                    RH,
                    SA) :-

        ( Object_distance =:= 0 ->
            Axis_slope_angleP is 0,
            Iang_sin is Ray_height / Radius_of_curvature
        ;
            Axis_slope_angleP is Axis_slope_angle,
            Iang_sin is ((Object_distance - Radius_of_curvature) /
                         Radius_of_curvature) * sin(Axis_slope_angle)
        ),
        Iang is asin(Iang_sin),
        Rang_sin is (From_index / To_index) * Iang_sin,
        Axis_slope_anglePP is Axis_slope_angleP + Iang - asin(Rang_sin),
        Sagitta is sin((Axis_slope_angle + Iang) / 2.0),
        SagittaP is 2.0 * Radius_of_curvature * Sagitta * Sagitta,
        OD is ((Radius_of_curvature * sin(Axis_slope_angle + Iang)) *
           (1 / tan(Axis_slope_anglePP))) + SagittaP,
        ( Object_distance =:= 0 ->
            RH is Ray_height
        ;
            RH is Object_distance * Axis_slope_angleP
        ),
        SA is Axis_slope_anglePP
    .

    %   Perform ray trace for the design for a specific
    %   spectral line and ray height.  Arguments specify the
    %   desired spectral line, axial incidence, and ray height.
    %   We apply transit_surface to walk the ray through the
    %   surfaces, and then return the object distance and axis
    %   slope angle after traversing all surfaces.

    trace_line([], Line, Incidence, FromIndex, RayHeight, OD_i, SA_i, OD_o, SA_o) :-
        OD_o is OD_i,
        SA_o is SA_i
    .

    trace_line([[Radius_of_curvature, Index_of_refraction,
                Dispersion, Edge_Thickness] | Rest],
                Line, Incidence, FromIndex, RayHeight, OD_i, SA_i, OD_o, SA_o) :-

        spectral_line(Line, Wavelength),
        ( Index_of_refraction > 1 ->
            spectral_line(d, D_line),
            spectral_line(c, C_line),
            spectral_line(f, F_line),
            To_index is Index_of_refraction + ((D_line -
                Wavelength) / (C_line - F_line)) *
                ((Index_of_refraction - 1) / Dispersion)
        ;
            To_index is Index_of_refraction
        ),
        transit_surface(Incidence,
                        Radius_of_curvature,
                        OD_i,
                        RayHeight,
                        SA_i,
                        FromIndex,
                        To_index,
                        OD,
                        RH,
                        SA),
        ODP is OD - Edge_Thickness,
        trace_line(Rest, Line, Incidence, To_index, RH, ODP, SA, OD_o, SA_o)
    .

    %   Perform a ray trace on a given design with a specified clear
    %   aperture and return an evaluation which includes the results
    %   for the D line and calculation of spherical aberration,
    %   coma, and chromatic aberration, along with the conventional
    %   acceptable upper bounds for these quantities.

    evaluate_design(D_marginal_OD, D_marginal_SA,
                    D_paraxial_OD, D_paraxial_SA,
                    Aberr_lspher, Aberr_osc, Aberr_lchrom,
                    Max_lspher, Max_osc, Max_lchrom) :-

        wyldLens(Design),
        wyldClearAperture(Aperture),
        RayHeight is Aperture / 2,

        %   Trace the rays upon which we base the evaluation

        trace_line(Design, d, marginal, 1.0, RayHeight, 0.0, 0.0,
            D_marginal_OD, D_marginal_SA),
        trace_line(Design, d, paraxial, 1.0, RayHeight, 0.0, 0.0,
            D_paraxial_OD, D_paraxial_SA),
        trace_line(Design, c, marginal, 1.0, RayHeight, 0.0, 0.0,
            C_marginal_OD, C_marginal_SA),
        trace_line(Design, f, marginal, 1.0, RayHeight, 0.0, 0.0,
            F_marginal_OD, F_marginal_SA),

        %   Compute aberrations of the design

        %   The longitudinal spherical aberration is just the
        %   difference between where the D line comes to focus
        %   for paraxial and marginal rays.

        Aberr_lspher is D_paraxial_OD - D_marginal_OD,

        %   The offense against the sine condition is a measure
        %   of the degree of coma in the design.  We compute it
        %   as the lateral distance in the focal plane between
        %   where a paraxial ray and marginal ray in the D line
        %   come to focus.

        Aberr_osc is 1 - ((D_paraxial_OD) * (D_paraxial_SA)) /
                        (sin(D_marginal_SA) * (D_marginal_OD)),

        %  The axial chromatic aberration is the distance between
        %  where marginal rays in the C and F lines come to focus.

        Aberr_lchrom is F_marginal_OD - C_marginal_OD,

        %   Maximum longitudinal spherical aberration, which is
        %   also the maximum for axial chromatic aberration.
        %   This is computed for the D line.

        Max_lspher_s is sin(D_marginal_SA),
        Max_lspher is 0.0000926 / (Max_lspher_s * Max_lspher_s),
        Max_osc is 0.0025,
        Max_lchrom is Max_lspher
    .

    %   Take the output from evaluate_design and return a list
    %   of strings containing a primate-readable version of the
    %   evaluation. For the purposes of the benchmark, this
    %   function is not timed; it serves only to create textual
    %   results which can be compared against those from the
    %   reference implementation.

    evaluation_report(Dm_OD, Dm_SA, Dp_OD, Dp_SA,
        Aberr_lspher, Aberr_osc, Aberr_lchrom,
        Max_lspher, Max_osc, Max_lchrom,
        [L1, L2, L3, L4, L5, L6, L7, L8]) :-
        format(atom(L1), '   Marginal ray          ~11f   ~11f',
            [Dm_OD, Dm_SA]),
        format(atom(L2), '   Paraxial ray          ~11f   ~11f',
            [Dp_OD, Dp_SA]),
        format(atom(L3), 'Longitudinal spherical aberration:        ~11f',
            [Aberr_lspher]),
        format(atom(L4), '    (Maximum permissible):                 ~11f',
            [Max_lspher]),
        format(atom(L5), 'Offense against sine condition (coma):     ~11f',
            [Aberr_osc]),
        format(atom(L6), '    (Maximum permissible):                 ~11f',
            [Max_osc]),
        format(atom(L7), 'Axial chromatic aberration:                ~11f',
            [Aberr_lchrom]),
        format(atom(L8), '    (Maximum permissible):                 ~11f',
            [Max_lchrom])
    .

    %   This version of evaluation_report prints the report on
    %   the terminal, using the standard format/2 predicate.

    evaluation_report_ISO(Dm_OD, Dm_SA, Dp_OD, Dp_SA,
        Aberr_lspher, Aberr_osc, Aberr_lchrom,
        Max_lspher, Max_osc, Max_lchrom) :-
        format('   Marginal ray          ~11f   ~11f~n',
            [Dm_OD, Dm_SA]),
        format('   Paraxial ray          ~11f   ~11f~n',
            [Dp_OD, Dp_SA]),
        format('Longitudinal spherical aberration:        ~11f~n',
            [Aberr_lspher]),
        format('    (Maximum permissible):                 ~11f~n',
            [Max_lspher]),
        format('Offense against sine condition (coma):     ~11f~n',
            [Aberr_osc]),
        format('    (Maximum permissible):                 ~11f~n',
            [Max_osc]),
        format('Axial chromatic aberration:                ~11f~n',
            [Aberr_lchrom]),
        format('    (Maximum permissible):                 ~11f~n',
            [Max_lchrom])
    .

    %   Print an evaluation report on the console (used only for debugging)

    print_evaluation([]).

    print_evaluation([Head | Tail]) :-
        write(Head), nl,
        print_evaluation(Tail)
    .

    %   The validateResults rule compares the primate-readable
    %   eReport from evaluation_report with the archival results
    %   from the reference implementation (which all language
    %   implementations must reproduce character-by-character
    %   [apart from trivia such as end of line conventions and
    %   trailing white space]).  It returns a Boolean indicating
    %   whether the results matched.

    validate_results([], [], _, Errors) :-
        ( Errors > 0 ->
            format('~d error(s) detected in results.  This is VERY SERIOUS.~n', Errors),
            false
        ;
            true
        )
    .

    validate_results([ExpH | ExpT], [RcvH | RcvT], Line, Errors) :-
        ( ExpH \= RcvH ->
            format('Validation failed on line ~d:~n', Line),
            format('  Expected: ''(~w)''~n', ExpH),
            format('  Received: ''(~w)''~n', RcvH),
            ErrorsN is Errors + 1
        ;
            ErrorsN is Errors
        ),
        validate_results(ExpT, RcvT, Line + 1, ErrorsN)
    .

    validate_results(Report) :-
        %   Reference results.  These happen to be derived from a
        %   run on Microsoft Quick BASIC on the IBM PC/AT.
        ExpectedResults = [
            '   Marginal ray          47.09479120920   0.04178472683',
            '   Paraxial ray          47.08372160249   0.04177864821',
            'Longitudinal spherical aberration:        -0.01106960671',
            '    (Maximum permissible):                 0.05306749907',
            'Offense against sine condition (coma):     0.00008954761',
            '    (Maximum permissible):                 0.00250000000',
            'Axial chromatic aberration:                0.00448229032',
            '    (Maximum permissible):                 0.05306749907'
        ],
        validate_results(ExpectedResults, Report, 1, 0)
    .

    %   Run the benchmark for the specified number of iterations

    %   This version performs its own internal validation of the
    %   results, which requires the format/3 predicate provided
    %   by SWI-Prolog.  If your Prolog system does not support
    %   this predicate, use run_benchmark_ISO below which uses
    %   the standard format/2 predicate to write the results of
    %   the last iteration to the terminal where they can be
    %   manually checked for accuracy.

    run_benchmark(0) :- !.

    run_benchmark(N) :-
        evaluate_design(Dm_OD, Dm_SA, Dp_OD, Dp_SA,
            Aberr_lspher, Aberr_osc, Aberr_lchrom,
            Max_lspher, Max_osc, Max_lchrom),
        Y is N - 1,

        %   If this is the final iteration, create an
        %   evaluation_report and pass it to validate_results
        %   to verify the values computed were as expected.

        (Y > 0 ->
            true
        ;
            evaluation_report(Dm_OD, Dm_SA, Dp_OD, Dp_SA,
                Aberr_lspher, Aberr_osc, Aberr_lchrom,
                Max_lspher, Max_osc, Max_lchrom, Report),
            validate_results(Report)
        ),
        run_benchmark(Y)
    .

   %    This is the version of the benchmark which doesn't
   %    require SWI-Prolog's format/3.  It prints the results on
   %    the terminal where it can be manually validated.

    run_benchmark_ISO(0) :- !.

    run_benchmark_ISO(N) :-
        evaluate_design(Dm_OD, Dm_SA, Dp_OD, Dp_SA,
            Aberr_lspher, Aberr_osc, Aberr_lchrom,
            Max_lspher, Max_osc, Max_lchrom),
        Y is N - 1,

        %   If this is the final iteration, print an
        %   evaluation_report.

        (Y > 0 ->
            true
        ;
            evaluation_report_ISO(Dm_OD, Dm_SA, Dp_OD, Dp_SA,
                Aberr_lspher, Aberr_osc, Aberr_lchrom,
                Max_lspher, Max_osc, Max_lchrom)
        ),
        run_benchmark_ISO(Y)
    .
