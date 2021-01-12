{

        John Walker's Floating Point Benchmark, derived from...

        Marinchip Interactive Lens Design System

                                     John Walker   December 1980

        By John Walker
           http://www.fourmilab.ch/

        This program may be used, distributed, and modified freely as
        long as the origin information is preserved.

        This is a complete optical design raytracing algorithm,
        stripped of its user interface and recast into Pascal. It not
        only determines execution speed on an extremely floating point
        (including trig function) intensive real-world application, it
        allows checking accuracy on an algorithm that is exquisitely
        sensitive to errors.  The performance of this program is
        typically far more sensitive to changes in the efficiency of
        the trigonometric library routines than the average floating
        point program.

        Ported from the Ada language implementation in September 2007
        by John Walker.
}

program fbench(input, output);

    const
        OUTER = 5753;
	INNER = 5753;
        max_surfaces = 10;
        current_surfaces = 4;
        clear_aperture = 4.0;

    type
        Axial_Incidence = ( Marginal_Ray, Paraxial_Ray );

        { Wavelengths of standard spectral lines in Angstroms
                 (Not all are used in this program) }

        Spectral_Line_Name = (
            A_line, B_line, C_line, D_line, E_line, F_line,
            Gprime_line, H_line
        );

        Surface_Property = (
            Curvature_Radius,       { Radius of curvature }
                                    {   (+ if convex to light source }
            Index_Of_Refraction,    { Index of refraction (1 for air space }
            Dispersion,             { Dispersion (Abbe number (V)) }
            Edge_Thickness          { Edge thickness (0 for last surface) }
        );

    var
        paraxial: Axial_Incidence;
        aberr_lspher, aberr_osc, aberr_lchrom,
        max_lspher, max_osc, max_lchrom,
        radius_of_curvature, object_distance, ray_height,
        axis_slope_angle, from_index, to_index,
        od_cline, od_fline: real;

        testcase: array[1..4, Surface_Property] of real;
        spectral_line: array[Spectral_Line_Name] of real;
        s: array[1..max_surfaces, Surface_Property] of real;
        od_sa: array[Axial_Incidence, 0..1] of real;

        i, itercount, thousand: integer;
        sp: Surface_Property;
        p: Axial_Incidence;

    {  The arcsin function defined in terms of arctan and sqrt  }

    function arcsin(x: real) : real;
        begin
            arcsin := arctan(x / sqrt(1 - sqr(x)));
        end;

    {  The tan function defined in terms of sin and cos  }

    function tan(x: real) : real;
        begin
            tan := sin(x) / cos(x);
        end;

    {  The cot function defined in terms of tan  }

    function cot(x: real) : real;
        begin
            cot := 1.0 / tan(x);
        end;

    procedure initialise;
    begin

        {  Wavelengths of standard spectral lines in Angstroms
           Not all are used in this program) }

        spectral_line[A_line] := 7621.0;         { A }
        spectral_line[B_line] := 6869.955;       { B }
        spectral_line[C_line] := 6562.816;       { C }
        spectral_line[D_line] := 5895.944;       { D }
        spectral_line[E_line] := 5269.557;       { E }
        spectral_line[F_line] := 4861.344;       { F }
        spectral_line[Gprime_line] := 4340.477;  { G' }
        spectral_line[H_line] := 3968.494;       { H }

        { The  test case used in this program is the design for a 4
          inch f/12 achromatic telescope objective used as the example
          in Wyld's classic work on ray tracing by hand, given in
          Amateur Telescope Making, Volume 3 (Volume 2 in the 1996
          reprint edition). }

        testcase[1][Curvature_Radius] := 27.05;
        testcase[1][Index_Of_Refraction] := 1.5137;
        testcase[1][Dispersion] := 63.6;
        testcase[1][Edge_Thickness] := 0.52;

        testcase[2][Curvature_Radius] := -16.68;
        testcase[2][Index_Of_Refraction] := 1;
        testcase[2][Dispersion] := 0;
        testcase[2][Edge_Thickness] := 0.138;

        testcase[3][Curvature_Radius] := -16.68;
        testcase[3][Index_Of_Refraction] := 1.6164;
        testcase[3][Dispersion] := 36.7;
        testcase[3][Edge_Thickness] := 0.38;

        testcase[4][Curvature_Radius] := -78.1;
        testcase[4][Index_Of_Refraction] := 1;
        testcase[4][Dispersion] := 0;
        testcase[4][Edge_Thickness] := 0;
    end;


    {     Calculate passage through surface

          If the variable paraxial is Paraxial_Ray, the trace through the
          surface will be done using the paraxial approximations.
          Otherwise, the normal trigonometric trace will be done.

          This subroutine takes the following global inputs:

          radius_of_curvature     Radius of curvature of surface
                                  being crossed.  If 0, surface is
                                  plane.

          object_distance         Distance of object focus from
                                  lens vertex.  If 0, incoming
                                  rays are parallel and
                                  the following must be specified:

          ray_height              Height of ray from axis.  Only
                                  relevant if $object_distance == 0

          axis_slope_angle        Angle incoming ray makes with axis
                                  at intercept

          from_index              Refractive index of medium being left

          to_index                Refractive index of medium being
                                  entered.

          The outputs are the following global variables:

          object_distance         Distance from vertex to object focus
                                  after refraction.

          axis_slope_angle        Angle incoming ray makes with axis
                                  at intercept after refraction. }

        procedure transit_surface;
        var
            iang,                   {  Incidence angle  }
            rang,                   {  Refraction angle  }
            iang_sin,               {  Incidence angle sin  }
            rang_sin,               {  Refraction angle sin  }
            old_axis_slope_angle, sagitta : real;
        begin
            if paraxial = Paraxial_Ray then begin
                if radius_of_curvature <> 0.0 then begin
                    if object_distance = 0.0  then begin
                        axis_slope_angle := 0.0;
                        iang_sin := ray_height / radius_of_curvature;
                    end else begin
                        iang_sin := ((object_distance -
                        radius_of_curvature) / radius_of_curvature) *
                        axis_slope_angle;
                    end;
                    rang_sin := (from_index / to_index) * iang_sin;
                    old_axis_slope_angle := axis_slope_angle;
                    axis_slope_angle := axis_slope_angle +
                        iang_sin - rang_sin;
                    if object_distance <> 0.0 then begin
                        ray_height := object_distance * old_axis_slope_angle;
                    end;
                    object_distance := ray_height / axis_slope_angle;
                end else begin
                    object_distance := object_distance * (to_index / from_index);
                    axis_slope_angle := axis_slope_angle * (from_index / to_index);
                end;
            end else begin
                if radius_of_curvature <> 0.0 then begin
                    if object_distance = 0.0 then begin
                        axis_slope_angle := 0.0;
                        iang_sin := ray_height / radius_of_curvature;
                    end else begin
                        iang_sin := ((object_distance -
                            radius_of_curvature) / radius_of_curvature) *
                            Sin(axis_slope_angle);
                    end;
                    iang := Arcsin(iang_sin);
                    rang_sin := (from_index / to_index) * iang_sin;
                    old_axis_slope_angle := axis_slope_angle;
                    axis_slope_angle := axis_slope_angle +
                        iang - Arcsin(rang_sin);
                    sagitta := Sin((old_axis_slope_angle + iang) / 2.0);
                    sagitta := 2.0 * radius_of_curvature * sagitta * sagitta;
                    object_distance := ((radius_of_curvature *
                        Sin(old_axis_slope_angle + iang)) *
                        Cot(axis_slope_angle)) + sagitta;
                end else begin
                    rang := -Arcsin((from_index / to_index) *
                        Sin(axis_slope_angle));
                    object_distance := object_distance * ((to_index *
                        Cos(-rang)) / (from_index *
                        Cos(axis_slope_angle)));
                    axis_slope_angle := -rang;
                end;
            end;
        end;


        {  Perform ray trace in specific spectral line  }

        procedure trace_line (line : Spectral_Line_Name; ray_h : real);
        var
            i: integer;
        begin
            object_distance := 0.0;
            ray_height := ray_h;
            from_index := 1.0;

            for i := 1 to current_surfaces do begin
                radius_of_curvature := s[i, Curvature_Radius];
                to_index := s[i, Index_Of_Refraction];
                if to_index > 1.0 then begin
                    to_index := to_index + ((spectral_line[D_line] -
                        spectral_line[line]) /
                        (spectral_line[C_line] - spectral_line[F_line])) *
                        ((s[i, Index_Of_Refraction] - 1.0) / s[i, Dispersion]);
                end;
                transit_surface;
                from_index := to_index;
                if i < current_surfaces then begin
                    object_distance := object_distance - s[i, Edge_Thickness];
                end;
            end;
        end;

    {  Main program  }

    begin

        {  Load test case into working array  }

        initialise;
        for i := 1 to current_surfaces do begin
            for sp := Curvature_Radius to Edge_Thickness do begin
                s[i, sp] := testcase[i, sp];
            end;
        end;

        write('Press return to begin benchmark: ');
        readln;

        {  Timing begins here  }

        for itercount := 1 to OUTER do begin
            for thousand := 1 to INNER do begin
                for p := Marginal_Ray to Paraxial_Ray do begin

                    {  Do main trace in D light  }

                    paraxial := p;
                    trace_line(D_line, clear_aperture / 2.0);
                    od_sa[paraxial, 0] := object_distance;
                    od_sa[paraxial, 1] := axis_slope_angle;
                end;
                paraxial := Marginal_Ray;

                {  Trace marginal ray in C  }

                trace_line(C_line, clear_aperture / 2.0);
                od_cline := object_distance;

                {  Trace marginal ray in F  }

                trace_line(F_line, clear_aperture / 2.0);
                od_fline := object_distance;

                aberr_lspher := od_sa[Paraxial_Ray, 0] - od_sa[Marginal_Ray, 0];
                aberr_osc := 1.0 - (od_sa[Paraxial_Ray, 0] * od_sa[Paraxial_Ray, 1]) /
                   (Sin(od_sa[Marginal_Ray, 1]) * od_sa[Marginal_Ray, 0]);
                aberr_lchrom := od_fline - od_cline;
                max_lspher := Sin(od_sa[Marginal_Ray, 1]);

                {  D light  }

                max_lspher := 0.0000926 / (max_lspher * max_lspher);
                max_osc := 0.0025;
                max_lchrom := max_lspher;
            end;
        end;

        {  Timing ends here  }

        write('Stop the timer: ');
        readln;
        writeln;

        writeln('   Marginal ray   ', od_sa[Marginal_Ray, 0]:21:11,
            '  ', od_sa[Marginal_Ray, 1]:14:11);
        writeln('   Paraxial ray   ', od_sa[Paraxial_Ray, 0]:21:11,
            '  ', od_sa[Paraxial_Ray, 1]:14:11);
        writeln('Longitudinal spherical aberration:      ', aberr_lspher:16:11);
        writeln('    (Maximum permissible):              ', max_lspher:16:11);
        writeln('Offense against sine condition (coma):  ', aberr_osc:16:11);
        writeln('    (Maximum permissible):              ', max_osc:16:11);
        writeln('Axial chromatic aberration:             ', aberr_lchrom:16:11);
        writeln('    (Maximum permissible):              ', max_lchrom:16:11);

    end.
