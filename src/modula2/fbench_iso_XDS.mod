(*

    John Walker's Floating Point Benchmark, derived from...

         Marinchip Interactive Lens Design System
                John Walker   December 1980

                    by John Walker
                http://www.fourmilab.ch/

    This program may be used, distributed, and modified freely
    as long as the origin information is preserved.

    This is a complete optical design raytracing algorithm,
    stripped of its user interface and recast into Modula-2.  It
    not only determines execution speed on an extremely floating
    point (including trig function) intensive real-world
    application, it checks accuracy on an algorithm that is
    exquisitely sensitive to errors. The performance of this
    program is typically far more sensitive to changes in the
    efficiency of the trigonometric library routines than the
    average floating point program.

    This program, in ISO Modula-2, is based upon the original
    Modula-2 port for J.P.I. TopSpeed v3.1 by Philippe Guiochon
    (fbench_JPI.mod in this directory) which, in turn, was based
    upon the FreeBASIC version of the benchmark (directory
    freebasic in the archive), which itself was based upon Jim
    White's (mathimagics@yahoo.co.uk) Visual Basic version,
    which may be found in the visualbasic/vb6 directory.  I
    have removed the INTRIG code from the benchmark, as we're
    no longer interested in comparing re-implementation of
    the mathematical function library.

    This version of the program was adapted by Philippe Guiochon
    from the GNU Modula-2 implementation (fbench.mod in this
    directory) to be compatible with the Excelsior XDS Modula-2
    compiler.  The XDS compiler defines the REAL type to be
    32-bit single-precision floating point, which is
    insufficient for the accuracy needs of this algorithm, and
    requires variables to be declared as LONGREAL to use 64-bit
    double precision (which then requires importing the
    trigonometric functions from the LongMath module).  In
    addition, this compiler does not permit enumeration types to
    be used as array indices or bounds in FOR statements and
    requires them to be explicitly converted to CARDINAL with
    the ORD() function.  Using the ORD() function is ugly but
    compatible with GNU Modula-2, but declaring variables as
    LONGREAL causes that compiler to use the Intel 80-bit
    floating point type (GCC's "long double") instead of the
    standard 64-bit data type it uses for REAL.  Since this may
    affect the timing (and still other compilers may define REAL
    and LONGREAL in yet other ways), I have kept this version
    separate and used the GNU Modula-2 fbench.mod for the
    cross-language timing tests.

*)

MODULE fbench_iso_XDS;

FROM ProgramArgs IMPORT ArgChan, IsArgPresent, NextArg;
FROM WholeIO IMPORT ReadCard;

FROM LongMath IMPORT sin, cos, tan, arcsin;
FROM STextIO IMPORT WriteString, WriteLn;
FROM SLongIO IMPORT WriteFixed;

CONST
    max_surfaces = 10;

TYPE
    Axial_Incidence = (marginal, paraxial);
    Spectral_Line = (A_line, B_line, C_line, D_line, E_line,
                     F_line, Gprime_line, H_line);

VAR
    current_surfaces    : CARDINAL;
    clear_aperture      : LONGREAL;
    aberr_lspher        : LONGREAL;
    aberr_osc           : LONGREAL;
    aberr_lchrom        : LONGREAL;
    max_lspher          : LONGREAL;
    max_osc             : LONGREAL;
    max_lchrom          : LONGREAL;
    radius_of_curvature : LONGREAL;
    object_distance     : LONGREAL;
    ray_height          : LONGREAL;
    axis_slope_angle    : LONGREAL;
    from_index          : LONGREAL;
    to_index            : LONGREAL;

    spectral_line       : ARRAY[0..8] OF LONGREAL;
    od_sa               : ARRAY[0..1], [0..1] OF LONGREAL;
    WyldLens            : ARRAY[0..3], [0..3] OF LONGREAL;
    surf                : ARRAY[0..max_surfaces - 1], [0..4] OF LONGREAL;

(*
    transit_surface propagates a ray across one surface in
    the design, calculating its refraction.  The ray and
    surface are described by the following global variables.

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

    The result of a call on transit_surface updates the
    following variables modified to reflect the geometry of
    the ray as it exits the surface.

    object_distance         Distance from vertex to object focus
                            after refraction.

    ray_height              Height of ray from axis.

    axis_slope_angle        Angle incoming ray makes with axis
                            at intercept after refraction.
*)

PROCEDURE transit_surface(axial_incidence: Axial_Incidence);
VAR
    iang                 : LONGREAL;      (*  Incidence angle  *)
    rang                 : LONGREAL;      (*  Refraction angle  *)
    iang_sin             : LONGREAL;      (*  Incidence angle sin  *)
    rang_sin             : LONGREAL;      (*  Refraction angle sin  *)
    old_axis_slope_angle : LONGREAL;
    sagitta              : LONGREAL;
BEGIN
    IF axial_incidence = paraxial THEN

        (*  Paraxial ray  *)

        IF radius_of_curvature <> 0.0 THEN
            (*  Curved surface  *)
            IF object_distance = 0.0 THEN
                    axis_slope_angle := 0.0;
                    iang_sin := ray_height / radius_of_curvature;
                ELSE
                    iang_sin := ((object_distance - radius_of_curvature) /
                        radius_of_curvature) * axis_slope_angle;
                END;
                rang_sin := (from_index / to_index) * iang_sin;
                old_axis_slope_angle := axis_slope_angle;
                axis_slope_angle := axis_slope_angle + iang_sin - rang_sin;
                IF object_distance <> 0.0 THEN
                    ray_height := object_distance * old_axis_slope_angle;
                END;
                object_distance := ray_height / axis_slope_angle;
                RETURN;
            END;

            (*  Flat surface  *)
            object_distance := object_distance * (to_index / from_index);
            axis_slope_angle := axis_slope_angle * (from_index / to_index);
            RETURN;
    END;

    (*  Marginal ray  *)

    IF radius_of_curvature <> 0.0 THEN
        (*  Curved surface  *)
        IF object_distance = 0.0 THEN
            axis_slope_angle := 0.0;
            iang_sin := ray_height / radius_of_curvature;
        ELSE
            iang_sin := ((object_distance - radius_of_curvature) /
                radius_of_curvature) * sin(axis_slope_angle);
        END;
        iang := arcsin(iang_sin);
        rang_sin := (from_index / to_index) * iang_sin;
        old_axis_slope_angle := axis_slope_angle;
        axis_slope_angle := axis_slope_angle + iang - arcsin(rang_sin);
        sagitta := sin((old_axis_slope_angle + iang) / 2.0);
        sagitta := 2.0 * radius_of_curvature * sagitta * sagitta;
        object_distance := ((radius_of_curvature *
            sin(old_axis_slope_angle + iang)) *
            (1.0 / tan(axis_slope_angle))) + sagitta;
        RETURN;
    END;

    (*  Flat surface  *)
    rang := -arcsin((from_index / to_index) * sin(axis_slope_angle));
    object_distance := object_distance * ((to_index * cos(-rang)) /
        (from_index * cos(axis_slope_angle)));
    axis_slope_angle := -rang;
END transit_surface;

(*
    trace_line propagates a ray with wavelength line and
    an initial ray height of ray_h through the surfaces
    and leaves the result of the trace in object_distance
    and axis_slope_angle.
*)

PROCEDURE trace_line(incidence: Axial_Incidence;
                     line : Spectral_Line; ray_h : LONGREAL);
VAR
    i : CARDINAL;
BEGIN
    object_distance := 0.0;
    ray_height := ray_h;
    from_index := 1.0;

    FOR i := 0 TO current_surfaces - 1 DO
        radius_of_curvature := surf[i, 0];
        to_index := surf[i, 1];
        (*  Adjust effective index of refraction according to
            dispersion and the wavelength of the ray being
            traced.  *)
        IF to_index > 1.0 THEN
            to_index := to_index +
                ((spectral_line[ORD(D_line)] - spectral_line[ORD(line)]) /
                                (spectral_line[ORD(C_line)] - spectral_line[ORD(F_line)])) *
                                ((surf[i, 1] - 1.0) / surf[i, 2]);
        END;
        transit_surface(incidence);
        from_index := to_index;
        IF (i < current_surfaces) THEN
            object_distance := object_distance - surf[i, 3];
        END;
    END;
END trace_line;


CONST
    mp = "    (Maximum permissible):              ";
    defaultniter = 1000000;

VAR
    itercount           : CARDINAL;
    niter               : CARDINAL;

VAR
    i, j: CARDINAL;
    od_fline, od_cline: LONGREAL;
    incidence: Axial_Incidence;

BEGIN

    niter := 0;

    (*   Allow iteration count to be specified on the command line *)

    NextArg();
    IF IsArgPresent() THEN
        ReadCard(ArgChan(), niter);
        IF niter = 0 THEN
            WriteString("Invalid iteration count on command line.");
            WriteLn;
            RETURN;
        END;
    END;
    IF niter < 1 THEN
        niter := defaultniter;
    END;

    (*
        The test case used in this program is the design for a 4
        inch achromatic telescope objective used as the example
        in Wyld's classic work on ray tracing by hand, given in
        Amateur Telescope Making, Volume 3 (Volume 2 in the 1996
        reprint edition).
    *)

    (*   Surface 1   *)
    WyldLens[0, 0] := 27.05;      (*   Radius of curvature   *)
    WyldLens[0, 1] := 1.5137;     (*   Index of Refraction   *)
    WyldLens[0, 2] := 63.6;       (*   Dispersion            *)
    WyldLens[0, 3] := 0.52;       (*   Edge Thickness        *)

    (*   Surface 2   *)
    WyldLens[1, 0] := -16.68;     (*   Radius of curvature   *)
    WyldLens[1, 1] := 1.0;        (*   Index of Refraction   *)
    WyldLens[1, 2] := 0.0;        (*   Dispersion            *)
    WyldLens[1, 3] := 0.138;      (*   Edge Thickness        *)

    (*   Surface 3   *)
    WyldLens[2, 0] := -16.68;     (*   Radius of curvature   *)
    WyldLens[2, 1] := 1.6164;     (*   Index of Refraction   *)
    WyldLens[2, 2] := 36.7;       (*   Dispersion            *)
    WyldLens[2, 3] := 0.38;       (*   Edge Thickness        *)

    (*   Surface 4   *)
    WyldLens[3, 0] := -78.1;      (*   Radius of curvature   *)
    WyldLens[3, 1] := 1.0;        (*   Index of Refraction   *)
    WyldLens[3, 2] := 0.0;        (*   Dispersion            *)
    WyldLens[3, 3] := 0.0;        (*   Edge Thickness        *)

    (*
        Wavelengths of standard spectral lines in Angstroms
                 (Not all are used in this program)
    *)

    spectral_line[ORD(A_line)]      := 7621.0;
    spectral_line[ORD(B_line)]      := 6869.955;
    spectral_line[ORD(C_line)]      := 6562.816;
    spectral_line[ORD(D_line)]      := 5895.944;
    spectral_line[ORD(E_line)]      := 5269.557;
    spectral_line[ORD(F_line)]      := 4861.344;
    spectral_line[ORD(Gprime_line)] := 4340.477;
    spectral_line[ORD(H_line)]      := 3968.494;

    (* Load test case into working ARRAY surf *)

    clear_aperture := 4.0;
    current_surfaces := 4;

    FOR i := 0 TO current_surfaces - 1 DO
        FOR j := 0 TO 3 DO
            surf[i, j] := WyldLens[i, j];
        END;
    END;

    (* Perform ray trace the specified number of times. *)

    FOR itercount := 1 TO niter DO

        FOR incidence := marginal TO paraxial DO

            (* Do main trace in D light *)

            trace_line(incidence, D_line, clear_aperture / 2.0);
            od_sa[ORD(incidence), 0] := object_distance;
            od_sa[ORD(incidence), 1] := axis_slope_angle;
        END;

        (* Trace marginal ray in C *)

        trace_line(marginal, C_line, clear_aperture / 2.0);
        od_cline := object_distance;

        (* Trace marginal ray in F *)

        trace_line(marginal, F_line, clear_aperture / 2.0);
        od_fline := object_distance;

        (*
            Compute aberrations of the design

            The longitudinal spherical aberration is just the
            difference between where the D line comes to focus
            for paraxial and marginal rays.
        *)
        aberr_lspher := od_sa[1, 0] - od_sa[0, 0];

        (*
            The offense against the sine condition is a measure
            of the degree of coma in the design.  We compute it
            as the lateral distance in the focal plane between
            where a paraxial ray and marginal ray in the D line
            come to focus.
        *)
        aberr_osc := 1.0 - (od_sa[1, 0] * od_sa[1, 1]) / (sin(od_sa[0, 1]) * od_sa[0, 0]);

        (*
            The axial chromatic aberration is the distance between
            where marginal rays in the C and F lines come to focus.
        *)
        aberr_lchrom := od_fline - od_cline;
        max_lspher := sin(od_sa[0, 1]);

        (*
            Compute maximum acceptable values for each aberration

            Maximum longitudinal spherical aberration, which is
            also the maximum for axial chromatic aberration.  This
            is computed for the D line.
        *)
        max_lspher := 0.0000926 / (max_lspher * max_lspher);
        max_osc := 0.0025;
        max_lchrom := max_lspher;   (* Same criterion *)

    END;

    (* Print the results of the trace in the standard format. *)

    WriteString("   Marginal ray        ");
    WriteFixed(od_sa[0, 0], 11, 16);
    WriteFixed(od_sa[0, 1], 11, 16);
    WriteLn;

    WriteString("   Paraxial ray        ");
    WriteFixed(od_sa[1, 0], 11, 16);
    WriteFixed(od_sa[1, 1], 11, 16);
    WriteLn;

    WriteString("Longitudinal spherical aberration:      ");
    WriteFixed(aberr_lspher, 11, 16);
    WriteLn;
    WriteString(mp);
    WriteFixed(max_lspher, 11, 16);
    WriteLn;

    WriteString("Offense against sine condition (coma):  ");
    WriteFixed(aberr_osc, 11, 16);
    WriteLn;
    WriteString(mp);
    WriteFixed(max_osc, 11, 16);
    WriteLn;

    WriteString("Axial chromatic aberration:             ");
    WriteFixed(aberr_lchrom, 11, 16);
    WriteLn;
    WriteString(mp);
    WriteFixed(max_lchrom, 11, 16);
    WriteLn;

END fbench_iso_XDS.
