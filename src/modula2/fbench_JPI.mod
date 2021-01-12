(*

John Walker's Floating Point Benchmark, derived from...

    Marinchip Interactive Lens Design System
           John Walker   December 1980

               by John Walker
           http://www.fourmilab.ch/

This program may be used, distributed, and modified
freely as long as the origin information is preserved.

This is a complete optical design raytracing algorithm,
stripped of its user interface and recast into
FreeBASIC.  It not only determines execution speed on an
extremely floating point (including trig function)
intensive real-world application, it checks accuracy on
an algorithm that is exquisitely sensitive to errors.
The performance of this program is typically far more
sensitive to changes in the efficiency of the
trigonometric library routines than the average floating
point program.

Implemented in November 2017 by John Walker.

This program is based upon the Visual Basic 6
implementation written in July 2005 by Jim White,
mathimagics@yahoo.co.uk, which may be found in the
visualbasic/vb6 directory of this distribution. I have
removed Visual Basic-specific code, restructured the
program, and added comments in the interest of
readability.

(* ------------------------------------------------------------ *)

An _almost_ line for line Q&D port to Modula-2 in November 2017 ;
Sources are FBENCH.BAS for FreeBasic by John Walker (see _supra_ for details)
and Jim White's VB6 port;
this Modula-2 port and adaptation culprit is Philippe Guiochon ;
compiler is J.P.I. TopSpeed v3.1 for DOS ;
memory model is small ; target CPU is 486 ; all runtime checks are off.

Code is ugly but original is not pretty either. :)

Same license as original, of course.

(* ------------------------------------------------------------ *)

Results (in seconds) on 3 home-made systems (1,000,000 iterations normalized to 1,000) :

0.11968, Pentium II 233 MHz, FBENCH, Novell DOS 7.15, INTRIG FALSE
0.17823, Pentium II 233 MHz, FBENCH, Novell DOS 7.15, INTRIG TRUE
0.12089, Pentium II 233 MHz, FBENCH, Windows 98SE console, INTRIG FALSE
0.17988, Pentium II 233 MHz, FBENCH, Windows 98SE console, INTRIG TRUE
0.06547, Pentium II 233 MHz, FBENCH32, Windows 98SE console, INTRIG FALSE
0.06542, Pentium II 233 MHz, FBENCH32, Windows 98SE console, INTRIG TRUE

0.01582, Pentium IV 2.67 GHz, FBENCH, Novell DOS 7.15, INTRIG FALSE
0.02093, Pentium IV 2.67 GHz, FBENCH, Novell DOS 7.15, INTRIG TRUE
0.01538, Pentium IV 2.67 GHz, FBENCH, Windows 98SE console, INTRIG FALSE
0.02120, Pentium IV 2.67 GHz, FBENCH, Windows 98SE console, INTRIG TRUE
0.01027, Pentium IV 2.67 GHz, FBENCH32, Windows 98SE console, INTRIG FALSE
0.01027, Pentium IV 2.67 GHz, FBENCH32, Windows 98SE console, INTRIG TRUE
0.01632, Pentium IV 2.67 GHz, FBENCH, Windows XP SP2 console, INTRIG FALSE
0.02208, Pentium IV 2.67 GHz, FBENCH, Windows XP SP2 console, INTRIG TRUE
0.00944, Pentium IV 2.67 GHz, FBENCH32, Windows XP SP2 console, INTRIG FALSE
0.00943, Pentium IV 2.67 GHz, FBENCH32, Windows XP SP2 console, INTRIG TRUE

0.00708, Pentium i3 540 3.07 GHz, FBENCH, Novell DOS 7.15, INTRIG FALSE
0.00867, Pentium i3 540 3.07 GHz, FBENCH, Novell DOS 7.15, INTRIG TRUE
0.00620, Pentium i3 540 3.07 GHz, FBENCH, Windows XP SP3 console, INTRIG FALSE
0.00687, Pentium i3 540 3.07 GHz, FBENCH, Windows XP SP3 console, INTRIG TRUE
0.00414, Pentium i3 540 3.07 GHz, FBENCH32, Windows XP SP3 console, INTRIG FALSE
0.00414, Pentium i3 540 3.07 GHz, FBENCH32, Windows XP SP3 console, INTRIG TRUE

*)

(* ------------------------------------------------------------ *)
(* ------------------------------------------------------------ *)

MODULE FBench;

IMPORT Lib,Str;
FROM MATHLIB IMPORT Sin,Cos,Tan,ASin;
FROM IO IMPORT WrStr,WrLn;

(* ------------------------------------------------------------ *)
(* ------------------------------------------------------------ *)

CONST
    INTRIG = FALSE; (* default is FALSE *)
TYPE
    str80 = ARRAY[0..80-1] OF CHAR;

PROCEDURE lr2str (v:LONGREAL):str80;
CONST
    fix = 11;
    wi  = 3;
VAR
    R:str80;
    ok:BOOLEAN;
    i,p:CARDINAL;
BEGIN
    Str.FixRealToStr(v,fix,R,ok);
    IF ok THEN
        p:=Str.CharPos(R,".");
        IF p # MAX(CARDINAL) THEN
            FOR i:=p+1 TO wi DO Str.Prepend(R," ");END;
        END;
    ELSE
        R:="ERROR";
    END;
    RETURN R;
END lr2str;

PROCEDURE lc2a (v:LONGCARD):str80;
CONST
    wi = 9;
VAR
    R:str80;
    ok:BOOLEAN;
    i:CARDINAL;
BEGIN
    Str.CardToStr(v,R,10,ok);
    IF ok THEN
        FOR i:=Str.Length(R)+1 TO wi DO Str.Prepend(R," ");END;
    ELSE
        R:="ERROR";
    END;
    RETURN R;
END lc2a;

PROCEDURE gethms (  ):LONGCARD;
VAR
    h,m,s,ss:CARDINAL;
BEGIN
    Lib.GetTime(h,m,s,ss); (* accuracy of hundredths of a second is a joke here ! *)
    RETURN LONGCARD( LONGCARD(h)*3600+LONGCARD(m*60+s))*100+LONGCARD(ss);
END gethms;

PROCEDURE dmp (VAR Z:ARRAY OF CHAR;
              count:CARDINAL;val1,val2:LONGREAL);
BEGIN
    Str.Subst(Z,"~",lr2str(val1));
    IF count > 1 THEN Str.Subst(Z,"~",lr2str(val2)); END; (* yes, we could count "~" too *)
    WrStr(Z);WrLn;
END dmp;

(* ------------------------------------------------------------ *)
(* ------------------------------------------------------------ *)

(*%T INTRIG  *)

(* ripped and adapted from 2005 VB6 code written by Jim White mathimagics@yahoo.co.uk *)

CONST
    PI         = 3.1415926535897932;
    piover2    = PI / 2.0;
    twopi      = PI * 2.0;
    piover4    = PI / 4.0;
    fouroverpi = 4.0 / PI;
VAR
    atanc : ARRAY [0..8] OF LONGREAL;

PROCEDURE initIntrig(jw:BOOLEAN);
BEGIN
    (* Coefficients for ATAN evaluation *)
    IF jw THEN
        atanc[0] := 0.0;
        atanc[1] := 0.4636476090008061165;
        atanc[2] := 0.7853981633974483094;
        atanc[3] := 0.98279372324732906714;
        atanc[4] := 1.1071487177940905022;
        atanc[5] := 1.1902899496825317322;
        atanc[6] := 1.2490457723982544262;
        atanc[7] := 1.2924966677897852673;
        atanc[8] := 1.3258176636680324644;
    ELSE
        atanc[0] := 0.0;
        atanc[1] := 0.463647609000806;
        atanc[2] := 0.785398163397448;
        atanc[3] := 0.982793723247329;
        atanc[4] := 0.982793723247329;
        atanc[5] := 1.19028994968253;
        atanc[6] := 1.24904577239825;
        atanc[7] := 1.29249666778979;
        atanc[8] := 1.32581766366803;
    END;
END initIntrig;

PROCEDURE fabs( x : LONGREAL) : LONGREAL;
BEGIN
   IF x < 0.0 THEN x := -x; END;
   RETURN x;
END fabs;

PROCEDURE aInt( x : LONGREAL) : LONGREAL;
VAR
    l : LONGINT;
BEGIN
    (*
    aint(x)   Return integer part of number.  Truncates towards 0

    Note that this routine cannot handle the full floating point
    number range.  This PROCEDURE should be in the machine-dependent
    floating point library!
    *)
    l := LONGINT(x);
    IF (LONGINT(-0.5) <> 0) AND (l < 0) THEN INC(l);END;
    RETURN LONGREAL(l);
END aInt;

PROCEDURE fSin(x : LONGREAL) : LONGREAL;
VAR
    sign : BOOLEAN;
    y,r,z: LONGREAL;
BEGIN
    (* sin(x)    Return sine, x in radians *)

    sign := (x < 0.0);

    IF sign THEN x := -x; END;

    IF x > twopi THEN x := x - (aInt(x / twopi) * twopi); END;

    IF (x > PI) THEN x := x - PI;  sign := NOT( sign);END;

    IF (x > piover2) THEN x := PI - x;END;

    IF (x < piover4) THEN
        y := x * fouroverpi;
        z := y * y;
        r := y * (((((((-2.02253129293E-14 * z + 6.9481520350522E-12) * z - 1.75724741761708E-09) * z + 3.13361688917325E-07) * z - 3.65762041821464E-05) * z + 2.49039457019272E-03) * z - 8.07455121882808E-02) * z + 0.785398163397448)
    ELSE
        y := (piover2 - x) * fouroverpi;
        z := y * y;
        r := ((((((-3.8577620372E-13 * z + 1.1500497024263E-10) * z - 2.461136382637E-08) * z + 3.59086044588582E-06) * z - 3.25991886926688E-04) * z + 1.58543442438154E-02) * z - 0.308425137534042) * z + 1.0;
    END;
    IF sign THEN r:=-r;END;
    RETURN r;
END fSin;

PROCEDURE fCos( x : LONGREAL) : LONGREAL;
BEGIN
    (* cos(x)    Return cosine, x in radians, by identity *)
    IF x < 0.0 THEN x := -x; END;

    IF x > twopi THEN                       (* Do range reduction here to limit *)
        x := x - (aInt(x / twopi) * twopi); (* roundoff on add of PI/2 *)
    END;
    RETURN fSin(x + piover2);
END fCos;

PROCEDURE fTan( x : LONGREAL) : LONGREAL;
BEGIN
    RETURN fSin(x) / fCos(x);
END fTan;

PROCEDURE fSqrt( x : LONGREAL) : LONGREAL;
VAR
    c,cl,y:LONGREAL;
    n:CARDINAL;
BEGIN
    (* sqrt(x)   Return square root.  Initial guess, THEN Newton-  Raphson refinement *)

    IF x = 0.0 THEN RETURN 0.0; END;
    IF x < 0.0 THEN HALT; END; (* "Invalid argument for fSqrt()" *)

    y := (0.154116 + 1.893872 * x) / (1.0 + 1.047988 * x);

    c := (y - x / y) / 2.0;
    cl := 0.0;
    n := 50;
    LOOP (* WHILE c <> cl *)
        IF c = cl THEN EXIT; END;
        y := y - c;
        cl := c;
        c := (y - x / y) / 2.0;
        DEC(n);
        IF n = 0 THEN EXIT; END;
    END;
    RETURN y;
END fSqrt;

PROCEDURE fATan( x : LONGREAL) : LONGREAL;
VAR
    sign:BOOLEAN;
    l: BOOLEAN;
    a,b,z,y:LONGREAL;
BEGIN
    (* atan(x)   Return arctangent in radians, range -pi/2 to pi/2 *)
    sign := (x < 0.0);
    IF sign THEN x := -x; END;
    l := FALSE;

    IF x >= 4.0 THEN
        l := TRUE;
        x := 1.0 / x;
        y := 0.0;
    ELSIF x >= 0.25 THEN
        y := aInt(x / 0.5);
        z := y * 0.5;
        x := (x - z) / (x * z + 1.0);
    ELSE
        y := 0.0;
    END;

    z := x * x;
    b := ((((893025.0 * z + 49116375.0) * z + 425675250.0) * z + 1277025750.0) * z + 1550674125.0) * z + 654729075.0;
    a := (((13852575.0 * z + 216602100.0) * z + 891080190.0) * z + 1332431100.0) * z + 654729075.0;
    a := (a / b) * x + atanc[ CARDINAL(y) ];
    IF l THEN a := piover2 - a; END;
    IF sign THEN a := -a;END;
    RETURN a;
END fATan;

PROCEDURE fATan2(y,x: LONGREAL) : LONGREAL;
VAR
    temp:LONGREAL;
BEGIN
    (* atan2(y,x)   Return arctangent in radians of y/x, range -pi to pi *)

    IF x = 0.0 THEN
        IF y = 0.0 THEN RETURN 0.0;END; (* Special case: atan2(0,0) = 0 *)
        IF (y > 0.0) THEN
            RETURN piover2;
        ELSE
            RETURN -piover2;
        END;
    END;

    temp := fATan(y / x);

    IF x < 0.0 THEN
        IF y >= 0.0 THEN
            temp := temp + PI;
        ELSE
            temp := temp - PI;
        END;
    END;
    RETURN temp;
END fATan2;

PROCEDURE fASin( x : LONGREAL) : LONGREAL;
BEGIN
    (* asin(x)   Return arcsine in radians of x *)
    IF fabs(x) > 1.0 THEN HALT; END; (* "Invalid argument for aSin()" *)

    RETURN fATan2(x, fSqrt(1.0 - x * x));
END fASin;

(*%E *)

(* ------------------------------------------------------------ *)
(* ------------------------------------------------------------ *)

CONST
    max_surfaces = 10;

VAR (* ugly global variables as in original ;-) *)
    current_surfaces    : CARDINAL;
    paraxial            : CARDINAL;
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

    (* 1-based *)
    spectral_line       : ARRAY[1..9] OF LONGREAL;
    surf                : ARRAY[1..max_surfaces],[1..5] OF LONGREAL;

    (* 0-based *)
    od_sa               : ARRAY[0..1],[0..1] OF LONGREAL;
    testcase            : ARRAY[0..3],[0..3] OF LONGREAL;

(* Define Cot in terms of Tan *)

PROCEDURE Cot(x : LONGREAL) : LONGREAL;
BEGIN
(*%F INTRIG  *)
    RETURN 1.0 / Tan(x);
(*%E *)
(*%T INTRIG  *)
    RETURN 1.0 / fTan(x);
(*%E *)
END Cot;

(* ------------------------------------------------------------ *)

(*

transit_surface propagates a ray across one surface in
the design, calculating its refraction.  The ray and
surface are described by the following global variables.

paraxial                Axial incidence of ray

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

transit_surface returns false when the last surface has been
traversed and true if more surfaces remain to be traced.

*)

PROCEDURE transit_surface();
VAR
    iang                 : LONGREAL;      (*   Incidence angle *)
    rang                 : LONGREAL;      (*   Refraction angle *)
    iang_sin             : LONGREAL;      (*   Incidence angle sin *)
    rang_sin             : LONGREAL;      (*   Refraction angle sin *)
    old_axis_slope_angle : LONGREAL;
    sagitta              : LONGREAL;
BEGIN
    IF paraxial=1 THEN (* //FIXME *)
        IF radius_of_curvature <> 0.0 THEN
            IF object_distance = 0.0 THEN
                    axis_slope_angle := 0.0;
                    iang_sin := ray_height / radius_of_curvature;
                ELSE
                    iang_sin := ((object_distance - radius_of_curvature) / radius_of_curvature) * axis_slope_angle;
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

            object_distance := object_distance * (to_index / from_index);
            axis_slope_angle := axis_slope_angle * (from_index / to_index);
            RETURN;
    END;

    (*   Marginal ray *)

    IF radius_of_curvature <> 0.0 THEN
        IF object_distance = 0.0 THEN
            axis_slope_angle := 0.0;
            iang_sin := ray_height / radius_of_curvature;
        ELSE
(*%F INTRIG  *)
            iang_sin := ((object_distance - radius_of_curvature) / radius_of_curvature) * Sin(axis_slope_angle)
(*%E *)
(*%T INTRIG  *)
            iang_sin := ((object_distance - radius_of_curvature) / radius_of_curvature) * fSin(axis_slope_angle)
(*%E *)
        END;
(*%F INTRIG  *)
        iang := ASin(iang_sin);
(*%E  *)
(*%T INTRIG  *)
        iang := fASin(iang_sin);
(*%E  *)
        rang_sin := (from_index / to_index) * iang_sin;
        old_axis_slope_angle := axis_slope_angle;
(*%F INTRIG  *)
        axis_slope_angle := axis_slope_angle + iang - ASin(rang_sin);
        sagitta := Sin((old_axis_slope_angle + iang) / 2.0);
(*%E *)
(*%T INTRIG  *)
        axis_slope_angle := axis_slope_angle + iang - fASin(rang_sin);
        sagitta := fSin((old_axis_slope_angle + iang) / 2.0);
(*%E *)
        sagitta := 2.0 * radius_of_curvature * sagitta * sagitta;
(*%F INTRIG  *)
        object_distance := ((radius_of_curvature * Sin(old_axis_slope_angle + iang)) * Cot(axis_slope_angle)) + sagitta;
(*%E *)
(*%T INTRIG  *)
        object_distance := ((radius_of_curvature * fSin(old_axis_slope_angle + iang)) * Cot(axis_slope_angle)) + sagitta;
(*%E *)
        RETURN;
    END;

(*%F INTRIG  *)
    rang := -ASin((from_index / to_index) * Sin(axis_slope_angle));
    object_distance := object_distance * ((to_index * Cos(-rang)) / (from_index * Cos(axis_slope_angle)));
(*%E  *)
(*%T INTRIG  *)
    rang := -fASin((from_index / to_index) * fSin(axis_slope_angle));
    object_distance := object_distance * ((to_index * fCos(-rang)) / (from_index * fCos(axis_slope_angle)));
(*%E  *)
    axis_slope_angle := -rang;
END transit_surface;

(*
trace_line propagates a ray with wavelength Linex and
an initial ray height of ray_h through the surfaces
and leaves the result of the trace in object_distance
and axis_slope_angle.
*)

PROCEDURE trace_line(Linex : CARDINAL; ray_h : LONGREAL);
VAR
    i : CARDINAL;
BEGIN
    object_distance := 0.0;
    ray_height := ray_h;
    from_index := 1.0;

    FOR i := 1 TO current_surfaces DO
        radius_of_curvature := surf[i, 1];
        to_index := surf[i, 2];
        IF to_index > 1.0 THEN
            to_index := to_index + ((spectral_line[4] - spectral_line[Linex]) / (spectral_line[3] - spectral_line[6])) * ((surf[i, 2] - 1.0) / surf[i, 3]);
        END;
        transit_surface();
        from_index := to_index;
        IF (i < current_surfaces) THEN
            object_distance := object_distance - surf[i, 4];
        END;
    END;
END trace_line;

(* ------------------------------------------------------------ *)
(* ------------------------------------------------------------ *)

CONST
    mp = "    (Maximum permissible):               ~";
    defaultniter = 1000;
    info =  ":: "; (* easily removed with NEWLINE utility *)
VAR
    itercount           : LONGCARD;
    niter               : LONGCARD;
VAR
    i,j,k:CARDINAL;
    od_fline, od_cline : LONGREAL;
    Z:str80;
    lc:LONGCARD;
    ok:BOOLEAN;
    chronostart,chronostop,chronodelta:LONGCARD;
    val1,val2:LONGREAL;
BEGIN
    WrLn;

    niter := defaultniter;
(*%T INTRIG *)
    initIntrig(TRUE);
(*%E *)

    (*   Allow iteration count to be specified on the command line *)
    IF Lib.ParamCount() # 0 THEN
        Lib.ParamStr(Z,1);
        lc:=Str.StrToCard(Z,10,ok);
        IF ( ok AND (lc # 0) ) THEN
            niter := lc;
        ELSE
            WrStr("FBENCH (J.P.I. TopSpeed v3.1 Modula-2 compiler for DOS)");WrLn;
            WrLn;
            WrStr("Syntax : FBENCH [iterations]");WrLn;
            WrLn;
            WrStr("a) This is John Walker's FBENCH program (found at fourmilab.ch Web site)");WrLn;
            WrStr("   quickly ported by PhG to J.P.I. TopSpeed v3.1 Modula-2 compiler for DOS.");WrLn;
            WrStr("b) Default value for iterations is 1000.");WrLn;
            WrStr("c) Please note program should be run for about 300 seconds :");WrLn;
            WrStr("   considering timer accuracy, rather high an iteration count is recommended,");WrLn;
            WrStr("   especially on faster systems.");WrLn;
            WrStr("d) This program was compiled with INTRIG constant set to ");
            IF INTRIG THEN
                WrStr("TRUE");
            ELSE
                WrStr("FALSE");
            END;
            WrStr(".");WrLn;
            WrStr("e) Reference results from a run on Microsoft Quick BASIC on the IBM PC/AT :");WrLn;
            WrLn;
            WrStr("   Marginal ray          47.09479120920   0.04178472683");WrLn;
            WrStr("   Paraxial ray          47.08372160249   0.04177864821");WrLn;
            WrStr("Longitudinal spherical aberration:        -0.01106960671");WrLn;
            WrStr("    (Maximum permissible):                 0.05306749907");WrLn;
            WrStr("Offense against sine condition (coma):     0.00008954761");WrLn;
            WrStr("    (Maximum permissible):                 0.00250000000");WrLn;
            WrStr("Axial chromatic aberration:                0.00448229032");WrLn;
            WrStr("    (Maximum permissible):                 0.05306749907");WrLn;
            Lib.SetReturnCode(0);
            HALT;
        END;
    END;

    chronostart:=gethms();

    (*
    The test case used in this program is the design for a 4
    inch achromatic telescope objective used as the example
    IN Wyld's classic work on ray tracing by hand, given in
    Amateur Telescope Making, Volume 3 (Volume 2 in the 1996
    reprint edition).
    *)

    (*   Surface 1 *)
    testcase[0, 0] := 27.05;      (*   Radius of curvature *)
    testcase[0, 1] := 1.5137;     (*   Index of Refraction *)
    testcase[0, 2] := 63.6;       (*   Dispersion          *)
    testcase[0, 3] := 0.52;       (*   Edge Thickness      *)

    (*   Surface 2 *)
    testcase[1, 0] := -16.68;     (*   Radius of curvature *)
    testcase[1, 1] := 1.0;        (*   Index of Refraction *)
    testcase[1, 2] := 0.0;        (*   Dispersion          *)
    testcase[1, 3] := 0.138;      (*   Edge Thickness      *)

    (*   Surface 3 *)
    testcase[2, 0] := -16.68;     (*   Radius of curvature *)
    testcase[2, 1] := 1.6164;     (*   Index of Refraction *)
    testcase[2, 2] := 36.7;       (*   Dispersion          *)
    testcase[2, 3] := 0.38;       (*   Edge Thickness      *)

    (*   Surface 4 *)
    testcase[3, 0] := -78.1;      (*   Radius of curvature *)
    testcase[3, 1] := 1.0;        (*   Index of Refraction *)
    testcase[3, 2] := 0.0;        (*   Dispersion          *)
    testcase[3, 3] := 0.0;        (*   Edge Thickness      *)

    (*
    Wavelengths of standard spectral lines in Angstroms
    (Not all are used in this program)
    *)

    spectral_line[1] := 7621.0;    (*   A  *)
    spectral_line[2] := 6869.955;  (*   B  *)
    spectral_line[3] := 6562.816;  (*   C  *)
    spectral_line[4] := 5895.944;  (*   D  *)
    spectral_line[5] := 5269.557;  (*   E  *)
    spectral_line[6] := 4861.344;  (*   F  *)
    spectral_line[7] := 4340.477;  (*   G' *)
    spectral_line[8] := 3968.494;  (*   H  *)

    (* Load test case into working ARRAY *)

    clear_aperture := 4.0;
    current_surfaces := 4;

    FOR i := 0 TO current_surfaces - 1 DO
        FOR j := 0 TO 3 DO
            surf[i + 1, j + 1] := testcase[i, j];
        END;
    END;

    (* Perform ray trace the specified number of times. *)

    FOR itercount := 1 TO niter DO

        FOR paraxial := 0 TO 1 DO

            (* Do main trace in D light *)

            trace_line(4, clear_aperture / 2.0);
            od_sa[paraxial, 0] := object_distance;
            od_sa[paraxial, 1] := axis_slope_angle;
        END;

        paraxial := 0;

        (* Trace marginal ray in C *)

        trace_line(3, clear_aperture / 2.0);
        od_cline := object_distance;

        (*  Trace marginal ray in F *)

        trace_line(6, clear_aperture / 2.0);
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
(*%F INTRIG  *)
        aberr_osc := 1.0 - (od_sa[1, 0] * od_sa[1, 1]) / (Sin(od_sa[0, 1]) * od_sa[0, 0]);
(*%E  *)
(*%T INTRIG  *)
        aberr_osc := 1.0 - (od_sa[1, 0] * od_sa[1, 1]) / (fSin(od_sa[0, 1]) * od_sa[0, 0]);
(*%E  *)
        (*
        The axial chromatic aberration is the distance between
        where marginal rays in the C and F lines come to focus.
        *)

        aberr_lchrom := od_fline - od_cline;
(*%F INTRIG  *)
        max_lspher := Sin(od_sa[0, 1]);
(*%E *)
(*%T INTRIG  *)
        max_lspher := fSin(od_sa[0, 1]);
(*%E *)

        (*
        Compute maximum acceptable values for each aberration

        Maximum longitudinal spherical aberration, which is
        also the maximum for axial chromatic aberration.  This
        is computed for the D line.
        *)

        max_lspher := 0.0000926 / (max_lspher * max_lspher);
        max_osc := 0.0025;
        max_lchrom := max_lspher;     (* Same criterion *)

    END;

    chronostop:=gethms();

    (*   Print the result in the standard format *)

    FOR i:=1 TO 8 DO
        k:=1;
        val2:=0.0;
        CASE i OF
        | 1 : Z:="   Marginal ray         ~ ~";
              val1:=od_sa[0, 0];
              val2:=od_sa[0, 1];
              k:=2;
        | 2 : Z:="   Paraxial ray         ~ ~";
              val1:=od_sa[1, 0];
              val2:=od_sa[1, 1];
              k:=2;
        | 3 : Z:="Longitudinal spherical aberration:       ~";
              val1:=aberr_lspher;
        | 4 : Z:=mp;
              val1:=max_lspher;
        | 5 : Z:="Offense against sine condition (coma):   ~";
              val1:=aberr_osc;
        | 6 : Z:=mp;
              val1:=max_osc;
        | 7 : Z:="Axial chromatic aberration:              ~";
              val1:=aberr_lchrom;
        | 8 : Z:=mp;
              val1:=max_lchrom;
        ELSE
            ; (* keep brain-damaged stupid ISO newer standard happy :-( *)
        END;
        dmp(Z,k,val1,val2);
    END;

    IF chronostart > chronostop THEN (* handing ONE midnight should do *)
        INC(chronostop,24*3600*100);
    END;
    chronodelta:=chronostop-chronostart;

    WrLn;

    Z:=info+"Time taken for ~ iterations : ~ seconds";
    Str.Subst(Z,"~",lc2a(niter));
    Str.Subst(Z,"~",lr2str(LONGREAL(chronodelta)/100.0));
    WrStr(Z);WrLn;

    IF niter # defaultniter THEN
        val1:=LONGREAL(niter)/LONGREAL(defaultniter);
        Z:=info+"Normalized time for 1000 iterations : ~ seconds";
        dmp(Z,1,(LONGREAL(chronodelta)/100.0)/val1,0.0);
    END;

    WrLn;
    WrStr(info+"Note INTRIG constant is set to ");
    IF INTRIG THEN
        WrStr("TRUE");
    ELSE
        WrStr("FALSE");
    END;
    WrStr(".");WrLn;

    Lib.SetReturnCode(0);
    HALT;
END FBench.
