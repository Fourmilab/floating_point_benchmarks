
    '   John Walker's Floating Point Benchmark, derived from...
    '
    '       Marinchip Interactive Lens Design System
    '              John Walker   December 1980
    '
    '                  by John Walker
    '              http://www.fourmilab.ch/
    '
    '   This program may be used, distributed, and modified
    '   freely as long as the origin information is preserved.
    '
    '   This is a complete optical design raytracing algorithm,
    '   stripped of its user interface and recast into
    '   FreeBASIC.  It not only determines execution speed on an
    '   extremely floating point (including trig function)
    '   intensive real-world application, it checks accuracy on
    '   an algorithm that is exquisitely sensitive to errors.
    '   The performance of this program is typically far more
    '   sensitive to changes in the efficiency of the
    '   trigonometric library routines than the average floating
    '   point program.
    '
    '   Implemented in November 2017 by John Walker.
    '
    '   This program is based upon the Visual Basic 6
    '   implementation written in July 2005 by Jim White,
    '   mathimagics@yahoo.co.uk, which may be found in the
    '   visualbasic/vb6 directory of this distribution. I have
    '   removed Visual Basic-specific code, restructured the
    '   program, and added comments in the interest of
    '   readability.

    Const max_surfaces = 10

    Dim itercount                  As Long
    Dim niter                      As Long

    Dim Shared current_surfaces    As Integer
    Dim Shared paraxial            As Integer
    Dim Shared clear_aperture      As Double
    Dim Shared aberr_lspher        As Double
    Dim Shared aberr_osc           As Double
    Dim Shared aberr_lchrom        As Double
    Dim Shared max_lspher          As Double
    Dim Shared max_osc             As Double
    Dim Shared max_lchrom          As Double
    Dim Shared radius_of_curvature As Double
    Dim Shared object_distance     As Double
    Dim Shared ray_height          As Double
    Dim Shared axis_slope_angle    As Double
    Dim Shared from_index          As Double
    Dim Shared to_index            As Double

    Dim Shared spectral_line(9)    As Double
    Dim Shared s(max_surfaces, 5)  As Double

    Dim        od_sa(2, 2)         As Double
    Dim        testcase(4, 4)      As Double

    '   Define Cot in terms of Tan

    Function Cot(x as Double) As Double
        Cot =  1.0 / Tan(x)
    End Function

    '   transit_surface propagates a ray across one surface in
    '   the design, calculating its refraction.  The ray and
    '   surface are described by the following global variables.
    '
    '   paraxial                Axial incidence of ray
    '
    '   radius_of_curvature     Radius of curvature of surface
    '                           being crossed.  If 0, surface is
    '                           plane.
    '
    '   object_distance         Distance of object focus
    '                           from lens vertex.  If 0,
    '                           incoming rays are parallel
    '                           and the following must be
    '                           specified:
    '
    '   ray_height              Height of ray from axis.  Only
    '                           relevant if object_distance == 0
    '
    '   axis_slope_angle        Angle incoming ray makes with axis
    '                           at intercept
    '
    '   from_index              Refractive index of medium being left
    '
    '   to_index                Refractive index of medium being
    '                           entered.
    '
    '   The result of a call on transit_surface updates the
    '   following variables modified to reflect the geometry of
    '   the ray as it exits the surface.
    '
    '   object_distance         Distance from vertex to object focus
    '                           after refraction.
    '
    '   ray_height              Height of ray from axis.
    '
    '   axis_slope_angle        Angle incoming ray makes with axis
    '                           at intercept after refraction.
    '
    '   transit_surface returns false when the last surface has been
    '   traversed and true if more surfaces remain to be traced.

    Sub transit_surface()
        Dim iang                 As Double      '   Incidence angle
        Dim rang                 As Double      '   Refraction angle
        Dim iang_sin             As Double      '   Incidence angle sin
        Dim rang_sin             As Double      '   Refraction angle sin
        Dim old_axis_slope_angle As Double
        Dim sagitta              As Double

        If paraxial Then
            If radius_of_curvature <> 0.0 Then
                If object_distance = 0.0 Then
                       axis_slope_angle = 0.0
                       iang_sin = ray_height / radius_of_curvature
                    Else
                        iang_sin = ((object_distance - radius_of_curvature) / _
                            radius_of_curvature) * axis_slope_angle
                    End If
                    rang_sin = (from_index / to_index) * iang_sin
                    old_axis_slope_angle = axis_slope_angle
                    axis_slope_angle = axis_slope_angle + iang_sin - rang_sin
                    If object_distance <> 0.0 Then ray_height = object_distance * _
                        old_axis_slope_angle
                    object_distance = ray_height / axis_slope_angle
                    Exit Sub
                End If

                object_distance = object_distance * (to_index / from_index)
                axis_slope_angle = axis_slope_angle * (from_index / to_index)
                Exit Sub
        End If

        '   Marginal ray

        If radius_of_curvature <> 0.0 Then
            If object_distance = 0.0 Then
                axis_slope_angle = 0#
                iang_sin = ray_height / radius_of_curvature
            Else
                iang_sin = ((object_distance - radius_of_curvature) / _
                    radius_of_curvature) * Sin(axis_slope_angle)
            End If
            iang = aSin(iang_sin)
            rang_sin = (from_index / to_index) * iang_sin
            old_axis_slope_angle = axis_slope_angle
            axis_slope_angle = axis_slope_angle + iang - aSin(rang_sin)
            sagitta = Sin((old_axis_slope_angle + iang) / 2.0)
            sagitta = 2.0 * radius_of_curvature * sagitta * sagitta
            object_distance = ((radius_of_curvature * _
                Sin(old_axis_slope_angle + iang)) * _
                Cot(axis_slope_angle)) + sagitta
            Exit Sub
        End If

        rang = -aSin((from_index / to_index) * Sin(axis_slope_angle))
        object_distance = object_distance * ((to_index * Cos(-rang)) / _
            (from_index * Cos(axis_slope_angle)))
        axis_slope_angle = -rang
    End Sub

    '   trace_line propagates a ray with wavelength Linex and
    '   an initial ray height of ray_h through the surfaces
    '   and leaves the result of the trace in object_distance
    '   and axis_slope_angle.

    Sub trace_line(Linex As Integer, ray_h As Double)
        Dim i As Integer
        object_distance = 0.0
        ray_height = ray_h
        from_index = 1.0

        For i = 1 To current_surfaces
            radius_of_curvature = s(i, 1)
            to_index = s(i, 2)
            If to_index > 1.0 Then
                to_index = to_index + ((spectral_line(4) - spectral_line(Linex)) / _
                    (spectral_line(3) - spectral_line(6))) * ((s(i, 2) - 1.0) / s(i, 3))
            End If
            transit_surface()
            from_index = to_index
            If (i < current_surfaces) Then object_distance = object_distance - s(i, 4)
        Next
    End Sub

    '   The test case used in this program is the design for a 4
    '   inch achromatic telescope objective used as the example
    '   in Wyld's classic work on ray tracing by hand, given in
    '   Amateur Telescope Making, Volume 3 (Volume 2 in the 1996
    '   reprint edition).

    '   Surface 1
    testcase(0, 0) = 27.05      '   Radius of curvature
    testcase(0, 1) = 1.5137     '   Index of Refraction
    testcase(0, 2) = 63.6       '   Dispersion
    testcase(0, 3) = 0.52       '   Edge Thickness

    '   Surface 2
    testcase(1, 0) = -16.68     '   Radius of curvature
    testcase(1, 1) = 1          '   Index of Refraction
    testcase(1, 2) = 0          '   Dispersion
    testcase(1, 3) = 0.138      '   Edge Thickness

    '   Surface 3
    testcase(2, 0) = -16.68     '   Radius of curvature
    testcase(2, 1) = 1.6164     '   Index of Refraction
    testcase(2, 2) = 36.7       '   Dispersion
    testcase(2, 3) = 0.38       '   Edge Thickness

    '   Surface 4
    testcase(3, 0) = -78.1      '   Radius of curvature
    testcase(3, 1) = 1          '   Index of Refraction
    testcase(3, 2) = 0          '   Dispersion
    testcase(3, 3) = 0          '   Edge Thickness

    '  Wavelengths of standard spectral lines in Angstroms
    '         (Not all are used in this program)

    spectral_line(1) = 7621.0    '   A
    spectral_line(2) = 6869.955  '   B
    spectral_line(3) = 6562.816  '   C
    spectral_line(4) = 5895.944  '   D
    spectral_line(5) = 5269.557  '   E
    spectral_line(6) = 4861.344  '   F
    spectral_line(7) = 4340.477  '   G'
    spectral_line(8) = 3968.494  '   H

    ' Load test case into working array

    clear_aperture = 4.0
    current_surfaces = 4

    Dim as integer i, j
    For i = 0 To current_surfaces - 1
       For j = 0 To 3
          s(i + 1, j + 1) = testcase(i, j)
          Next
       Next

    Dim as integer k
    Dim od_fline As Double, od_cline As Double

    niter = 1000

    '   Allow iteration count to be specified on the command line
    If __FB_ARGC__ > 1 Then
        niter = CLng(Command(1))
    End If

    ' Perform ray trace the specified number of times.

    For itercount = 1 To niter

        For paraxial = 0 To 1

             ' Do main trace in D light

             trace_line(4, clear_aperture / 2.0)
             od_sa(paraxial, 0) = object_distance
             od_sa(paraxial, 1) = axis_slope_angle
        Next

        paraxial = 0

        ' Trace marginal ray in C

        trace_line(3, clear_aperture / 2.0)
        od_cline = object_distance

        ' Trace marginal ray in F

        trace_line(6, clear_aperture / 2.0)
        od_fline = object_distance

        '   Compute aberrations of the design

        '   The longitudinal spherical aberration is just the
        '   difference between where the D line comes to focus
        '   for paraxial and marginal rays.

        aberr_lspher = od_sa(1, 0) - od_sa(0, 0)

        '   The offense against the sine condition is a measure
        '   of the degree of coma in the design.  We compute it
        '   as the lateral distance in the focal plane between
        '   where a paraxial ray and marginal ray in the D line
        '   come to focus.

        aberr_osc = 1.0 - (od_sa(1, 0) * od_sa(1, 1)) / _
            (Sin(od_sa(0, 1)) * od_sa(0, 0))

        '   The axial chromatic aberration is the distance between
        '   where marginal rays in the C and F lines come to focus.

        aberr_lchrom = od_fline - od_cline
        max_lspher = Sin(od_sa(0, 1))

        '   Compute maximum acceptable values for each aberration

        '   Maximum longitudinal spherical aberration, which is
        '   also the maximum for axial chromatic aberration.  This
        '   is computed for the D line.

        max_lspher = 0.0000926 / (max_lspher * max_lspher)
        max_osc = 0.0025
        max_lchrom = max_lspher     ' Same criterion

    Next

    '   Print the result in the standard format

    Const mp = "    (Maximum permissible):               ###.###########"
    Print using "   Marginal ray         ###.########### ###.###########"; _
        od_sa(0, 0); od_sa(0, 1)
    Print using "   Paraxial ray         ###.########### ###.###########"; _
        od_sa(1, 0); od_sa(1, 1)
    Print using "Longitudinal spherical aberration:       ###.###########"; _
        aberr_lspher
    Print using mp; max_lspher
    Print using "Offense against sine condition (coma):   ###.###########"; _
        aberr_osc
    Print using mp;  max_osc
    Print using "Axial chromatic aberration:              ###.###########"; _
        aberr_lchrom
    Print using mp; max_lchrom

    End
