Option Strict Off
Option Explicit On
Module modFBENCH
	
	Dim time1, time2 As Integer
	Dim data1, data2 As Double
	Dim itns As Integer
	
	'===================================================================
	'  FBENCH for Visual Basic 6
	'
	'    Conversion by Jim White, mathimagics@yahoo.co.uk
	'       July 2005
	'
	'===================================================================
	'
	'  Converted from fbench.c (John Walker's Floating Point Benchmark,
	'    'C' version), by
	'
	'                John Walker   December 1980
	'                Marinchip Interactive Lens Design System
	'                http://www.fourmilab.ch/
	'
	'  This  program may be used, distributed, and modified freely as
	'  long as the origin information is preserved.
	'
	'  This  is  a  complete  optical   design   raytracing  algorithm,
	'  stripped of its user interface and recast into portable C.  It
	'  not only determines execution speed on an  extremely  floating
	'  point (including   trig   function) intensive   real-world
	'  application, it  checks  accuracy  on  an  algorithm  that  is
	'  exquisitely  sensitive  to  errors.   The  performance of this
	'  program is typically far more  sensitive  to  changes  in  the
	'  efficiency  of the  trigonometric  library  routines than the
	'  average floating point program.
	'
	'===================================================================
	
	Const max_surfaces As Short = 10
	
	Dim itercount As Integer
	Dim niter As Integer
	
	Dim current_surfaces As Short
	Dim paraxial As Short
	Dim clear_aperture As Double
	Dim aberr_lspher As Double
	Dim aberr_osc As Double
	Dim aberr_lchrom As Double
	Dim max_lspher As Double
	Dim max_osc As Double
	Dim max_lchrom As Double
	Dim radius_of_curvature As Double
	Dim object_distance As Double
	Dim ray_height As Double
	Dim axis_slope_angle As Double
	Dim from_index As Double
	Dim to_index As Double
	
	Dim spectral_line(9) As Double
	Dim s(10, 5) As Double ' [max_surfaces][5]
	Dim od_sa(2, 2) As Double
	Dim testcase(4, 4) As Double
	
	Dim tstart, tfinish As Integer ' GetTickCount start/end values
	' (not the most accurate clock for small
	'  interval measuring, but it's fine if
	'  NITER is big enough to make the
	'  test take at least 5 seconds)
	
	'
	'===================================================================
	' Reference data (64-bit mode):
	'
	'       "   Marginal ray          47.09479120920   0.04178472683",
	'       "   Paraxial ray          47.08372160249   0.04177864821",
	'       "Longitudinal spherical aberration:        -0.01106960671",
	'       "    (Maximum permissible):                 0.05306749907",
	'       "Offense against sine condition (coma):     0.00008954761",
	'       "    (Maximum permissible):                 0.00250000000",
	'       "Axial chromatic aberration:                0.00448229032",
	'       "    (Maximum permissible):                 0.05306749907"
	
	'  The   test  case  used  in  this program is the  design for a 4 inch
	'  achromatic telescope  objective  used  as  the  example  in  Wyld's
	'  classic  work  on  ray  tracing by hand, given in Amateur Telescope
	'  Making, Volume 3.
	'===================================================================
	
    Private Declare Function GetTickCount Lib "kernel32" () As Integer

	
	Private Sub Initialise()
		Dim i, j As Integer
		testcase(0, 0) = 27.05 : testcase(0, 1) = 1.5137 : testcase(0, 2) = 63.6 : testcase(0, 3) = 0.52
		testcase(1, 0) = -16.68 : testcase(1, 1) = 1 : testcase(1, 2) = 0 : testcase(1, 3) = 0.138
		testcase(2, 0) = -16.68 : testcase(2, 1) = 1.6164 : testcase(2, 2) = 36.7 : testcase(2, 3) = 0.38
		testcase(3, 0) = -78.1 : testcase(3, 1) = 1
		
		spectral_line(1) = 7621# '   /* A */
		spectral_line(2) = 6869.955 '   /* B */
		spectral_line(3) = 6562.816 '   /* C */
		spectral_line(4) = 5895.944 '   /* D */
		spectral_line(5) = 5269.557 '   /* E */
		spectral_line(6) = 4861.344 '   /* F */
		spectral_line(7) = 4340.477 '   /* G'*/
		spectral_line(8) = 3968.494 '   /* H */
		
		' Load test case into working array
		
		clear_aperture = 4#
		current_surfaces = 4
		For i = 0 To current_surfaces - 1
			For j = 0 To 3
				s(i + 1, j + 1) = testcase(i, j)
			Next 
		Next 
        Form1.DefInstance.IsInitializing = False
    End Sub


    Sub fbench_run()
        Dim k, i, j, errors As Integer
        Dim od_fline, od_cline As Double

        Initialise()
        niter = Val(Form1.DefInstance.Text1.Text)

        Form1.DefInstance.List1.Items.Add(TimeString & "   Test iterations = " & niter)
        Form1.DefInstance.List1.Items.Add(TimeString & "   Starting test")
        System.Windows.Forms.Application.DoEvents()
        tstart = GetTickCount()
        '
        ' Perform ray trace the specified number of times.
        '

        For itercount = 1 To niter

            For paraxial = 0 To 1
                '
                ' Do main trace in D light
                '
                trace_line(4, clear_aperture / 2.0#)
                od_sa(paraxial, 0) = object_distance
                od_sa(paraxial, 1) = axis_slope_angle
            Next

            paraxial = 0

            ' Trace marginal ray in C

            trace_line(3, clear_aperture / 2.0#)
            od_cline = object_distance

            ' Trace marginal ray in F

            trace_line(6, clear_aperture / 2.0#)
            od_fline = object_distance

            aberr_lspher = od_sa(1, 0) - od_sa(0, 0)
#If INTRIG Then
			'UPGRADE_NOTE: #If #EndIf block was not upgraded because the expression INTRIG did not evaluate to True or was not evaluated. Click for more: 'ms-help://MS.VSCC/commoner/redir/redirect.htm?keyword="vbup1035"'
			aberr_osc = 1# - (od_sa(1, 0) * od_sa(1, 1)) / (Sin(od_sa(0, 1)) * od_sa(0, 0))
			aberr_lchrom = od_fline - od_cline
			max_lspher = Sin(od_sa(0, 1))
#Else
            aberr_osc = 1.0# - (od_sa(1, 0) * od_sa(1, 1)) / (fSin(od_sa(0, 1)) * od_sa(0, 0))
            aberr_lchrom = od_fline - od_cline
            max_lspher = fSin(od_sa(0, 1))
#End If
            ' D light

            max_lspher = 0.0000926 / (max_lspher * max_lspher)
            max_osc = 0.0025
            max_lchrom = max_lspher
        Next
        tfinish = GetTickCount()

        With Form1.DefInstance.List1
            .Items.Add(TimeString & "   Test completed, " & Format((tfinish - tstart) / 1000.0#, "#.000") & " secs")
            .Items.Add("------------------------------------------")
            .Items.Add("Marginal ray   " & Format(od_sa(0, 0), "#.00000000000") & "     " & Format(od_sa(0, 1), "#.00000000000"))
            .Items.Add("Paraxial ray   " & Format(od_sa(1, 0), "#.00000000000") & "     " & Format(od_sa(1, 1), "#.00000000000"))
            .Items.Add("Longitudinal spherical aberration:     " & Format(aberr_lspher, "#.00000000000"))
            .Items.Add("    (max permissible)            :      " & Format(max_lspher, "#.00000000000"))
            .Items.Add("Offense against sine condition (coma):  " & Format(aberr_osc, "#.00000000000"))
            .Items.Add("    (max permissible)            :      " & Format(max_osc, "#.00000000000"))
            .Items.Add("Axial chromatic aberration:             " & Format(aberr_lchrom, "#.00000000000"))
            .Items.Add("    (max permissible)            :      " & Format(max_lchrom, "#.00000000000"))
            .Items.Add("------------------------------------------")
            .Items.Add("")
            .Items.Add("Time per 1000 iterations = " & Format((tfinish - tstart) / niter, "#.0000") & " secs")
        End With
    End Sub

    Sub trace_line(ByVal Line As Short, ByVal ray_h As Double)
        '
        '  Perform ray trace in specific spectral line
        '
        Dim i As Short
        object_distance = 0.0#
        ray_height = ray_h
        from_index = 1.0#

        For i = 1 To current_surfaces
            radius_of_curvature = s(i, 1)
            to_index = s(i, 2)
            If to_index > 1.0# Then
                to_index = to_index + ((spectral_line(4) - spectral_line(Line)) / (spectral_line(3) - spectral_line(6))) * ((s(i, 2) - 1.0#) / s(i, 3))
            End If
            transit_surface()
            from_index = to_index
            If (i < current_surfaces) Then object_distance = object_distance - s(i, 4)
        Next
    End Sub
#If INTRIG Then
	'UPGRADE_NOTE: #If #EndIf block was not upgraded because the expression INTRIG did not evaluate to True or was not evaluated. Click for more: 'ms-help://MS.VSCC/commoner/redir/redirect.htm?keyword="vbup1035"'
	Sub transit_surface()

	Dim iang                 As Double      ' Incidence angle
	Dim rang                 As Double      ' Refraction angle
	Dim iang_sin             As Double      ' Incidence angle sin
	Dim rang_sin             As Double      ' Refraction angle sin
	Dim old_axis_slope_angle As Double
	Dim sagitta              As Double

	' Calculate passage through surface
	'
	'   If  the variable PARAXIAL is true, the trace through the
	'   surface will be done using the paraxial  approximations.
	'   Otherwise,  the normal trigonometric trace will be done.
	'
	' This routine takes the following inputs:
	'
	'   RADIUS_OF_CURVATURE    Radius of curvature of surface
	'           being crossed.  If 0, surface is a plane.
	'
	'   OBJECT_DISTANCE        Distance of object focus from
	'           lens vertex. If 0, incoming rays are parallel
	'           and the following must be specified:
	'
	'      RAY_HEIGHT       Height of ray from axis.  Only
	'                      relevant if OBJECT.DISTANCE == 0
	'
	'         AXIS_SLOPE_ANGLE       Angle incoming ray makes with axis
	'                 at intercept
	'
	'         FROM_INDEX       Refractive index of medium being left
	'
	'         TO_INDEX         Refractive index of medium being
	'                 entered.
	'
	'         The outputs are the following variables:
	'
	'         OBJECT_DISTANCE        Distance from vertex to object focus
	'                 after refraction.
	'
	'         AXIS_SLOPE_ANGLE       Angle incoming ray makes with axis
	'                 at intercept after refraction.'



	If paraxial Then
	If radius_of_curvature <> 0# Then

	If object_distance = 0# Then
	axis_slope_angle = 0#
	iang_sin = ray_height / radius_of_curvature
	Else
	iang_sin = ((object_distance - radius_of_curvature) / radius_of_curvature) * axis_slope_angle
	End If
	rang_sin = (from_index / to_index) * iang_sin
	old_axis_slope_angle = axis_slope_angle
	axis_slope_angle = axis_slope_angle + iang_sin - rang_sin
	If object_distance <> 0# Then ray_height = object_distance * old_axis_slope_angle
	object_distance = ray_height / axis_slope_angle
	Exit Sub
	End If
	object_distance = object_distance * (to_index / from_index)
	axis_slope_angle = axis_slope_angle * (from_index / to_index)
	Exit Sub
	End If

	If radius_of_curvature <> 0# Then
	If object_distance = 0# Then
	axis_slope_angle = 0#
	iang_sin = ray_height / radius_of_curvature
	Else
	iang_sin = ((object_distance - radius_of_curvature) / radius_of_curvature) * Sin(axis_slope_angle)
	End If
	iang = aSin(iang_sin)
	rang_sin = (from_index / to_index) * iang_sin
	old_axis_slope_angle = axis_slope_angle
	axis_slope_angle = axis_slope_angle + iang - aSin(rang_sin)
	sagitta = Sin((old_axis_slope_angle + iang) / 2#)
	sagitta = 2# * radius_of_curvature * sagitta * sagitta
	object_distance = ((radius_of_curvature * Sin(old_axis_slope_angle + iang)) * Cot(axis_slope_angle)) + sagitta
	Exit Sub
	End If

	rang = -aSin((from_index / to_index) * Sin(axis_slope_angle))
	object_distance = object_distance * ((to_index * Cos(-rang)) / (from_index * Cos(axis_slope_angle)))
	axis_slope_angle = -rang
	End Sub

#Else
    Sub transit_surface()

        Dim iang As Double ' Incidence angle
        Dim rang As Double ' Refraction angle
        Dim iang_sin As Double ' Incidence angle sin
        Dim rang_sin As Double ' Refraction angle sin
        Dim old_axis_slope_angle As Double
        Dim sagitta As Double

        If paraxial Then
            If radius_of_curvature <> 0.0# Then

                If object_distance = 0.0# Then
                    axis_slope_angle = 0.0#
                    iang_sin = ray_height / radius_of_curvature
                Else
                    iang_sin = ((object_distance - radius_of_curvature) / radius_of_curvature) * axis_slope_angle
                End If
                rang_sin = (from_index / to_index) * iang_sin
                old_axis_slope_angle = axis_slope_angle
                axis_slope_angle = axis_slope_angle + iang_sin - rang_sin
                If object_distance <> 0.0# Then ray_height = object_distance * old_axis_slope_angle
                object_distance = ray_height / axis_slope_angle
                Exit Sub
            End If
            object_distance = object_distance * (to_index / from_index)
            axis_slope_angle = axis_slope_angle * (from_index / to_index)
            Exit Sub
        End If

        If radius_of_curvature <> 0.0# Then
            If object_distance = 0.0# Then
                axis_slope_angle = 0.0#
                iang_sin = ray_height / radius_of_curvature
            Else
                iang_sin = ((object_distance - radius_of_curvature) / radius_of_curvature) * fSin(axis_slope_angle)
            End If
            iang = aSin(iang_sin)
            rang_sin = (from_index / to_index) * iang_sin
            old_axis_slope_angle = axis_slope_angle
            axis_slope_angle = axis_slope_angle + iang - aSin(rang_sin)
            sagitta = fSin((old_axis_slope_angle + iang) / 2.0#)
            sagitta = 2.0# * radius_of_curvature * sagitta * sagitta
            object_distance = ((radius_of_curvature * fSin(old_axis_slope_angle + iang)) * Cot(axis_slope_angle)) + sagitta
            Exit Sub
        End If

        rang = -aSin((from_index / to_index) * fSin(axis_slope_angle))
        object_distance = object_distance * ((to_index * System.Math.Cos(-rang)) / (from_index * fCos(axis_slope_angle)))
        axis_slope_angle = -rang
    End Sub
#End If
End Module