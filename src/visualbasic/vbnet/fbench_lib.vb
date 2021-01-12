Option Strict Off
Option Explicit On
Module modFBLIB
	
    '  VB6 intrinsic trig functions are used when INTRIG is true.  The
	'    functions supplied here may be enabled to obtain timings that reflect
	'    the machine 's floating point performance rather than the speed  of
	'    its trig function evaluation. To use these functions set INTRIG = False
	
#Const INTRIG = True
	
	Const pi As Double = 3.14159265358979
	Public Const piover2 As Double = pi / 2#
	
	Const twopi As Double = pi * 2#
	Const piover4 As Double = pi / 4#
	Const fouroverpi As Double = 4# / pi
	
	
	Dim atanc(8) As Double
	
	Function fabs(ByVal x As Double) As Double
		If x < 0# Then fabs = -x Else fabs = x
	End Function
	
	Sub LibInitialise()
		'
		'  Coefficients for ATAN evaluation
		'
		atanc(1) = 0.463647609000806 '0.4636476090008061165
		atanc(2) = 0.785398163397448 '0.7853981633974483094,
		atanc(3) = 0.982793723247329 '0.98279372324732906714
		atanc(4) = 0.982793723247329 '1.1071487177940905022
		atanc(5) = 1.19028994968253 '1.1902899496825317322
		atanc(6) = 1.24904577239825 '1.2490457723982544262
		atanc(7) = 1.29249666778979 '1.2924966677897852673
		atanc(8) = 1.32581766366803 '1.3258176636680324644
	End Sub
	
	
	Function aInt(ByVal x As Double) As Double
		'
		' aint(x)   Return integer part of number.  Truncates towards 0
		'
		Dim l As Integer
		'
		' Note that this routine cannot handle the full floating point
		'   number range.  This function should be in the machine-dependent
		'   floating point library!
		'
		l = x
		If (CInt(-0.5) <> 0) And (l < 0) Then l = l + 1
		aInt = l
	End Function
	
	Function fSin(ByVal x As Double) As Double
		'
		' sin(x)    Return sine, x in radians
		'
		Dim sign As Boolean
		Dim r, Y, z As Double
		
		sign = x < 0#
		
		If sign Then x = -x
		
		If x > twopi Then x = x - (aInt(x / twopi) * twopi)
		
		If (x > pi) Then x = x - pi : sign = Not sign
		
		If (x > piover2) Then x = pi - x
		
		If (x < piover4) Then
			Y = x * fouroverpi
			z = Y * Y
			r = Y * (((((((-2.02253129293E-14 * z + 6.9481520350522E-12) * z - 1.75724741761708E-09) * z + 3.13361688917325E-07) * z - 3.65762041821464E-05) * z + 2.49039457019272E-03) * z - 8.07455121882808E-02) * z + 0.785398163397448)
		Else
			Y = (piover2 - x) * fouroverpi
			z = Y * Y
			r = ((((((-3.8577620372E-13 * z + 1.1500497024263E-10) * z - 2.461136382637E-08) * z + 3.59086044588582E-06) * z - 3.25991886926688E-04) * z + 1.58543442438154E-02) * z - 0.308425137534042) * z + 1#
		End If
		If sign Then fSin = -r Else fSin = r
	End Function
	
	
	Function fCos(ByVal x As Double) As Double
		'
		'  cos(x)    Return cosine, x in radians, by identity
		'
		If x < 0# Then x = -x
		
		If x > twopi Then ' Do range reduction here to limit
			x = x - (aInt(x / twopi) * twopi) ' roundoff on add of PI/2
		End If
		fCos = fSin(x + piover2)
	End Function
	
	
	Function fTan(ByVal x As Double) As Double
		fTan = fSin(x) / fCos(x)
	End Function
	
	
	
	Function fSqrt(ByVal x As Double) As Double
		'
		'  sqrt(x)   Return square root.  Initial guess, then Newton-  Raphson refinement
		
		Dim cl, C, Y As Double
		Dim n As Integer
		
		If x = 0# Then Exit Function
		
		If x < 0# Then Err.Raise(14,  , "Invalid argument for fSqr()")
		
		Y = (0.154116 + 1.893872 * x) / (1# + 1.047988 * x)
		
		C = (Y - x / Y) / 2#
		cl = 0#
		n = 50
		Do While C <> cl
			Y = Y - C
			cl = C
			C = (Y - x / Y) / 2#
			n = n - 1
			If n = 0 Then Exit Do
		Loop 
		fSqrt = Y
	End Function
	
	
	Function aTan(ByVal x As Double) As Double
		'
		' atan(x)   Return arctangent in radians, range -pi/2 to pi/2
		'
		Dim sign As Short
		Dim l, Y As Integer
		Dim b, a, z As Double
		
		sign = (x < 0#)
		If sign Then x = -x
		
		If x >= 4# Then
			l = -1
			x = 1# / x
			GoTo atl
		ElseIf (x < 0.25) Then 
			GoTo atl
		End If
		
		Y = aInt(x / 0.5)
		z = Y * 0.5
		x = (x - z) / (x * z + 1)
		
atl: 
		z = x * x
		b = ((((893025# * z + 49116375#) * z + 425675250#) * z + 1277025750#) * z + 1550674125#) * z + 654729075#
		a = (((13852575# * z + 216602100#) * z + 891080190#) * z + 1332431100#) * z + 654729075#
		a = (a / b) * x + atanc(Y)
		If (l) Then a = piover2 - a
		If sign Then aTan = -a Else aTan = a
	End Function
	
	
	Function aTan2(ByVal Y As Double, ByRef x As Double) As Double
		'
		'  atan2(y,x)   Return arctangent in radians of y/x,
		'     range -pi to pi
		'
		Dim temp As Double
		
		If x = 0# Then
			If Y = 0# Then Exit Function '  Special case: atan2(0,0) = 0
			If (Y > 0) Then aTan2 = piover2 Else aTan2 = -piover2
			Exit Function
		End If
		
		temp = aTan(Y / x)
		
		If x < 0# Then
			If Y >= 0# Then
				temp = temp + pi
			Else
				temp = temp - pi
			End If
		End If
		aTan2 = temp
	End Function
	
	Function aSin(ByVal x As Double) As Double
		'
		'  asin(x)   Return arcsine in radians of x
		'
		If fabs(x) > 1# Then Err.Raise(14,  , "Invalid argument for aSin()")
		
#If INTRIG Then
		aSin = aTan2(x, System.Math.Sqrt(1 - x * x))
#Else
		'UPGRADE_NOTE: #If #EndIf block was not upgraded because the expression Else did not evaluate to True or was not evaluated. Click for more: 'ms-help://MS.VSCC/commoner/redir/redirect.htm?keyword="vbup1035"'
		aSin = aTan2(x, fSqrt(1 - x * x))
#End If
	End Function
	
	
	Function Cot(ByVal x As Double) As Double
#If INTRIG Then
		Cot = (1# / System.Math.Tan(x))
#Else
		'UPGRADE_NOTE: #If #EndIf block was not upgraded because the expression Else did not evaluate to True or was not evaluated. Click for more: 'ms-help://MS.VSCC/commoner/redir/redirect.htm?keyword="vbup1035"'
		Cot = (1# / fTan(x))
#End If
	End Function
End Module