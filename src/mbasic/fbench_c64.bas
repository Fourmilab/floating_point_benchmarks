1000 rem  john walker's optical bench
1010 rem
1020 rem  microsoft basic version
1022 rem
1024 rem  configured for commodore 62 basic v2
1030 rem
1040 rem  this must be run on a basic which provides double
1050 rem  precision, and in which all trig functions are done
1060 rem  in double precision.  if not, the results are utterly
1070 rem  useless.  when the test case is run, the results should
1080 rem  agree with those that follow to the last digit, or the
1090 rem  accuracy of the host system is suspect.
1100 rem
1110 rem                             focal length       angle to axis
1120 rem marginal ray               47.09479120920      0.04178472683
1130 rem paraxial ray               47.08372160249      0.04177864821
1140 rem
1150 rem longitudinal spherical aberration:       -0.01106960671
1160 rem     (maximum permissible):                0.05306749907    acceptable
1170 rem
1180 rem  offense against sine condition (coma):   0.00008954761
1190 rem     (maximum permissible):                0.00250000000    acceptable
1200 rem
1210 rem axial chromatic aberration:               0.00448229032
1220 rem     (maximum permissible):                0.05306749907    acceptable
1230 rem
1240 rem  correspondence of variable names in this program to
1250 rem  variable names in the original qbasic program.
1260 rem
1270 rem s1     spectral.line
1280 rem z      zline
1290 rem s3     od.sa
1300 rem d$     dname$
1310 rem i1     icurrent.surfaces
1320 rem c1     clear.aperture
1330 rem p      iparaxial
1340 rem h      ray.h
1350 rem h2     ray.height
1360 rem o1     object.distance
1370 rem a2     axis.slope.angle
1380 rem o2     od.cline
1390 rem o3     od.fline
1400 rem b1     aberr.lspher
1410 rem b2     aberr.lchrom
1420 rem b3     aberr.osc
1430 rem e1     max.lspher
1440 rem e2     max.lchrom
1450 rem e3     max.osc
1460 rem x1     from.index
1470 rem x2     to.index
1480 rem v1     liang.sin
1490 rem v2     rang.sin
1500 rem v3     liang
1510 rem v4     rang
1520 rem c2     radius.of.curvature
1530 rem g1     sagitta
1540 rem a3     old.axis.slope.angle
1550 dim s1(9)
1560 dim s(10, 5)
1570 dim s3(2,2)
1580 rem
1590 def fna(x) = atn(x / sqr(1 - x * x))
1600 rem
1610 rem  initialise when called the first time
1620 rem
1630 s1(1) = 7621      :  rem a
1640 s1(2) = 6869.955  :  rem b
1650 s1(3) = 6562.816  :  rem c
1660 s1(4) = 5895.944  :  rem d
1670 s1(5) = 5269.557  :  rem e
1680 s1(6) = 4861.344  :  rem f
1690 s1(7) = 4340.477  :  rem g
1700 s1(8) = 3968.494  :  rem h
1710 rem
1720 read d$
1730 read i1
1740 read c1
1750 for i = 1 to i1
1760    read s(i, 1), s(i, 2), s(i, 3), s(i, 4)
1770 next i
1780 rem
1790 print "surface         radius       ref. index    dispersion   thickness"
1800 print
1810 for i = 1 to i1
1820    print i, s(i,1), s(i,2), s(i,3), s(i,4)
1830 next i
1832 input "how many iterations"; i9
1840 rem
1841 print "press return to begin: ";
1842 get s9$ : if s9$ ="" then 1842
1843 print
1844 for i8 = 1 to i9
1850 for p = 0 to 1
1860    rem do main trace in d light
1870    z = 4
1880    h = c1 / 2
1890    gosub 2410
1900    s3(p,0) = o1
1910    s3(p,1) = a2
1920 next p
1930 p = 0
1940 rem
1950 rem trace marginal ray in c
1960 z = 3
1970 h = c1 / 2
1980 gosub 2410
1990 o2 = o1
2000 rem
2010 rem trace marginal ray in f
2020 z = 6
2030 h = c1 / 2
2040 gosub 2410
2050 o3 = o1
2060 rem
2070 b1 = s3(1,0) - s3(0,0)
2080 b3 = 1 - (s3(1,0) * s3(1,1)) / (sin(s3(0,1)) * s3(0,0))
2090 b2 = o3 - o2
2100 e1 = sin(s3(0,1))
2110 rem d light
2120 e1 = .0000926 / (e1 * e1)
2130 e3 = .0025
2140 e2 = e1
2141 next i8
2143 print "stop the timer: ";
2144 get s9$ : if s9$ = "" then 2144
2145 print
2150 rem
2170 rem
2180 print
2190 print "            ", "focal length", "angle to axis"
2200 print "marginal ray", s3(0,0), s3(0,1)
2210 print "paraxial ray", s3(1,0), s3(1,1)
2220 print
2230 print "longitudinal spherical aberration: ", b1
2240 print "  (maximum permissible): ", e1;
2250 if abs(b1) <= e1 then print "acceptable" : goto 2260 
2252 print "** excessive **"
2260 print
2270 rem
2280 print "offense against sine condition (coma): ", b3
2290 print "  (maximum permissible): ", e3;
2300 if abs(b3) <= e3 then print "acceptable" : goto 2310 
2302 print "** excessive **"
2310 print
2320 rem
2330 print "axial chromatic aberration: ", b2
2340 print "  (maximum permissible): ", e2;
2350 if abs(b2) <= e2 then print "acceptable" : goto 2360 
2352 print "** excessive **"
2360 rem
2370 end
2380 rem
2390 rem  perform ray trace in specific spectral line
2400 rem
2410 rem trace.l(z, h)
2420 rem
2430         o1 = 0
2440         h2 = h
2450         x1 = 1
2460 rem
2470         for i = 1 to i1
2480            c2 = s(i,1)
2490            x2 = s(i,2)
2500            if x2 <= 1 then  2510
2502            x2=x2+((s1(4)-s1(z))/(s1(3)-s1(6)))*((s(i,2)-1)/s(i,3))
2510            gosub 2630
2520            x1 = x2
2530            if i < i1 then o1 = o1 - s(i,4)
2540         next i
2550         return
2560 rem
2570 rem            calculate passage through surface
2580 rem
2590 rem             if the variable p is true, the trace through the
2600 rem             surface will be done using the p approximations.
2610 rem             otherwise, the normal trigonometric trace will be done.
2620 rem
2630 rem transit.surface()
2640 rem               v3   incidence angle
2650 rem               v4   refraction angle
2660 rem               v1   incidence angle sin
2670 rem               v2   refraction angle sin
2680 rem
2690         if p = 0 then goto 2820
2700            if c2 = 0 then goto 2780
2710               if o1 = 0 then a2 = 0 : v1 = h2/c2 : goto 2720 
2712               v1 = ((o1 - c2) / c2) * a2
2720               v2 = (x1 / x2) * v1
2730               a3 = a2
2740               a2 = a2 + v1 - v2
2750               if o1 <> 0 then h2 = o1 * a3
2760               o1 = h2 / a2
2770               return
2780            o1 = o1 * (x2 / x1)
2790            a2 = a2 * (x1 / x2)
2800            return
2810 rem
2820         if c2 = 0 then goto 2930
2830            if  o1 = 0 then a2 = 0 : v1 = h2/c2 : goto 2840 
2832            v1 = ((o1 - c2) / c2) * sin(a2)
2840            v3 = fna(v1)
2850 rem
2860            v2 = (x1 / x2) * v1
2870            a3 = a2
2880            a2 = a2 + v3 - fna(v2)
2890            g1 = sin((a3 + v3) / 2)
2900            g1 = 2 * c2*g1*g1
2910            o1 = ((c2 * sin(a3 + v3)) * (1 / tan(a2))) + g1
2920            return
2930         v4 = -fna((x1 / x2) * sin(a2))
2940         o1 = o1 * ((x2 * cos(-v4)) / (x1 * cos(a2)))
2950         a2 = -v4
2960         return
2970 rem
2980 rem       the design is defined in the following data statements
2990 rem
3000 rem                  design name
3010         data "4 inch f/12 achromatic objective"
3020 rem          number of surfaces, clear aperture
3030         data 4, 4.0
3040 rem  for each surface:  radius of curvature (+ if convex to light source)
3050 rem                     index of refraction (1 for air space)
3060 rem                     dispersion (abbe number (v))
3070 rem                     edge thickness (0 for last surface)
3080         data  27.05, 1.5137, 63.6, .52
3090         data -16.68, 1, 0, .138
3100         data -16.68, 1.6164, 36.7, .38
3110         data -78.1, 1, 0, 0
3120 end
