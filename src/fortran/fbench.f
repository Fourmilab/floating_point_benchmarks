c   John Walker's Floating Point Benchmark, derived from...
c
c        Marinchip Interactive Lens Design System
c               John Walker   December 1980
c
c                   by John Walker
c               http://www.fourmilab.ch/
c
c    This program may be used, distributed, and modified freely as
c    long as the origin information is preserved.
c
c    This is a complete optical design raytracing algorithm,
c    stripped of its user interface and recast into Fortran 77.
c    It not only determines execution speed on an extremely
c    floating point (including trig function) intensive
c    real-world application, it checks accuracy on an algorithm
c    that is exquisitely sensitive to errors. The performance of
c    this program is typically far more sensitive to changes in
c    the efficiency of the trigonometric library routines than the
c    average floating point program.
c
c    When the test case is run, the results should
c    agree with those that follow to the last digit, or the
c    accuracy of the host system is suspect.
c
c                           Focal Length          Angle to Axis
c    Marginal ray         47.09479120920           .04178472683
c    Paraxial ray         47.08372160249           .04177864821
c
c    Longitudinal spherical aberration:         -.01106960671
c        (Maximum permissible):                  .05306749907
c                                                  Acceptable
c
c    Offense against sine condition (coma):      .00008954761
c        (Maximum permissible):                  .00250000000
c                                                  Acceptable
c
c    Axial chromatic  aberration:                 .00448229032
c        (Maximum permissible):                   .05306749907
c                                                  Acceptable
c
c
c    This version was based upon a Microsoft BASIC version of the
c    benchmark.  It is logically equivalent to the original C version,
c    and generates the same results to the last decimal place.  The
c    output differs from the C version in that Fortran formatted
c    output of real and double precision numbers with absolute
c    value less than one does not display a zero to the left of
c    the decimal point.
c
c    The original version was created in July 2005.  In June 2014
c    the code was cleaned up, with gnarly BASIC variable names
c    replaced with less cryptic (within the constraints of Fortran)
c    identifiers.

c       Axial incidence codes
        parameter (Margin = 1, Parax = 2)
c       Spectral line indices
        parameter (LineA = 1, LineB = 2, LineC = 3, LineD = 4,
     1             LineE = 5, LineF = 6, LineGp = 7, LineH = 8)

        double precision clearap, rheight, tlod, tlsa
        double precision cmargd, fmargd, lspher
        double precision lchrom, abosc, mxlsph, maxosc
        double precision  mxcrom, findex, indref, radcur
        integer axiali

        double precision odsa(2,2)
        integer tbfr(80)
        integer niter

        common /dat/ rheight, tlod, tlsa,
     1               findex, indref,
     2               radcur,
     3               axiali

c       Clear aperture of the lens design.  The design is defined
c       within subroutine tline.
        data clearap / 4.0d0 /
c
c       The variable niter specifies the number of iterations to be
c       executed.  Change it and recompile to adjust for reasonable
c       timing on your system.  For archival timings, set the iteration
c       count so the benchmark runs around five minutes.
c
        niter = 10000
c        niter = 200 000 000

        print 1111, niter
1111    format ('Press return to begin ', i8, ' iterations:')
        read (5, 1000) tbfr
1000    format (80A1)
        do 102 i8 = 1, niter

c       Trace marginal and paraxial rays in D light
        do 101 axiali = Margin, Parax
           call tline(LineD, clearap / 2)
           odsa(axiali, 1) = tlod
           odsa(axiali, 2) = tlsa
101     continue
        axiali = Margin

c       Trace marginal ray in C

        call tline(LineC, clearap / 2)
        cmargd = tlod

c       Trace marginal ray in F

        call tline(LineF, clearap / 2)
        fmargd = tlod

c       Compute aberrations of the design

c       The longitudinal spherical aberration is just the
c       difference between where the D line comes to focus
c       for paraxial and marginal rays.
        lspher = odsa(Parax, 1) - odsa(Margin, 1)

c       The offense against the sine condition is a measure
c       of the degree of coma in the design.  We compute it
c       as the lateral distance in the focal plane between
c       where a paraxial ray and marginal ray in the D line
c       come to focus.
        abosc = 1 - (odsa(Parax, 1) * odsa(Parax, 2)) /
     1       (sin(odsa(Margin, 2)) * odsa(Margin, 1))

c       The axial chromatic aberration is the distance between
c       where marginal rays in the C and F lines come to focus.
        lchrom = fmargd - cmargd

c       Compute maximum acceptable values for each aberration

c       Maximum longitudinal spherical aberration, which is
c       also the maximum for axial chromatic aberration.  This
c       is computed for the D line.
        mxlsph = .0000926d0 / (sin(odsa(Margin, 2)) ** 2)
        mxcrom = .0025d0
        maxosc = mxlsph
102     continue
        print *, 'Stop the timer:'
        read (5, 1000) tbfr

        print 1108
1108    format (23x, 'Focal Length          Angle to Axis')
        print 1100, odsa(1,1), odsa(1,2)
1100    format ('Marginal ray  ',f21.11,'  ',f21.11)
        print 1101, odsa(2,1), odsa(2,2)
1101    format ('Paraxial ray  ',f21.11,'  ',f21.11,/)
        print 1102, lspher
1102    format ('Longitudinal spherical aberration:      ', f16.11)
        print 1103, mxlsph
1103    format ('    (Maximum permissible):              ', f16.11)
1109    format (42x, '   Acceptable', /)
1110    format (42x, '** Excessive **', /)
        if (abs(lspher) .le. mxlsph) then
           print 1109
        else
           print 1110
        end if

        print 1104, abosc
1104    format ('Offense against sine condition (coma):  ', f16.11)
        print 1105, mxcrom
1105    format ('    (Maximum permissible):              ', f16.11)
        if (abs(abosc) .le. mxcrom) then
           print 1109
        else
           print 1110
        end if

        print 1106, lchrom
1106    format ('Axial chromatic aberration:             ', f16.11)
        print 1107, maxosc
1107    format ('    (Maximum permissible):              ', f16.11)
        if (abs(lchrom) .le. maxosc) then
           print 1109
        else
           print 1110
        end if

        end

c    Perform ray trace in specific spectral line

        subroutine tline(nline, rayh)
        integer nline
        double precision rayh

c       Spectral line indices
        parameter (LineA = 1, LineB = 2, LineC = 3, LineD = 4,
     1             LineE = 5, LineF = 6, LineGp = 7, LineH = 8)

c       Fields in the design array
        parameter (Dcurad = 1, Dindex = 2, Ddisper = 3, Dethick = 4)

        double precision rheight, tlod, tlsa
        double precision findex, indref
        double precision radcur
        integer nsurf, axiali

        double precision slines(8)
        double precision design(4, 4)

        common /dat/ rheight, tlod, tlsa,
     1               findex, indref,
     2               radcur,
     3               axiali

c       The  test case used in this program is the design for a 4 inch
c       f/12 achromatic telescope objective used as the example in Wyld's
c       classic work on ray tracing by hand, given in Amateur Telescope
c       Making, Volume 3 (Volume 2 in the 1996 reprint edition).
        data nsurf / 4 /
        data design /
     1        27.05d0, 1.5137d0, 63.6d0, .52d0,
     2       -16.68d0, 1d0, 0d0, .138d0,
     3       -16.68d0, 1.6164d0, 36.7d0, .38d0,
     4       -78.1d0, 1d0, 0d0, 0d0 /

c       Wavelengths of standard spectral lines in Angstroms
c              (Not all are used in this program)
        data slines / 7621d0, 6869.955d0, 6562.816d0, 5895.944d0,
     1             5269.557d0, 4861.344d0, 4340.477d0, 3968.494d0 /

        tlod = 0
        rheight = rayh
        findex = 1

        do 201 i = 1, nsurf
           radcur = design(Dcurad, i)
           indref = design(Dindex, i)
           if (indref .gt. 1) then
              indref = indref + ((slines(LineD) - slines(nline)) /
     1            (slines(LineC) - slines(LineF))) *
     2            ((design(Dindex, i) - 1) / design(Ddisper, i))
           end if
           call tsurf
           findex = indref
           if (i .lt. nsurf) tlod = tlod - design(Dethick, i)
201     continue
        return
        end

c              Calculate passage through surface
c
c               If the variable axiali is true, the trace through the
c               surface will be done using the paraxial approximations.
c               Otherwise, the normal trigonometric trace will be done.
c
c                 iang    Incidence angle
c                 rang    Refraction angle
c                 iasin   Incidence angle sin
c                 rasin   Refraction angle sin

        subroutine tsurf

        parameter (Margin = 1, Parax = 2)

        double precision rheight, tlod, tlsa
        double precision findex, indref, iasin, rasin, iang
        double precision rang, radcur, sagita, oasa
        integer axiali

        common /dat/ rheight, tlod, tlsa,
     1               findex, indref,
     2               radcur,
     3               axiali

        if (axiali .ne. Margin) then

c           Paraxial ray
           if (radcur .ne. 0) then
              if (tlod .eq. 0) then
                 tlsa = 0
                 iasin = rheight / radcur
              else
                 iasin = ((tlod - radcur) / radcur) * tlsa
              end if
              rasin = (findex / indref) * iasin
              oasa = tlsa
              tlsa = tlsa + iasin - rasin
              if (tlod .ne. 0) rheight = tlod * oasa
              tlod = rheight / tlsa
              return
           end if
           tlod = tlod * (indref / findex)
           tlsa = tlsa * (findex / indref)
           return

        end if

c       Marginal ray
        if (radcur .ne. 0) then
           if (tlod .eq. 0) then
              tlsa = 0
              iasin = rheight / radcur
           else
              iasin = ((tlod - radcur) / radcur) * sin(tlsa)
           end if
           iang = asin(iasin)

           rasin = (findex / indref) * iasin
           oasa = tlsa
           tlsa = tlsa + iang - asin(rasin)
           sagita = 2 * radcur * ((sin((oasa + iang) / 2)) ** 2)
           tlod = ((radcur * sin(oasa + iang)) *
     1            (1 / tan(tlsa))) + sagita
           return
        end if
        rang = - asin((findex / indref) * sin(tlsa))
        tlod = tlod * ((indref * cos(-rang)) / (findex * cos(tlsa)))
        tlsa = -rang
        return
        end
