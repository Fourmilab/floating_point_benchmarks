begin

    comment

        John Walker's Floating Point Benchmark, derived from...

        Marinchip Interactive Lens Design System

        By John Walker
           http://www.fourmilab.ch/

        This program may be used, distributed, and modified freely as
        long as the origin information is preserved.

        This is a complete optical design raytracing algorithm,
        stripped of its user interface and recast into ALGOL 60.
        It not only determines execution speed on an extremely
        floating point (including trig function) intensive
        real-world application, it checks accuracy on an algorithm
        that is exquisitely sensitive to errors. The performance of
        this program is typically far more sensitive to changes in
        the efficiency of the trigonometric library routines than the
        average floating point program.

        Ported from the Algol 68 language implementation in May 2014
        by John Walker.
    ;

    comment
        Variables and constants global to all procedures.
    ;
    
    integer a line, b line, c line, d line,
            e line, f line, g prime line,
            h line;
    real array spectral line[1:8];
    
    integer current surfaces;
    real array test case[1:4, 1:4];
    integer curvature radius, index of refraction, dispersion,
        edge thickness;
    
    real clear aperture, aberr lspher, aberr osc, aberr lchrom,
         max lspher, max osc, max lchrom, radius of curvature,
         object distance, ray height, axis slope angle,
         from index, to index;
         
    integer number of iterations, iteration;
         
    integer paraxial, marginal ray, paraxial ray;
    
    comment
        The od sa array holds the result from the main trace of
        paraxial and marginal rays in D light.  The first subscript
        indicates whether the ray is marginal or paraxial and the
        second selects the object distance or axis slope angle.
    ;
    
    integer od field, sa field;
    real array od sa[1:2, 1:2];
    
    real od cline, od fline;
    
    comment
        The cot function defined in terms of sin and cos
    ;
    
    real procedure cot(x);
        value x;
        real x; 
    begin
        cot := 1 / (sin(x) / cos(x))
    end cot;
    
    comment
        The arc sin function defined in terms of arc tan
    ;
    
    real procedure arc sin(x);
        value x;
        real x;
    begin
        arc sin := 2 * arc tan(x / (1 + sqrt(1 - x^2)))
    end arc sin;
    
    comment
        The built-in procedure outreal() uses the equivalent of
        a C format "%.12g", which does not provide sufficient
        precision for the largest numbers we edit.  This
        brutal hack gets around this limitation by defining
        an outbigreal() procedure which uses the marst
        inline() mechanism to print its argument with a format
        of "%.11f" as used by other implementations of the
        benchmark.  This is only needed, and only used, when
        printing the object distance for the marginal and
        paraxial rays in D light.
        
        If you have trouble with this, just replace the body of
        this procedure with:
        
            outreal(channel, x)
            
        but then you'll have to manually verify these values are
        rounded correctly.
        
        Usually, I would avoid using a compiler-specific feature
        like this in a benchmark, but since ALGOL 60 has no
        I/O facilities defined in the standard, they differ
        among almost all compilers, so any working version of the
        program is necessarily non-portable without modifying the
        I/O.
    ;
    
    procedure outbigreal(channel, x);
        value channel, x;
        integer channel;
        real x;
    begin
        real fraction;
        
        fraction := x;
        inline("printf(\"%.11f\", dsa_1->fraction_10);")
    end outbigreal;
    
    comment
    
        Calculate passage through surface
  
        If the variable paraxial is paraxial ray, the trace through the
        surface will be done using the paraxial approximations.
        Otherwise, the normal trigonometric trace will be done.
  
        This procedure takes the following global inputs:
  
        radius of curvature     Radius of curvature of surface
                                being crossed.  If 0, surface is
                                plane.
  
        object distance         Distance of object focus from
                                lens vertex.  If 0, incoming
                                rays are parallel and
                                the following must be specified:
  
        ray height              Height of ray from axis.  Only
                                relevant if object distance = 0
  
        axis slope angle        Angle incoming ray makes with axis
                                at intercept
  
        from index              Refractive index of medium being left
  
        to index                Refractive index of medium being
                                entered.
  
        The outputs are the following global variables:
  
        object distance         Distance from vertex to object focus
                                after refraction.
  
        axis slope angle        Angle incoming ray makes with axis
                                at intercept after refraction.
    ;

    procedure transit surface;
    begin
        real iang, rang, iang sin, rang sin,
             old axis slope angle, sagitta;
             
        if paraxial = paraxial ray then begin
            if radius of curvature != 0 then begin
                if object distance = 0 then begin
                    axis slope angle := 0;
                    iang sin := ray height / radius of curvature
                end else
                    iang sin := ((object distance -
                        radius of curvature) / radius of curvature) *
                        axis slope angle;
                rang sin := (from index / to index) * iang sin;
                old axis slope angle := axis slope angle;
                axis slope angle := axis slope angle +
                    iang sin - rang sin;
                if object distance != 0 then
                    ray height := object distance * old axis slope angle;
                object distance := ray height / axis slope angle
            end else begin
                object distance := object distance * (to index / from index);
                axis slope angle := axis slope angle * (from index / to index)
            end
        end else begin
            if radius of curvature != 0 then begin
                if object distance = 0 then begin
                    axis slope angle := 0;
                    iang sin := ray height / radius of curvature
                end else
                    iang sin := ((object distance -
                        radius of curvature) / radius of curvature) *
                        sin(axis slope angle);
                iang := arc sin(iang sin);
                rang sin := (from index / to index) * iang sin;
                old axis slope angle := axis slope angle;
                axis slope angle := axis slope angle +
                    iang - arcsin(rang sin);
                sagitta := sin((old axis slope angle + iang) / 2);
                sagitta := 2 * radius of curvature * (sagitta ^ 2);
                object distance := ((radius of curvature *
                    sin(old axis slope angle + iang)) *
                    cot(axis slope angle)) + sagitta
            end else begin
                rang := -arcsin((from index / to index) *
                    sin(axis slope angle));
                object distance := object distance * ((to index *
                    cos(-rang)) / (from index *
                    cos(axis slope angle)));
                axis slope angle := -rang
           end

        end
    end transit surface;
    
    comment
        Perform ray trace for a given design for a specific
        spectral line and ray height.  The caller specifies the
        desired spectral line and ray height.  The global
        object distance is updated based upon tracing this
        ray.
    ;
    
    procedure trace line(line, ray h);
        value line, ray h;
        integer line;
        real ray h;
    begin
        integer i;
        
        object distance := 0;
        ray height := ray h;
        from index := 1;
        
        for i := 1 step 1 until current surfaces do begin
            radius of curvature := test case[i, curvature radius];
            to index := test case[i, index of refraction];
            if to index > 1 then
                to index := to index + ((spectral line[d line] -
                    spectral line[line]) /
                    (spectral line[c line] - spectral line[f line])) *
                    ((test case[i, index of refraction] - 1.0) /
                    test case[i, dispersion]);
            transit surface;
            from index := to index;
            if i < current surfaces then
                object distance := object distance -
                    test case[i, edge thickness]
        end
    end trace line;
    
    comment END GLOBAL DECLARATONS ;
    
    comment Spectral lines ;
        
    a line := 1; b line := 2; c line := 3; d line := 4;
        e line := 5; f line := 6; g prime line := 7;
        h line := 8;
    
    spectral line[a line]       := 7621.0; 
    spectral line[b line]       := 6869.955; 
    spectral line[c line]       := 6562.8160; 
    spectral line[d line]       := 5895.944; 
    spectral line[e line]       := 5269.557; 
    spectral line[f line]       := 4861.344; 
    spectral line[g prime line] := 4340.477; 
    spectral line[h line]       := 3968.494; 
    
    comment  Initialise the test case array
    
             The  test case used in this program is the design for a 4 inch
             f/12 achromatic telescope objective used as the example in Wyld's
             classic work on ray tracing by hand, given in Amateur Telescope
             Making, Volume 3 (Volume 2 in the 1996 reprint edition). 
    ;
    
    current surfaces := 4;
    
    curvature radius := 1; index of refraction := 2; dispersion := 3;
        edge thickness := 4;
        
    test case[1, curvature radius]    :=   27.05;
    test case[1, index of refraction] :=    1.5137;
    test case[1, dispersion]          :=   63.6;
    test case[1, edge thickness]      :=    0.52;
    
    test case[2, curvature radius]    :=  -16.68;
    test case[2, index of refraction] :=    1.0;
    test case[2, dispersion]          :=    0.0;
    test case[2, edge thickness]      :=    0.138;
    
    test case[3, curvature radius]    :=  -16.68;
    test case[3, index of refraction] :=    1.6164;
    test case[3, dispersion]          :=   36.7;
    test case[3, edge thickness]      :=    0.38;
    
    test case[4, curvature radius]    :=  -78.1;
    test case[4, index of refraction] :=    1.0;
    test case[4, dispersion]          :=    0.0;
    test case[4, edge thickness]      :=    0.0;
    
    marginal ray := 1; paraxial ray := 2;
    
    od field := 1; sa field := 2;
 
 comment       
    number of iterations:= 40263820;
    
    number of iterations:= 1;
    
    clear aperture := 4;
    
    for iteration := 1 step 1 until number of iterations do begin
    
        for paraxial:= marginal ray step 1 until paraxial ray do begin
            trace line(d line, clear aperture / 2);
            od sa[paraxial, od field] := object distance;
            od sa[paraxial, sa field] := axis slope angle
        end;

        comment Trace marginal ray in C ;

        paraxial := marginal ray;
        trace line(c line, clear aperture / 2);
        od cline := object distance;

        comment Trace marginal ray in F ;

        trace line(f line, clear aperture / 2);
        od fline := object distance;

        comment

            Compute aberrations of the design

            The longitudinal spherical aberration is just the
            difference between where the D line comes to focus
            for paraxial and marginal rays. ;

        aberr lspher := od sa[paraxial ray, od field] -
            od sa[marginal ray, od field];

        comment

            The offense against the sine condition is a measure
            of the degree of coma in the design.  We compute it
            as the lateral distance in the focal plane between
            where a paraxial ray and marginal ray in the D line
            come to focus. ;

        aberr osc := 1 - ((od sa[paraxial ray, od field] *
                           od sa[paraxial ray, sa field]) /
                           (sin(od sa[marginal ray, sa field]) *
                            od sa[marginal ray, od field]));

        comment

            The axial chromatic aberration is the distance between
            where marginal rays in the C and F lines come to focus. ;

        aberr lchrom := od fline - od cline;

        comment

            Compute maximum acceptable values for each aberration

            Maximum longitudinal spherical aberration, which is
            also the maximum for axial chromatic aberration.  This
            is computed for the D line. ;

        max lspher := 0.0000926 / sin(od sa[marginal ray, sa field]) ^ 2;
        max lchrom := max lspher;
        max osc := 0.0025
    end;
    
    comment Print the analysis of the ray trace ;
    
    outstring(1, "  Marginal ray    ");
    outbigreal(1, od sa[marginal ray, od field]);
    outstring(1, "   ");
    outreal(1, od sa[marginal ray, sa field]); outstring(1, "\n");
    
    outstring(1, "  Paraxial ray    ");
    outbigreal(1, od sa[paraxial ray, od field]);
    outstring(1, "   ");
    outreal(1, od sa[paraxial ray, sa field]); outstring(1, "\n");
    
    outstring(1, "Longitudinal spherical aberration: ");
    outreal(1, aberr lspher); outstring(1, "\n");
    outstring(1, "    (Maximum permissible): ");
    outreal(1, max lspher); outstring(1, "\n");
    
    outstring(1, "Offense against sine condition (coma): ");
    outreal(1, aberr osc); outstring(1, "\n");
    outstring(1, "    (Maximum permissible): ");
    outreal(1, max osc); outstring(1, "\n");
    
    outstring(1, "Axial chromatic aberration: ");
    outreal(1, aberr lchrom); outstring(1, "\n");
    outstring(1, "    (Maximum permissible): ");
    outreal(1, max lchrom); outstring(1, "\n")

end
