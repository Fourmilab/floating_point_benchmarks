%  John Walker's Floating Point Benchmark, derived from...
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
%   stripped of its user interface and recast into Erlang.
%   It not only determines execution speed on an extremely
%   floating point (including trig function) intensive
%   real-world application, it checks accuracy on an algorithm
%   that is exquisitely sensitive to errors. The performance of
%   this program is typically far more sensitive to changes in
%   the efficiency of the trigonometric library routines than the
%   average floating point program.
%
%   Implemented in June 2014 by John Walker.

    -module(fbench).

    -export([main/0]).

    -import(math, [sin/1, cos/1, tan/1, asin/1]).
    
    %       Wavelengths of standard spectral lines in Angstroms
    %               (Not all are used in this program)

    -define(A_line, 7621.0).            % A
    -define(B_line, 6869.955).          % B
    -define(C_line, 6562.816).          % C
    -define(D_line, 5895.944).          % D
    -define(E_line, 5269.557).          % E
    -define(F_line, 4861.344).          % F
    -define(G_prime_line, 4340.477).    % G'
    -define(H_line, 3968.494).          % H

    %   The cot function defined in terms of tan

    cot(X) ->
        1 / tan(X).

    %   The heavy lifting of the ray trace is performed by the
    %   transit_surface function.  It is called with a 7-tuple
    %   Trace Context consisting of the following items.  The
    %   first is an integer indicating whether the ray is
    %   marginal or paraxial and the balance are all real
    %   numbers.
    %
    %   Axial_incidence         This is either ?MarginalRay or
    %                           ?ParaxialRay depending on the
    %                           offset of the initial ray from
    %                           the axis.
    %
    %   Radius_of_curvature     Radius of curvature of surface
    %                           being crossed.  If 0, surface is
    %                           plane.
    %   
    %   Object_distance         Distance of object focus from
    %                           lens vertex.  If 0, incoming
    %                           rays are parallel and
    %                           the following must be specified:
    %   
    %   Ray_height              Height of ray from axis.  Only
    %                           relevant if object_distance == 0
    %   
    %   Axis_slope_angle        Angle incoming ray makes with axis
    %                           at intercept
    %   
    %   From_index              Refractive index of medium being left
    %   
    %   
    %   To_index                Refractive index of medium being
    %                           entered.
    %   
    %   The result of a call on transit_surface is a Trace Context
    %   with the following components modified to reflect the geometry
    %   of the ray as it exits the surface.
    %   
    %   Object_distance         Distance from vertex to object focus
    %                           after refraction.
    %   
    %   Ray_height              Height of ray from axis.
    %   
    %   Axis_slope_angle        Angle incoming ray makes with axis
    %                           at intercept after refraction.

    
    -define(MarginalRay, 0).
    -define(ParaxialRay, 1).
    
    %   We break the function down into four cases depending
    %   upon whether the ray is marginal or paraxial and
    %   whether the surface is flat or curved.
 
    %   Paraxial ray with flat radius of curvature
    
    transit_surface({Ray, Radius, Odist,
                    Rheight, Asa, FromIndex, ToIndex}) 
        when ((Ray == ?ParaxialRay) and (Radius == 0)) ->
        { Ray, 0,
          (Odist * (ToIndex / FromIndex)),
          Rheight,
          (Asa * (FromIndex / ToIndex)),
          FromIndex, ToIndex }
    ;
    
    %   Paraxial ray with curved surface
    
    transit_surface({Ray, Radius, Odist,
                    Rheight, Asa, FromIndex, ToIndex})
        when (Ray == ?ParaxialRay) ->
        
        Odz = Odist == 0,
        
        AsaPrime = if Odz -> 0; true -> Asa end,
        
        IangSin = if Odz -> Rheight / Radius;
                     true -> ((Odist - Radius) / Radius) * Asa
                  end,
                  
        RangSin = (FromIndex / ToIndex) * IangSin,
        
        AsaDoublePrime = AsaPrime + IangSin - RangSin,
        
        RayHeightPrime = if Odz -> Rheight;
                            true -> Odist * AsaPrime end,
                            
        ObjectDistancePrime = RayHeightPrime / AsaDoublePrime,
        
        { Ray, Radius,
          ObjectDistancePrime,
          RayHeightPrime,
          AsaDoublePrime,
          FromIndex, ToIndex }
    ;

    
    %   Marginal ray with flat radius of curvature
    
    transit_surface({Ray, Radius, Odist,
                    Rheight, Asa, FromIndex, ToIndex}) 
        when ((Ray == ?MarginalRay) and (Radius == 0)) ->
        
        Rang = -(asin((FromIndex / ToIndex))) *
                    sin(Asa),
                    
        { Ray, 0,
          (Odist * ((ToIndex *
                    cos((-Rang))) / (FromIndex *
                    cos(Asa)))),
          Rheight,
          -Rang,
          FromIndex, ToIndex }
    ;
    
    %   Marginal ray with curved surface
    
    transit_surface({Ray, Radius, Odist,
                    Rheight, Asa, FromIndex, ToIndex}) ->
        Odz = Odist == 0,
        
        AsaPrime = if Odz -> 0; true -> Asa end,

        IangSin = if Odz -> Rheight / Radius;
                     true -> ((Odist -
                                Radius) /
                                Radius) *
                                sin(Asa) end,
        Iang = asin(IangSin),
        RangSin = (FromIndex / ToIndex) * IangSin,
        AsaDoublePrime = AsaPrime + Iang - asin(RangSin),
        SinAsaIang = sin((AsaPrime + Iang) / 2),
        Sagitta = 2 * Radius * SinAsaIang * SinAsaIang,
        Ray_heightPrime = if Odz -> Rheight;
                             true -> Odist * AsaPrime end,
        Object_distancePrime = ((Radius *
                               sin(AsaPrime + Iang)) *
                               cot(AsaDoublePrime)) + Sagitta,
                            
        { Ray, Radius,
          Object_distancePrime,
          Ray_heightPrime,
          AsaDoublePrime,
          FromIndex, ToIndex
        }
    .
    
    %   The trace_line function walks a given ray through the
    %   surfaces of a design, invoking transit_surface for each
    %   surface.  It is called with the lens design, the
    %   spectral line to be traced, and a Trace Context 7-tuple
    %   as described above for transit_surface.  It walks through
    %   the design list and uses transit_surface to follow the ray
    %   through the surfaces.
    
    trace_line([], _, Context) ->
        Context;
        
    trace_line([{Radius, Index, Dispersion, Edge} | TDesign], SpectralLine, Context) ->
        {AxialIncidence, _, ObjectDistance,
         RayHeight, AxisSlopeAngle, FromIndex, _} = Context,
        ToIndexPrime = if (Index > 1) ->
            Index +
                ((?D_line - SpectralLine) /
                 (?C_line - ?F_line)) * ((Index - 1) / Dispersion);
        true -> Index end,
        
        ContextPrime = transit_surface({AxialIncidence, Radius,
            ObjectDistance, RayHeight, AxisSlopeAngle,
            FromIndex, ToIndexPrime}),
        
        {AxialIncidenceP, RadiusOfCurvatureP, ObjectDistanceP,
         RayHeightP, AxisSlopeAngleP, _, ToIndexP} = ContextPrime,
        trace_line(TDesign, SpectralLine,
            {AxialIncidenceP, RadiusOfCurvatureP, ObjectDistanceP - Edge,
             RayHeightP, AxisSlopeAngleP, ToIndexP, 0}).
        
%       The trace_lens function is a little bit of syntactic sugar
%       to simplify invoking trace_line.  It creates an initial
%       Trace Context for trace_line which specifies the axial_incidence
%       and ray_height (computed from the clear aperture of the
%       design), invokes trace_line with the supplied design and
%       spectral line, then returns a tuple containing the
%       object_distance and axis_slope_angle resulting from the ray
%       trace.          

    trace_lens(Design, ClearAperture, SpectralLine, AxialIncidence) ->
        Context = trace_line(Design, SpectralLine,
            {AxialIncidence, 0, 0, ClearAperture / 2, 0, 1, 0}),
        {_, _, ObjectDistance,
         _, AxisSlopeAngle, _, _} = Context,
        {ObjectDistance, AxisSlopeAngle}
    .

%       The evaluate_design function performs a ray trace on a given design
%       with a specified clear aperture and returns a design evaluation
%       which includes the results for the D line and calculation of
%       spherical aberration, coma, and chromatic aberration, along with
%       the conventional acceptable upper bounds for these quantities.

    evaluate_design(Design, ClearAperture) ->
    
        %   Trace the rays upon which the evaluation will be based
        {Dm_od, Dm_sa} = trace_lens(Design, ClearAperture, ?D_line, ?MarginalRay),
        {Dp_od, Dp_sa} = trace_lens(Design, ClearAperture, ?D_line, ?ParaxialRay),
        {Cm_od, _} = trace_lens(Design, ClearAperture, ?C_line, ?MarginalRay),
        {Fm_od, _} = trace_lens(Design, ClearAperture, ?F_line, ?MarginalRay),

        %   Compute aberrations of the design
        
        %   The longitudinal spherical aberration is just the
        %   difference between where the D line comes to focus
        %   for paraxial and marginal rays.
        AberrLspher = Dp_od - Dm_od,
        
        %   The offense against the sine condition is a measure
        %   of the degree of coma in the design.  We compute it
        %   as the lateral distance in the focal plane between
        %   where a paraxial ray and marginal ray in the D line
        %   come to focus.
        AberrOSC = 1 - (Dp_od * Dp_sa) / (sin(Dm_sa) * Dm_od),
        
        %   The axial chromatic aberration is the distance between
        %   where marginal rays in the C and F lines come to focus.
        AberrLchrom = Fm_od - Cm_od,
        
        %   Compute maximum acceptable values for each aberration

        %   Maximum longitudinal spherical aberration, which is
        %   also the maximum for axial chromatic aberration.  This
        %   is computed for the D line.
        Sin_Dm_sa = sin(Dm_sa),
        MaxLspher = 0.0000926 / (Sin_Dm_sa * Sin_Dm_sa),
        
        {
            Dm_od, Dm_sa,           % D marginal ray
            Dp_od, Dp_sa,           % D paraxial ray
            AberrLspher,            % Longitudinal spherical aberration
            AberrOSC,               % Offense against the sine condiion (coma)
            AberrLchrom,            % Axial chromatic aberration
            MaxLspher,              % Maximum longitudinal spherical aberration
            0.0025,                 % Maximum offense against the sine condition
            MaxLspher               % Maximum axial chromatic aberration
                                    %   (NOT an error: see above)
        }       
    .
    
%       The run_benchmark function runs the raytrace and analysis the
%       specified number of times.  The DesignEvaluation from the last
%       run is returned.

    run_benchmark(1, Design, ClearAperture) ->
        evaluate_design(Design, ClearAperture)
    ;

    run_benchmark(N, Design, ClearAperture) ->
        evaluate_design(Design, ClearAperture),
        run_benchmark(N - 1, Design, ClearAperture)
    .
    
%       The evaluation_report function takes a DesignEvaluation and returns a
%       list of strings containing a primate-readable version of the evaluation.
%       For the purposes of the benchmark, this function is not timed; it serves
%       only to create textual results which can be compared against those from
%       the reference implementation.

    evaluation_report({Dm_od, Dm_sa, Dp_od, Dp_sa,
                       AberrLspher, AberrOSC, AberrLchrom,
                       MaxLspher, MaxOSC, MaxLchrom}) ->
        io:format("   Marginal ray        ~16.11f  ~14.11f~n", [Dm_od, Dm_sa]),
        io:format("   Paraxial ray        ~16.11f  ~14.11f~n", [Dp_od, Dp_sa]),
        io:format("Longitudinal spherical aberration:      ~16.11f~n", [AberrLspher]),
        io:format("    (Maximum permissible):              ~16.11f~n", [MaxLspher]),
        io:format("Offense against sine condition (coma):  ~16.11f~n", [AberrOSC]),
        io:format("    (Maximum permissible):              ~16.11f~n", [MaxOSC]),
        io:format("Axial chromatic aberration:             ~16.11f~n", [AberrLchrom]),
        io:format("    (Maximum permissible):              ~16.11f~n", [MaxLchrom])
    .
    
%   Get the current time in seconds and a fraction.

    time_seconds() ->
        {Tm, Ts, Tu} = now(),
        (Tm * 1.0E6) + Ts + (Tu * 1.0E-6)
    .
        
main() ->

    %   How many iterations to run.  For archival timing, adjust for
    %   a run time of around five  minutes.
%    Iterations = 21386789,
%    Iterations = 53152290,
    Iterations = 10000,
        
%   The  test case used in this program is the design for a 4 inch
%   f/12 achromatic telescope objective used as the example in Wyld's
%   classic work on ray tracing by hand, given in Amateur Telescope
%   Making, Volume 3 (Volume 2 in the 1996 reprint edition).

    
    WyldClearAperture = 4,
    
    WyldLens = [ %    CurRad   Index    Disp   Edge
                   {   27.05,  1.5137,  63.6,  0.52  },
                   {  -16.68,  1.0,      0.0,  0.138 },
                   {  -16.68,  1.6164,  36.7,  0.38  },
                   {  -78.1,   1.0,      0.0,  0.0   }
               ],
               
    
%    Evaluation = evaluate_design(WyldLens, WyldClearAperture),

    StartT = time_seconds(),
    Evaluation = run_benchmark(Iterations, WyldLens, WyldClearAperture),
    EndT = time_seconds(),
    evaluation_report(Evaluation),
    
    io:format("Run time for ~B iterations is ~.6f seconds.~n", [Iterations, EndT - StartT])
.


