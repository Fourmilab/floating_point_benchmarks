--[[

    John Walker's Floating Point Benchmark, derived from...

        Marinchip Interactive Lens Design System
               John Walker   December 1980

                   by John Walker
               http://www.fourmilab.ch/

    This program may be used, distributed, and modified freely as
    long as the origin information is preserved.

    This is a complete optical design raytracing algorithm,
    stripped of its user interface and recast into Lua.
    It not only determines execution speed on an extremely
    floating point (including trig function) intensive
    real-world application, it checks accuracy on an algorithm
    that is exquisitely sensitive to errors. The performance of
    this program is typically far more sensitive to changes in
    the efficiency of the trigonometric library routines than the
    average floating point program.

    Implemented in July 2014 by John Walker.

]]

--  Number of iterations to run for the benchmark

    niter = 9585971
    niter = 86984367
niter = 100

--      Wavelengths of standard spectral lines in Angstroms
--              (Not all are used in this program)

    a_line  = 7621.0                -- A
    b_line  = 6869.955              -- B
    c_line  = 6562.816              -- C
    d_line  = 5895.944              -- D
    e_line  = 5269.557              -- E
    f_line  = 4861.344              -- F
    g_prime_line = 4340.477         -- G'
    h_line  = 3968.494              -- H

--[[    The  test case used in this program is the design for a 4 inch
        f/12 achromatic telescope objective used as the example in Wyld's
        classic work on ray tracing by hand, given in Amateur Telescope
        Making, Volume 3 (Volume 2 in the 1996 reprint edition).
]]

    wyldClearAperture = 4     -- Clear aperture of the following lens design

    wyldLens = {
                  { curRad =  27.05, index = 1.5137, disp = 63.6, edge = 0.52 },
                  { curRad = -16.68, index = 1.0,    disp = 0.0,  edge = 0.138 },
                  { curRad = -16.68, index = 1.6164, disp = 36.7, edge = 0.38 },
                  { curRad = -78.1,  index = 1.0,    disp = 0.0,  edge = 0.0 }
               }

    --  Define cot in terms of tan

    function cot(x)
        return 1 / math.tan(x)
    end

--[[    Calculate passage through surface

        If the variable paraxial is true, the trace through the
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
                                relevant if object_distance == 0

        axis_slope_angle        Angle incoming ray makes with axis
                                at intercept

        from_index              Refractive index of medium being left

        to_index                Refractive index of medium being
                                entered.

        The outputs are the following global variables:

        object_distance         Distance from vertex to object focus
                                after refraction.

        axis_slope_angle        Angle incoming ray makes with axis
                                at intercept after refraction.
]]

    function transit_surface()
        if paraxial then

            --  Paraxial ray

            if radius_of_curvature ~= 0 then
                if object_distance == 0 then
                    axis_slope_angle = 0
                    iang_sin = ray_height / radius_of_curvature
                else
                    iang_sin = ((object_distance -
                            radius_of_curvature) / radius_of_curvature) *
                            axis_slope_angle
                end
                rang_sin = (from_index / to_index) * iang_sin
                old_axis_slope_angle = axis_slope_angle
                axis_slope_angle = axis_slope_angle +
                    iang_sin - rang_sin
                if object_distance ~= 0 then
                    ray_height = object_distance * old_axis_slope_angle
                end
                object_distance = ray_height / axis_slope_angle
            else
                object_distance = object_distance * (to_index / from_index)
                axis_slope_angle = axis_slope_angle * (from_index / to_index)
            end

        else

            --  Marginal ray

            if radius_of_curvature ~= 0 then
                if object_distance == 0 then
                    axis_slope_angle = 0
                    iang_sin = ray_height / radius_of_curvature
                else
                    iang_sin = ((object_distance -
                        radius_of_curvature) / radius_of_curvature) *
                        math.sin(axis_slope_angle)
                end
                iang = math.asin(iang_sin)
                rang_sin = (from_index / to_index) * iang_sin
                old_axis_slope_angle = axis_slope_angle
                axis_slope_angle = axis_slope_angle +
                    iang - math.asin(rang_sin)
                sagitta = math.sin((old_axis_slope_angle + iang) / 2)
                sagitta = 2 * radius_of_curvature * sagitta * sagitta
                object_distance = ((radius_of_curvature *
                    math.sin(old_axis_slope_angle + iang)) *
                    cot(axis_slope_angle)) + sagitta
            else

                rang = -math.asin((from_index / to_index) *
                    math.sin(axis_slope_angle))
                object_distance = object_distance * ((to_index *
                    math.cos(-rang)) / (from_index *
                    math.cos(axis_slope_angle)))
                axis_slope_angle = -rang
            end

        end
    end

--[[    Perform ray trace for a given design for a specific
        spectral line and ray height.  The caller specifies the
        desired spectral line and ray height.  The global
        object distance is updated based upon tracing this
        ray.
]]

    function trace_line(line, ray_h)
        object_distance = 0
        ray_height = ray_h
        from_index = 1

        for i = 1, current_surfaces do
            radius_of_curvature = s[i].curRad
            to_index = s[i].index
            if to_index > 1 then
                to_index = to_index + ((d_line - line) /
                    (c_line - f_line)) *
                    ((s[i].index - 1) / s[i].disp)
            end
            transit_surface()
            from_index = to_index
            if i < current_surfaces then
                object_distance = object_distance - s[i].edge
            end
        end
    end

    s = wyldLens
    clear_aperture = wyldClearAperture
    current_surfaces = #s

    for iteration = 1, niter do

        --      Trace marginal ray in D light

        paraxial = false

        trace_line(d_line, clear_aperture / 2)
        od_d_marginal, sa_d_marginal = object_distance, axis_slope_angle

        --      Trace paraxial ray in D light

        paraxial = true

        trace_line(d_line, clear_aperture / 2)
        od_d_paraxial, sa_d_paraxial = object_distance, axis_slope_angle

        --      Trace marginal ray in C light

        paraxial = false
        trace_line(c_line, clear_aperture / 2)
        od_cline = object_distance

        --      Trace marginal ray in F light

        trace_line(f_line, clear_aperture / 2)
        od_fline = object_distance


        --[[     Compute aberrations of the design

          The longitudinal spherical aberration is just the
          difference between where the D line comes to focus
          for paraxial and marginal rays.
        ]]

        aberr_lspher = od_d_paraxial - od_d_marginal

        --[[ The offense against the sine condition is a measure
             of the degree of coma in the design.  We compute it
             as the lateral distance in the focal plane between
             where a paraxial ray and marginal ray in the D line
             come to focus.
        ]]

        aberr_osc = 1 - ((od_d_paraxial * sa_d_paraxial) /
                          (math.sin(sa_d_marginal) * od_d_marginal))

        --[[ The axial chromatic aberration is the distance between
             where marginal rays in the C and F lines come to focus.
        ]]

        aberr_lchrom = od_fline - od_cline

        --[[ Compute maximum acceptable values for each aberration

             Maximum longitudinal spherical aberration, which is
             also the maximum for axial chromatic aberration.  This
             is computed for the D line.
        ]]

        sin_dm_sa = math.sin(sa_d_marginal)
        max_lspher = 0.0000926 / (sin_dm_sa * sin_dm_sa)
        max_lchrom = max_lspher
        max_osc = 0.0025

    end

    --  Print evaluation of the design based upon the last ray trace

    io.write(string.format("   Marginal ray   %21.11f  %14.11f\n",
        od_d_marginal, sa_d_marginal))
    io.write(string.format("   Paraxial ray   %21.11f  %14.11f\n",
        od_d_paraxial, sa_d_paraxial))

    io.write(string.format("Longitudinal spherical aberration:      %16.11f\n",
        aberr_lspher))
    io.write(string.format("    (Maximum permissible):              %16.11f\n",
        max_lspher))

    io.write(string.format("Offense against sine condition (coma):  %16.11f\n",
        aberr_osc))
    io.write(string.format("    (Maximum permissible):              %16.11f\n",
        max_osc))

    io.write(string.format("Axial chromatic aberration:             %16.11f\n",
        aberr_lchrom))
    io.write(string.format("    (Maximum permissible):              %16.11f\n",
        max_lchrom))
