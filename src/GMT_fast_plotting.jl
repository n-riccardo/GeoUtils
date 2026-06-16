############################################################
### Useful functions for ternary plots
############################################################

function spherical_area(A, B, C)
    num = dot(A, cross(B, C))
    den = 1 + dot(A,B) + dot(B,C) + dot(A,C)
    return 2 * atan(num/den)
end

function spherical_barycentric(P, A, B, C)
    ΩABC = spherical_area(A, B, C)

    λA = spherical_area(P, B, C) / ΩABC
    λB = spherical_area(P, C, A) / ΩABC
    λC = spherical_area(P, A, B) / ΩABC

    return λA, λB, λC
end

function rgb_to_hex(colors)
    "#" .* [
        join(string.(round.(Int, 255 .* c), base=16, pad=2), "")
        for c in eachrow(colors)
    ]
end
function map_ternary_color_octant(lon,lat,cA,cB,cC)

    # longitude
    θ=deg2rad.(lon)
    # colatitude:
    φ=deg2rad.(90 .- lat)

    # Octant vertices (spherical triangle)
    A = normalize([1.0, 0.0, 0.0])
    B = normalize([0.0, 1.0, 0.0])
    C = normalize([0.0, 0.0, 1.0])

    Cmap = zeros(length(θ),3)
    λAvec = zeros(length(θ))
    λBvec = zeros(length(θ))
    λCvec = zeros(length(θ))

    for i in eachindex(θ)
        th = θ[i]
        ph = φ[i]

        # Cartesian coorinates
        x = cos(th)*sin(ph)
        y = sin(th)*sin(ph)
        z = cos(ph)

        P = [x, y, z]

        # Barycentric coordinates
        λA, λB, λC = spherical_barycentric(P, A, B, C)
        λAvec[i]=λA
        λBvec[i]=λB
        λCvec[i]=λC

        # Color interpolation
        color = λA .* cA .+ λB .* cB .+ λC .* cC
        Cmap[i,:] =  color
    end

    return Cmap, λAvec, λBvec, λCvec 

end
function map_ternary_color_octant(X,Y,Z,cA,cB,cC)
    
    # Octant vertices (spherical triangle)
    A = normalize([1.0, 0.0, 0.0])
    B = normalize([0.0, 1.0, 0.0])
    C = normalize([0.0, 0.0, 1.0])

    Cmap = zeros(length(X),3)
    λAvec = zeros(length(X))
    λBvec = zeros(length(X))
    λCvec = zeros(length(X))

    for i in eachindex(X)

        x = X[i]
        y = Y[i]
        z = Z[i]
    
        P = [x, y, z]
    
        # Barycentric coordinates
        λA, λB, λC = spherical_barycentric(P, A, B, C)
        λAvec[i]=λA
        λBvec[i]=λB
        λCvec[i]=λC

        # Color interpolation
        color = λA .* cA .+ λB .* cB .+ λC .* cC
        Cmap[i,:] =  color
    end
    
    return Cmap, λAvec, λBvec, λCvec 

end

"""
plot\\_vector\\_map_S:

*Author:* 
Riccardo Nucci (riccardo.nucci4@unibo.it)

*Description:*
Plot a vector velocity field on a map using GMT.jl. My favorite options are set as default.

*Required arguments:*
- range\\_plot::Vector{<:Real}
- scale\\_plot::Real
- velo\\_data1::Matrix{<:Real}: format is lon-lat-ve-vn-se-sn-corr

*Optional arguments:*
- legend\\_data::Matrix{<:Real}=[0.0 0.0 0.0 0.0]
- leg\\_offset::Vector{<:Real} = [0.0, 0.0]
- leg\\_font::String = "11p,Helvetica-Bold,black"
- projection\\_s::String = "M10c"
- frame\\_s::String = "a5f1"
- color::String ="blue"
- arrow::String = "0.2c+a45+p0.025c+e+n30/0.005+g"
- pen\\_coast::String = "0.01c,black"
- pen\\_velo::String = "0.012c,black"
- CI::String = "0.95"
- Are\\_you\\_Overwriting::Bool = false
- Sr::Bool = false
- transparency::String="50"

*Notes:*
- Use Are\\_you\\_Overwriting=false (default is true) to set basemap and coast
- Use Sr=true to plot the velocity field with the Sr GMT option for error ellipses

"""
function plot_vector_map_S(range_plot::Vector{<:Real}, scale_plot::Real, velo_data1::Matrix{<:Real}; legend_data::Matrix{<:Real}=[0.0 0.0 0.0 0.0], 
    leg_offset::Vector{<:Real} = [0.0, 0.0], leg_font::String = "11p,Helvetica-Bold,black", projection_s::String = "M15c", 
    frame_s::String = "af", color::String ="blue", arrow::String = "0.2c+a45+e+n30/0.005+g", 
    pen_coast::String = "0.01c,black", pen_velo::String = "0.012c,black", CI::String = "0.95", Are_you_Overwriting::Bool = false, Sr::Bool = false,
    transparency::String="50")

    if(!Are_you_Overwriting)

        gmtset(MAP_FRAME_TYPE="plain", PROJ_LENGTH_UNIT="c")
        basemap(region=range_plot, projection=projection_s, frame=frame_s)
        coast!(coast=pen_coast, area=1000)

    end

    scale_s = string(scale_plot) * "c"
	
    velo!(legend_data, pen=pen_velo,
        fill_wedges=color * "@"*transparency, outlines=true, Se=scale_s * "/" * CI,
        arrow=arrow * color)

    text!(x=legend_data[1] + leg_offset[1], y=legend_data[2] + leg_offset[2], text=string(legend_data[3])*" mm/yr", font=leg_font)

    if(Sr)
        velo!(velo_data1, pen=pen_velo,
            fill_wedges=color * "@"*transparency, outlines=true, Sr=scale_s * "/" * CI,
            arrow=arrow * color)
    else
        velo!(velo_data1, pen=pen_velo,
        fill_wedges=color * "@"*transparency, outlines=true, Se=scale_s * "/" * CI,
        arrow=arrow * color)
    end

end

"""
plot\\_fancy\\_scalar\\_field:

*Author:* 
Riccardo Nucci (riccardo.nucci4@unibo.it)

*Description:*
Plot a scalar field, e.g. topography, on a map using GMT.jl. My favorite options are set as default.

*Required arguments:*
- netcdf\\_path::String
- range\\_plot::Vector{<:Real}
- range\\_cpt::Vector{<:Real}

*Optional arguments:*

*BASEMAP and COAST arguments (active only if Are\\_you\\_Overwriting::Bool = true)*
- projection\\_s::String = "M15c"
- frame\\_s::String = "af"
- pen\\_coast::String = "0.01c,black"
- cpt_string::String="turbo"

*COLORBAR arguments (active only if colorbar\\_bool::Bool = true)*
- font\\_colorbar::String = "11p,Helvetica-Bold,black"
- position\\_colorbar::String="JMR+w10c/0.5c+v+o1c/0c"
- axis\\_colorbar::String="af"
- title\\_colorbar="z"

*Optionally set the color of water (superimposed to the grd image)*
- color\\_water::String=""

*SCALE (active only if scale\\_bool::Bool = true)*
- position\\_scale::Vector{<:Real}=[12,40]
- width\\_scale::Real=50

*SHADE (default is my_shade::Bool = true)*

*Notes:*
- Use Are\\_you\\_Overwriting=false (default is true) to set basemap and coast
"""
function plot_fancy_scalar_field(netcdf_path::String, range_plot::Vector{<:Real}, range_cpt::Vector{<:Real};
    Are_you_Overwriting::Bool = false, projection_s::String = "M15c", frame_s::String = "af", pen_coast::String = "0.01c,black",
    cpt_string::String="turbo",
    colorbar_bool::Bool=true, font_colorbar::String = "11p,Helvetica-Bold,black",  position_colorbar::String="JMR+w10c/0.5c+v+o1.5c/0c+e", axis_colorbar::String="af", title_colorbar="z",
    color_water::String="",
    scale_bool::Bool=true, position_scale::Vector{<:Real}=[12,40], width_scale::Real=50, my_shade::Bool=true, interp_::String="n")

    if(!Are_you_Overwriting)
        gmtset(MAP_FRAME_TYPE="plain", PROJ_LENGTH_UNIT="c")
        basemap(region=range_plot, projection=projection_s, frame=frame_s)
    end
    
    cpt = makecpt(cmap=cpt_string, range=string(range_cpt[1]) * "/" * string(range_cpt[2]))

    if(my_shade)
        grdimage!(netcdf_path, cmap=cpt, shade=my_shade)
    else
        grdimage!(netcdf_path, cmap=cpt, interp=interp_,dpi=:i )
    end

    if(colorbar_bool)
        gmtset(FONT_ANNOT_PRIMARY=font_colorbar)
        colorbar!(position=position_colorbar, frame=("x"*axis_colorbar*"+l"*title_colorbar))
    end
    
    if(!Are_you_Overwriting)
        coast!(coast=pen_coast, area=1000)
    end

    if(!isempty(color_water))   
        coast!(water=color_water)
        coast!(water=color_water)
        coast!(water=color_water)
    end

    if(scale_bool)
        basemap!(L="g"*string(position_scale[1])*"/"*string(position_scale[2])*"+w"*string(width_scale))
    end
end

function plot_fancy_scatterplot(x::Vector{<:Real}, y::Vector{<:Real}, z::Vector{<:Real}, range_plot::Vector{<:Real}, ptSzp::Vector{<:Real};
    Are_you_Overwriting::Bool = false, projection_s::String = "M15c", frame_s::String = "af", pen_coast::String = "0.01c,black",
    cpt_string::String="turbo",
    colorbar_bool::Bool=true, font_colorbar::String = "11p,Helvetica-Bold,black",  position_colorbar::String="JMR+w10c/0.5c+v+o1c/0c", axis_colorbar::String="af", title_colorbar="z",
    color_water::String="",
    scale_bool::Bool=true, position_scale::Vector{<:Real}=[12,40], width_scale::Real=50, shade::Bool=true)

    if(!Are_you_Overwriting)
        gmtset(MAP_FRAME_TYPE="plain", PROJ_LENGTH_UNIT="c")
        basemap(region=range_plot, projection=projection_s, frame=frame_s)
    end
    
    scatter!(x,y,markersize=ptSzp,C=cpt_string, zcolor=z)

    if(colorbar_bool)
        gmtset(FONT_ANNOT_PRIMARY=font_colorbar)
        colorbar!(position=position_colorbar, frame=("x"*axis_colorbar*"+l"*title_colorbar))
    end
    
    if(!Are_you_Overwriting)
        coast!(coast=pen_coast, area=1000)
    end

    if(!isempty(color_water))   
        coast!(water=color_water)
        coast!(water=color_water)
        coast!(water=color_water)
    end

    if(scale_bool)
        basemap!(L="g"*string(position_scale[1])*"/"*string(position_scale[2])*"+w"*string(width_scale))
    end
    
end

"""
matrix_strain must have:
1,2: longitude, latitude, of station (-: option interchanges order)

3: eps1, the most extensional eigenvalue of strain tensor, with extension taken positive.

4: eps2, the most compressional eigenvalue of strain tensor, with extension taken positive.

5: azimuth of eps2 in degrees CW from North.
"""
function PlotFancyHorStrainTensor(range_plot::Vector{<:Real}, scale_plot::Real, matrix_strain::Matrix{<:Real}; legend_data::Matrix{<:Real}=[0.0 0.0 0.0 0.0 0.0], 
    leg_offset::Vector{<:Real} = [0.0, 0.0], leg_font::String = "11p,Helvetica-Bold,black", projection_s::String = "M15c", 
    frame_s::String = "af", colors::Vector{String} =["blue","red"], arrow::String = "0.2c+a45+p0.025c+e+n30/0.005+g", 
    pen_coast::String = "0.01c,black", pen_velo::String = "0.012c", Are_you_Overwriting::Bool = false,
    transparency::String="50",legend_text::String="str/yr")

    matrix_strainExtension=copy(matrix_strain)
    matrix_strainExtension[:,4]=matrix_strain[:,4].*0
    matrix_strainCompression=copy(matrix_strain)
    matrix_strainCompression[:,3]=matrix_strain[:,3].*0

    legend_dataExtension=copy(legend_data)
    legend_dataExtension[:,4]=legend_data[:,4].*0
    legend_dataCompression=copy(legend_data)
    legend_dataCompression[:,3]=legend_data[:,3].*0

    if(!Are_you_Overwriting)

        gmtset(MAP_FRAME_TYPE="plain", PROJ_LENGTH_UNIT="c")
        basemap(region=range_plot, projection=projection_s, frame=frame_s)
        coast!(coast=pen_coast, area=1000)

    end

    scale_s = string(scale_plot) * "c"
	
    velo!(legend_dataCompression, pen=pen_velo*","*colors[1],
    outlines=true, Sx=scale_s,
    arrow=arrow * colors[1])

    velo!(legend_dataExtension, pen=pen_velo*","*colors[2],
    outlines=true, Sx=scale_s,
    arrow=arrow * colors[2])

    text!(x=legend_data[1] + leg_offset[1], y=legend_data[2] + leg_offset[2], text=string(legend_data[3])*legend_text, font=leg_font)

    velo!(matrix_strainCompression, pen=pen_velo*","*colors[1],
    outlines=true, Sx=scale_s,
    arrow=arrow * colors[1])

    velo!(matrix_strainExtension, pen=pen_velo*","*colors[2],
    outlines=true, Sx=scale_s,
    arrow=arrow * colors[2])


end

"""
https://www.generic-mapping-tools.org/GMTjl_doc/tutorials/ISC/isc.jl/

"""

function currentFunctions()
    # Get all functions defined in the current module
    all_symbols = names(Main, all=true)
    functions = filter(x -> isdefined(Main, x) && typeof(getfield(Main, x)) <: Function, all_symbols)
    
    println("Available functions:")
    for fn in functions
        println("  - ", fn)
    end
end
