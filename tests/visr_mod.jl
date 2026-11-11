using DelaunayTriangulation
using LinearAlgebra
using Statistics

"""
Single station record with the same information that the Fortran code read
from the solution file.
"""
mutable struct Station
    lon::Float64
    lat::Float64
    x::Float64
    y::Float64
    vx::Float64
    sigmax::Float64
    vy::Float64
    sigmay::Float64
    area::Float64
end

"""
Prepared input bundle.
"""
struct Parameters
    lon_mean::Float64
    lat_mean::Float64
    cutoff_distance::Float64
    min_tau::Float64
    max_tau::Float64
    tau_step::Float64
    weight_threshold::Float64
end

"""
Result of one VISR solve at one location.
`status == 0` means success. status,Ux,Uy,[exx,exy,eyy],omega,covariance,chisq,total_weight,rtau,nselected)VisrSolution(status, 0.0, 0.0, [0.0,0.0,0.0], 0.0, zeros{Float64}(6,6), 0.0, 0.0, 0.0, 0)
"""
struct VisrSolution
    status::Int
    Ux::Float64
    Uy::Float64
    strain_rate::Vector{Float64}
    omega::Float64
    covariance::Matrix{Float64}
    chisq::Float64
    maximum_gap::Float64
    total_weight::Float64
    rtau::Float64
    nselected::Int
    design::Matrix{Float64}
end

function llxy(center_lon::Float64, center_lat::Float64,
    longitudes::AbstractVector{Float64}, latitudes::AbstractVector{Float64})

    DEG2RAD = pi / 180.0

    # Same geodetic-to-local projection as the original Fortran routine.
    flattening = 1.0 / 298.2572
    eccentricity2 = 2.0 * flattening - flattening^2
    center_lat_rad = center_lat * DEG2RAD
    center_lon_rad = center_lon * DEG2RAD
    radius = 6378.137 * sqrt(1.0 - eccentricity2) / (1.0 - eccentricity2 * sin(center_lat_rad)^2)

    t11 =  sin(center_lat_rad) * cos(center_lon_rad)
    t12 =  sin(center_lat_rad) * sin(center_lon_rad)
    t13 = -cos(center_lat_rad)
    t21 = -sin(center_lon_rad)
    t22 =  cos(center_lon_rad)

    local_x=zeros(length(longitudes))
    local_y=zeros(length(longitudes))

    for i in eachindex(longitudes)
        lat_rad = latitudes[i] * DEG2RAD
        lon_rad = longitudes[i] * DEG2RAD
        v1 = cos(lat_rad) * cos(lon_rad)
        v2 = cos(lat_rad) * sin(lon_rad)
        v3 = sin(lat_rad)
        local_x[i] = radius * (t21 * v1 + t22 * v2)
        local_y[i] = -radius * (t11 * v1 + t12 * v2 + t13 * v3)
    end

    return local_x,local_y
end


function solve_point_local(parameters::Parameters,x_val::Float64, y_val::Float64, stations::Vector{Station})

    tau = parameters.min_tau
    nstations=length(stations)
    status=0
    
    while tau <= parameters.max_tau

        status=0
        rtau = copy(tau)
        total_weight = 0.0

        indx_selected = Int64[]
        azimuths_selected = Float64[]
        voronoi_area_selected=Float64[]

        for i in 1:nstations

            dx = stations[i].x - x_val
            dy = stations[i].y  - y_val
            distance = hypot(dx, dy)

            if distance / rtau > parameters.cutoff_distance
                continue
            end

            az = atan(dx, dy)
            if az < 0
                az = az + 2π
            end
            az_deg = rad2deg(az)

            area = stations[i].area

            push!(indx_selected,i)
            push!(azimuths_selected,az_deg)
            push!(voronoi_area_selected,area)
        end

        nselected=length(indx_selected)

        if(nselected < 3) # in this case not enough stations to perform lsq computation
            status=1
            tau += parameters.tau_step
            continue
        end

        sort!(azimuths_selected)
        gap_selected = diff(azimuths_selected)
        push!(gap_selected, (azimuths_selected[1] + 360) - azimuths_selected[end])
        maximum_gap = maximum(gap_selected)

        if(maximum_gap > 180) # in this case coverage of stations around gridpoint is poor
            status=2
            tau += parameters.tau_step
            continue
        end

        mean_area = sum(voronoi_area_selected) / nselected

        weights = zeros(nselected)
        for i in 1:nselected
            idx = indx_selected[i]
            weights[i] = stations[idx].area / mean_area
        end

        design = zeros(2 * nselected, 6)
        rhs = zeros(2 * nselected)
        for i in 1:nselected
            idx = indx_selected[i]
            dx = stations[idx].x - x_val
            dy = stations[idx].y - y_val
            distance = hypot(dx, dy)
            weight = exp(-((distance / rtau)^2)) * weights[i]

            total_weight += 2.0 * weight #two components

            sx = stations[idx].sigmax
            sy = stations[idx].sigmay

            root_weight = sqrt(weight)
            
            design[2i - 1, 1] = root_weight / sx
            design[2i - 1, 3] = root_weight * dx / sx
            design[2i - 1, 4] = root_weight * dy / sx
            design[2i - 1, 6] = root_weight * dy / sx
            design[2i, 2] = root_weight / sy
            design[2i, 4] = root_weight * dx / sy
            design[2i, 5] = root_weight * dy / sy
            design[2i, 6] = -root_weight * dx / sy
            
            rhs[2i - 1] = root_weight * stations[idx].vx / sx
            rhs[2i] = root_weight * stations[idx].vy / sy
        end

        if total_weight < parameters.weight_threshold # in this case not yet above weight_threshold
            status=3
            tau += parameters.tau_step
            continue
        end

        solution = design \ rhs
        #solution = (design' * design) \ (design' * rhs)

        N = (design')* design
        b=(design')*rhs
        I6 = Matrix(I, 6, 6)

        covariance = N \ I6
        chisq = (sum(abs2, rhs)) - ((b')*covariance*b)

        Ux=solution[1]
        Uy=solution[2]
        exx=solution[3] .* (10^(3)) #nstr/yr
        exy=solution[4] .* (10^(3))
        eyy=solution[5] .* (10^(3))
        omega=solution[6] .* (10^(3))

        return VisrSolution(status,Ux,Uy,[exx,exy,eyy],omega,covariance,chisq,maximum_gap,total_weight,rtau,nselected,design)

    end


    # if it attives here it means that computation failed:
    void_corr=zeros(Float64, 6, 6)
    void_design=zeros(Float64, 6, 6)

    return VisrSolution(status, 0.0, 0.0, [0.0,0.0,0.0], 0.0,void_corr, 0.0, 0.0, 0.0, 0.0, 0,void_design)

end

function circular_area_fallback(x::AbstractVector{Float64}, y::AbstractVector{Float64}, ith::Int)

    sample_points::Int=6
    coefficient::Float64=2.0

    # Same geometric fallback used by the Fortran code for hull cells.
    nearest_distances = fill(1.0e6, sample_points)
    largest_previous = 0.0

    for sample in 1:sample_points
        nearest_distances[sample] = 1.0e6
        for i in eachindex(x)
            if(i == ith)
                continue
            end
            distance = hypot(x[ith] - x[i], y[ith] - y[i])

            if (distance < nearest_distances[sample]) && (distance > largest_previous)
                nearest_distances[sample] = distance
            end
        end
        largest_previous = nearest_distances[sample]
    end

    mean_distance = sum(@view nearest_distances[1:sample_points]) / sample_points
    radius = coefficient * mean_distance / 2.0
    return (pi * radius^2), radius

end

function compute_voronoi_areas(x,y)

    points = [x'; y']
    triangulation = DelaunayTriangulation.triangulate(points)
    tessellation = DelaunayTriangulation.voronoi(triangulation)
    areas = fill(-1.0, length(x))

    polygons_x=Float64[]
    polygons_y=Float64[]
    unbounded = Set(DelaunayTriangulation.get_unbounded_polygons(tessellation))

    for i in 1:num_polygons(tessellation) 
        if(i ∉ unbounded)
            vertices_indices = get_polygon(tessellation, i)
            for indx in vertices_indices
    
                my_point = get_polygon_point(tessellation, indx)
    
                for indx in vertices_indices

                    push!(polygons_x,my_point[1])
                    push!(polygons_y,my_point[2])

                end
            end 
                push!(polygons_x,NaN)
                push!(polygons_y,NaN)
        end
    end

    unbounded = Set(DelaunayTriangulation.get_unbounded_polygons(tessellation))
    for i in eachindex(areas)
        if i in unbounded
            continue
        end
        areas[i] = DelaunayTriangulation.get_area(tessellation, i)
    end

    indices_fallback=Int64[]
    radii=Float64[]
    for i in eachindex(areas)
        fallback_area, radius = circular_area_fallback(x, y, i)
        if (areas[i] == -1.0) || (areas[i] > (2.0 * fallback_area))
            areas[i] = fallback_area
            push!(indices_fallback,i)
            push!(radii,radius)
        end
    end

    return areas,indices_fallback,radii,polygons_x,polygons_y

end

function visr_mod(parameters::Parameters, stations::Vector{Station}, lon_min::Float64, lon_max::Float64,
    lat_min::Float64, lat_max::Float64, dlon::Float64, dlat::Float64)

    x_stations=Float64[];
    y_stations=Float64[];
    for i in 1:length(stations)
        x_temp, y_temp = llxy(parameters.lon_mean, parameters.lat_mean, [stations[i].lon], [stations[i].lat])
        push!(x_stations,x_temp[1])
        push!(y_stations,y_temp[1])
    end

    areas, _,_,_,_=compute_voronoi_areas(x_stations,y_stations)

    for i in 1:length(stations)
        stations[i].x=x_stations[i]
        stations[i].y=y_stations[i]
        stations[i].area=areas[i]
    end

    nlon = Int(floor((lon_max - lon_min) / dlon + 1.01))
    nlat = Int(floor((lat_max - lat_min) / dlat + 1.01))

    results=Vector{VisrSolution}()
    lons=Float64[]
    lats=Float64[]

    for j in 1:nlat
        lat_val = lat_min + dlat * (j - 1)
        for i in 1:nlon
            lon_val = lon_min + dlon * (i - 1)
            x_val, y_val = llxy(parameters.lon_mean, parameters.lat_mean, [lon_val], [lat_val])
            x_val=x_val[1]
            y_val=y_val[1]
            result = solve_point_local(parameters, x_val, y_val, stations)
            push!(results,result)
            push!(lons,lon_val)
            push!(lats,lat_val)
        end
    end

    return lons,lats,results

end