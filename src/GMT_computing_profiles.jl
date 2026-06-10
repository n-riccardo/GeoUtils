#using GMT
#using LibGEOS
#using CSV
#using DataFrames
#using DelimitedFiles

const global where_my_functions_are = @__DIR__

include(where_my_functions_are*"/RhumbLinesCalculations.jl")

"""
compute\\_GNSS\\_profiles:

*Last Update:*
06-02-2025

*Author:* 
Riccardo Nucci (riccardo.nucci9@gmail.com)

*Description:*
This function compute the parallel and orthogonal velocities of GNSS stations along a profile (works with rhumblines).
Outputs are ready for GMT plotting.

*Required arguments:*
- cross\\_profile\\_half\\_width::Real  # half width of the profile in km
- start\\_end\\_points::Vector{<:Real}=[0.0, 0.0, 0.0, 0.0], # format is lon1,lat1,lon2,lat2
- velo\\_data1::Matrix{<:Real}=Array{Float64}(undef, 0, 0) # matrix in format lon,lat,ve,vn,vu,se,sn,su

*Optional arguments:*
- outer\\_output\\_folder::String # where intermediate results will be saved
- string\\_id::String             # a string to identify the profile
- string\\_id\\_GNSS::String=""   # a string to identify the GNSS velocity dataset
- num\\_points_track::Int64=500   # number of points along the track on which crossprofiles will be erected
- StationNames::Vector{String}    # vector containing the names of the stations

*Notes:*
- The function work in spherical Earth approximation
- Correlation between GNSS velocity components is ignored

"""
function compute_GNSS_profiles(cross_profile_half_width::Real, start_end_points::Vector{<:Real}, velo_data1::Matrix{<:Real}; outer_output_folder::String="", string_id::String="", string_id_GNSS::String="", num_points_track::Int64=500, StationNames::Vector{String}=String[])

    if isempty(velo_data1)
        print("\n-> Nothing to compute: please specify velo_data1::Matrix{<:Real} in the format lon-lat-ve-vn-vu-se-sn-su \n")
        return
    end

    output_folder=""

    if(outer_output_folder!="")

        output_folder = outer_output_folder * "/" * "Profile_"*string_id * "/GNSS_"*string_id_GNSS*"/"
        print("\n-> Output files will be saved in:" * output_folder * "\n")
        
        if (isdir(output_folder))
            rm(output_folder, recursive=true)
            print("\n-> Previous folder " * output_folder * " has been deleted\n")
        end

        mkpath(output_folder)

    end
        
    total_distance, azimuth_, dataOnTrack, distances, lon_rect, lat_rect = compute_tracks(cross_profile_half_width,start_end_points,num_points_track)

    # DataFrames relative to the SwathProfile
    Info=DataFrame(Lon_Start=[start_end_points[1]], Lat_Start=[start_end_points[2]], Lon_End=[start_end_points[3]], Lat_End=[start_end_points[4]], Total_Distance=[total_distance], Azimuth=[azimuth_], Width=[cross_profile_half_width*2])
    LonLatCentralTrack=DataFrame(Distances=distances, Lon=dataOnTrack[:,1], Lat=dataOnTrack[:,2])
    LonLatRect=DataFrame(Lon=lon_rect, Lat=lat_rect)

    if(output_folder!="")
        CSV.write(output_folder*"Info.csv", Info, writeheader=true)
        CSV.write(output_folder*"LonLatCentralTrack.csv", LonLatCentralTrack, writeheader=true)
        CSV.write(output_folder*"LonLatSwathProfile.csv", LonLatRect, writeheader=true)
    end

    # Fasten the computation, first select stations using a crude lon-lat rectangle:
    lon_min=minimum(lon_rect);
    lon_max=maximum(lon_rect);
    lat_min=minimum(lat_rect);
    lat_max=maximum(lat_rect);

    velo_data=copy(velo_data1);
    
    cond=((velo_data[:,1] .>= lon_min) .& (velo_data[:,1] .<= lon_max)) .& ((velo_data[:,2] .>= lat_min) .& (velo_data[:,2] .<= lat_max));
    
    velo_data=velo_data[cond,:];
    if(!isempty(StationNames))
        StationNames=StationNames[cond]
    end

    # Keep stations within the polygon:
    lon_lat_rect=[lon_rect lat_rect]
    polycontour=[lon_lat_rect[k,:] for k in 1:size(lon_lat_rect,1)];
    SwathPolygon=LibGEOS.Polygon([polycontour]);
    withinPoly=Bool[]
    for i=1:length(velo_data[:,1])
        onePoint=LibGEOS.Point(velo_data[i,1],velo_data[i,2]);
        push!(withinPoly,LibGEOS.intersects(onePoint,SwathPolygon));
    end

    velo_data=velo_data[withinPoly,:];
    if(!isempty(StationNames))
        StationNames=StationNames[withinPoly]
    end

    # Now compute the distance of each station from the beginning of the SwathProfile and compute parallel and orthogonal velocities
    proj_points_lon=Float64[];
    proj_points_lat=Float64[];
    disances_velo=Float64[];
    
    for i=1:length(velo_data[:,1])
        proj_point_lon,proj_point_lat=rhxrh(velo_data[i,1], velo_data[i,2], azimuth_, start_end_points[1], start_end_points[2], azimuth_-90)
        push!(proj_points_lon,proj_point_lon)
        push!(proj_points_lat,proj_point_lat)
        push!(disances_velo,rhumb_distance([velo_data[i,1], velo_data[i,2]], [proj_point_lon, proj_point_lat]))
    end
    
    dim_data_r = length(velo_data[:,1])
    velocity_parallel_data = zeros(dim_data_r)
    sigma_parallel_data = zeros(dim_data_r)
    velocity_orthogonal_data = zeros(dim_data_r)
    sigma_orthogonal_data = zeros(dim_data_r)

    # to avoid confusion, explicitly define the variables
    ve_data=velo_data[:,3];
    vn_data=velo_data[:,4];
    vu_data=velo_data[:,5];
    se_data=velo_data[:,6];
    sn_data=velo_data[:,7];
    su_data=velo_data[:,8];

    for i in 1:dim_data_r
        velocity_parallel_data[i] = ve_data[i] * sind(azimuth_) + vn_data[i] * cosd(azimuth_)
        velocity_orthogonal_data[i] = ve_data[i] * cosd(azimuth_) - vn_data[i] * sind(azimuth_) #positive on the right of profile direction
        sigma_parallel_data[i] = sqrt(sind(azimuth_)^2 * se_data[i]^2 + cosd(azimuth_)^2 * sn_data[i]^2)
        sigma_orthogonal_data[i] = sqrt(cosd(azimuth_)^2 * se_data[i]^2 + sind(azimuth_)^2 * sn_data[i]^2)
    end

    # DataFrames relative to the velocities
    if(isempty(StationNames))
        StationNames=string.(1:length(velo_data[:,1]))
    end
    Velocities=DataFrame(Station=StationNames, Lon=velo_data[:,1], Lat=velo_data[:,2], Distance=disances_velo,  Ve=ve_data, Vn=vn_data, Vu=vu_data, Se=se_data, Sn=sn_data, Su=su_data, Vel_Parallel=velocity_parallel_data, Sigma_Parallel=sigma_parallel_data, Vel_Orthogonal=velocity_orthogonal_data, Sigma_Orthogonal=sigma_orthogonal_data, Lon_Proj=proj_points_lon, Lat_Proj=proj_points_lat)

    # Sort the dataframe based on Distance:
    sort!(Velocities, :Distance)

    if(output_folder!="")
        CSV.write(output_folder*"Velocities.csv", Velocities, writeheader=true)
    end

    return Info,LonLatCentralTrack,LonLatRect,Velocities

end
### - ###
### - ###
### - ###
### - ###
"""
compute\\_profile\\_from\\_grid:

*Last Update:*
06-02-2025

*Author:* 
Riccardo Nucci (riccardo.nucci9@gmail.com)

*Description:*
This function uses GMT.jl functions to compute the scalar profile and crossprofile from a grid. 

*Required arguments:*
- cross\\_profile\\_half\\_width::Real 
- start\\_end_points::Vector{<:Real}=[0.0, 0.0, 0.0, 0.0], 
- grid\\_directory::String #a string of the path to the grid file in netcdf format

*Optional arguments:*
- outer\\_output\\_folder::String # where intermediate results will be saved
- string\\_id::String             # a string to identify the profile
- string\\_id\\_grid::String=""   # a string to identify the grid
- num\\_points\\_track::Int64=500
- sampling\\_dist::Real=1         # distance in km between points along the crossprofile

*Notes:*
- The function work in spherical Earth approximation
- Crossprofiles are erected using great circle lines
"""
function compute_profile_from_grid(cross_profile_half_width::Real, start_end_points::Vector{<:Real}=[0.0, 0.0, 0.0, 0.0],
    grid_directory::String=""; outer_output_folder::String="", string_id::String="", string_id_grid::String="",num_points_track::Int64=500, sampling_dist::Real=1)

    if isempty(grid_directory)
        print("\n-> Nothing to compute: please specify grid_directory::String file in netcdf format \n")
        return
    end

    output_folder=""

    if(outer_output_folder!="")

        output_folder = outer_output_folder * "/" * "Profile_"*string_id * "/Grid_"*string_id_grid*"/"
        print("\n-> Output files will be saved in:" * output_folder * "\n")
        
        if (isdir(output_folder))
            rm(output_folder, recursive=true)
            print("\n-> Previous folder " * output_folder * " has been deleted\n")
        end

        mkpath(output_folder)

    end

    total_distance, azimuth_, dataOnTrack, distances, lon_rect, lat_rect = compute_tracks(cross_profile_half_width,start_end_points,num_points_track)

    # DataFrames relative to the SwathProfile
    Info=DataFrame(Lon_Start=[start_end_points[1]], Lat_Start=[start_end_points[2]], Lon_End=[start_end_points[3]], Lat_End=[start_end_points[4]], Total_Distance=[total_distance], Azimuth=[azimuth_], Width=[cross_profile_half_width*2])
    LonLatCentralTrack=DataFrame(Distances=distances, Lon=dataOnTrack[:,1], Lat=dataOnTrack[:,2])
    LonLatRect=DataFrame(Lon=lon_rect, Lat=lat_rect)

    if(output_folder!="")
        CSV.write(output_folder*"Info.csv", Info, writeheader=true)
        CSV.write(output_folder*"LonLatCentralTrack.csv", LonLatCentralTrack, writeheader=true)
        CSV.write(output_folder*"LonLatSwathProfile.csv", LonLatRect, writeheader=true)
    end

    Lon_Lat_points=hcat(LonLatCentralTrack.Lon, LonLatCentralTrack.Lat);

    # Problem: if the point is outside the grid, it will not be included in the track
    GRD_track = grdtrack(Lon_Lat_points, grid=grid_directory)

    Lon_track = GRD_track[:, 1]
    Lat_track = GRD_track[:, 2]
    Values_track = GRD_track[:, 3]

    # Correction: look at the equal values of Lon-Lat
    Values_with_NaN = zeros(length(Lon_Lat_points[:, 1]))

    for i in 1:length(Lon_Lat_points[:, 1])

        my_index = findall(((Lon_track .== Lon_Lat_points[i, 1]) .& (Lat_track .== Lon_Lat_points[i, 2])))

        if (length(my_index) > 1)
            print("-> ERROR: multiple matches should not happen\n")
        elseif (length(my_index) == 1)
            Values_with_NaN[i] = Values_track[my_index[1]]
        else
            Values_with_NaN[i] = NaN
        end

    end

    # Save the results relative to the central track
    dfCentralTrack = DataFrame(
    Distances = LonLatCentralTrack.Distances,
    Longitude = Lon_Lat_points[:, 1],
    Latitude = Lon_Lat_points[:, 2],
    Values = Values_with_NaN
    )

    width_profile=cross_profile_half_width*2;
    cross_profile = grdtrack(Lon_Lat_points, grid=grid_directory, crossprofile="$(width_profile)k/$(sampling_dist)k")
    cross_profile_dataframe = DataFrame(cross_profile)

    dfCrossProfile=DataFrame(Longitude = Float64[], Latitude = Float64[], Values = Float64[])

    for i in 1:length(Lon_Lat_points[:,1])

        lon_temp = cross_profile_dataframe[i, 1]
        lat_temp = cross_profile_dataframe[i, 2]
        values_temp = cross_profile_dataframe[i, 5] #checked on GMT documentation

        for j in 1:length(lon_temp)
            push!(dfCrossProfile, hcat(lon_temp[j], lat_temp[j], values_temp[j]))
        end

        push!(dfCrossProfile, hcat(NaN,NaN,NaN))

    end

    if(output_folder!="")
        CSV.write(output_folder*"/ProfileCentralTrack.csv", dfCentralTrack, writeheader=true)
        CSV.write(output_folder*"/CrossProfileValues.csv", dfCrossProfile, writeheader=false)
    end

    return Info, LonLatCentralTrack, LonLatRect, dfCentralTrack, dfCrossProfile
    
end
### - ###
### - ###
### - ###
### - ###
"""
    compute\\_parorth\\_grid:

    *Author:*
    Riccardo Nucci (riccardo.nucci9@gmail.com)

    *Description:*
    It computes the parallel and orthogonal velocities for a given azimuth and return them

"""
function compute_parorth_grid(ve_model::Matrix{<:Real}, vn_model::Matrix{<:Real}, azimuth_::Real)


    # Velocity calculations
    velocity_parallel_model = ve_model .* sind.(azimuth_) .+ vn_model .* cosd.(azimuth_)
    velocity_orthogonal_model = ve_model .* cosd.(azimuth_) .- vn_model .* sind.(azimuth_)

    return velocity_parallel_model, velocity_orthogonal_model

end
### - ###
### - ###
### - ###
### - ###
"""
compute\\_seismicity\\_profiles:

*Last Update:*
06-02-2025

*Author:* 
Riccardo Nucci (riccardo.nucci9@gmail.com)

*Required arguments:*
- cross\\_profile\\_half\\_width::Real  # half width of the profile in km
- start\\_end\\_points::Vector{<:Real}=[0.0, 0.0, 0.0, 0.0], # format is lon1,lat1,lon2,lat2
- seismicity\\_data1::Matrix{<:Real}=Array{Float64}(undef, 0, 0) # matrix in format lon, lat, depth, Magnitude

*Optional arguments:*
- outer\\_output\\_folder::String # where intermediate results will be saved
- string\\_id::String             # a string to identify the profile
- string\\_id\\_seism::String=""   # a string to identify the seismicity dataset
- num\\_points_track::Int64=500   # number of points along the track on which crossprofiles will be erected

*Notes:*
- The function work in spherical Earth approximation
"""
function compute_seismicity_profiles(cross_profile_half_width::Real, start_end_points::Vector{<:Real}=[0.0, 0.0, 0.0, 0.0], seismicity_data1::Matrix{<:Real}=Array{Float64}(undef, 0, 0); outer_output_folder::String="", string_id::String="", string_id_seism::String="",num_points_track::Int64=500)

    if isempty(seismicity_data1)
        print("\n-> Nothing to compute: please specify seismicity::Matrix{<:Real} in the format lon-lat-depth-magnitude \n")
        return
    end

    output_folder=""

    if(outer_output_folder!="")

        output_folder = outer_output_folder * "/" * "Profile_"*string_id * "/Seism_"*string_id_seism*"/"
        print("\n-> Output files will be saved in:" * output_folder * "\n")
        
        if (isdir(output_folder))
            rm(output_folder, recursive=true)
            print("\n-> Previous folder " * output_folder * " has been deleted\n")
        end

        mkpath(output_folder)

    end

    total_distance, azimuth_, dataOnTrack, distances, lon_rect, lat_rect = compute_tracks(cross_profile_half_width,start_end_points,num_points_track)

    # DataFrames relative to the SwathProfile
    Info=DataFrame(Lon_Start=[start_end_points[1]], Lat_Start=[start_end_points[2]], Lon_End=[start_end_points[3]], Lat_End=[start_end_points[4]], Total_Distance=[total_distance], Azimuth=[azimuth_], Width=[cross_profile_half_width*2])
    LonLatCentralTrack=DataFrame(Distances=distances, Lon=dataOnTrack[:,1], Lat=dataOnTrack[:,2])
    LonLatRect=DataFrame(Lon=lon_rect, Lat=lat_rect)

    if(output_folder!="")
        CSV.write(output_folder*"Info.csv", Info, writeheader=true)
        CSV.write(output_folder*"LonLatCentralTrack.csv", LonLatCentralTrack, writeheader=true)
        CSV.write(output_folder*"LonLatSwathProfile.csv", LonLatRect, writeheader=true)
    end

    # Fasten the computation, first select events using a crude lon-lat rectangle:
    lon_min=minimum(lon_rect);
    lon_max=maximum(lon_rect);
    lat_min=minimum(lat_rect);
    lat_max=maximum(lat_rect);

    seismicity_data=copy(seismicity_data1);

    cond=((seismicity_data[:,1] .>= lon_min) .& (seismicity_data[:,1] .<= lon_max)) .& ((seismicity_data[:,2] .>= lat_min) .& (seismicity_data[:,2] .<= lat_max));
    
    seismicity_data=seismicity_data[cond,:];

    # Keep events within the polygon:
    lon_lat_rect=[lon_rect lat_rect]
    polycontour=[lon_lat_rect[k,:] for k in 1:size(lon_lat_rect,1)];
    SwathPolygon=LibGEOS.Polygon([polycontour]);
    withinPoly=Bool[]
    for i=1:length(seismicity_data[:,1])
        onePoint=LibGEOS.Point(seismicity_data[i,1],seismicity_data[i,2]);
        push!(withinPoly,LibGEOS.intersects(onePoint,SwathPolygon));
    end

    seismicity_data=seismicity_data[withinPoly,:];

    # Now compute the distance of each event from the beginning of the SwathProfile
    proj_points_lon=Float64[];
    proj_points_lat=Float64[];
    disances_seism=Float64[];
    
    for i=1:length(seismicity_data[:,1])
        proj_point_lon,proj_point_lat=rhxrh(seismicity_data[i,1], seismicity_data[i,2], azimuth_, start_end_points[1], start_end_points[2], azimuth_-90)
        push!(proj_points_lon,proj_point_lon)
        push!(proj_points_lat,proj_point_lat)
        push!(disances_seism,rhumb_distance([seismicity_data[i,1], seismicity_data[i,2]], [proj_point_lon, proj_point_lat]))
    end

    # DataFrames relative to the seismicity
    Seismicity=DataFrame(Lon=seismicity_data[:,1], Lat=seismicity_data[:,2], Distance=disances_seism,  Depth=seismicity_data[:,3], Magnitude=seismicity_data[:,4], Lon_Proj=proj_points_lon, Lat_Proj=proj_points_lat)

    # Sort the dataframe based on Distance:
    sort!(Seismicity, :Distance)

    if(output_folder!="")
        CSV.write(output_folder*"Seismicity.csv", Seismicity, writeheader=true)
    end

    return Seismicity

end
### - ###
### - ###
### - ###
### - ###
"""
percentile\\_from\\_cross\\_sections:

*Author:* 
Riccardo Nucci (riccardo.nucci9@gmail.com)

*Description:*
Automatically compute the desired percentile in orth directions from CrossProfileValues.

*Required arguments:*
- my\\_percentile::Real # The percentile to compute (one value only)

*Optional arguments (specify one of the two):*
- path\\_to\\_file::String       # Path of CrossProfileValues.csv
- CrossProfileValues::DataFrame  # DataFrame containing the crossprofile values

*Notes:*
- NaN values are ignored in the computation of the percentile
"""
function percentile_from_cross_sections(my_percentile::Real;path_to_file::String="",CrossProfileValues::DataFrame=DataFrame())

    if(isempty(path_to_file) && isempty(CrossProfileValues))
        print("\n-> Nothing to compute: please specify path_to_file::String or CrossProfileValues::DataFrame \n")
        return 0
    end

    if(isempty(path_to_file))
        path_to_file=where_my_functions_are*"/../temp/temp.csv"
        CSV.write(path_to_file,CrossProfileValues,writeheader=false)
    end

    # Read the file into an array
    data = readdlm(path_to_file, ',', Float64)  

    # Initialize storage for blocks
    blocks = []
    current_block = []

    # Separate data into blocks
    for row in eachrow(data)
        if all(isnan.(row))  # Check for NaN row as block separator
            if !isempty(current_block)
                push!(blocks, vcat(current_block...))  # Save current block
                current_block = []  # Reset for the next block
            end
        else
            push!(current_block, row')  # Add row to the current block
        end
    end

    # Add the last block if it's non-empty
    if !isempty(current_block)
        push!(blocks, vcat(current_block...))
    end

    # Compute the desired percentile for the third column in each block
    desired_percentile = my_percentile  
    percentile_results = Float64[]

    for block in blocks
        if !isempty(block)
            third_column = block[:, 3]  # Extract third column
			valid_values=third_column[.!(isnan.(third_column))]
			 if isempty(valid_values)
				push!(percentile_results, NaN)  # If no valid values, add NaN
			else
				push!(percentile_results, quantile(valid_values, desired_percentile / 100))
			end
        end
    end

    return percentile_results
end
### - ###
### - ###
### - ###
### - ###
"""
compute\\_tracks:

*Author:* 
Riccardo Nucci (riccardo.nucci9@gmail.com)

*Description:*
Compute azimuth distances, track and rectangle of a profile. It works with rhumb lines.

*Notes:*
- num_points_track defines the sampling of the distances
"""
function compute_tracks(HalfWidth,start_end_points,num_points_track)

    # General info on the profile:
    total_distance=rhumb_distance(start_end_points[1:2], start_end_points[3:4])
    azimuth_=Eval_azimuth(start_end_points[1:2], start_end_points[3:4])

    # Lon-lat coordinates of the "num_points_track" points on the profile:
    incremental_distance=total_distance/(num_points_track+1)
    dataOnTrackGMT=GMT.sample1d([start_end_points[1] start_end_points[2];start_end_points[3] start_end_points[4]], resample="R+l", inc=string(incremental_distance)*"k");
    dataOnTrack=Matrix(dataOnTrackGMT)

    # Distance of these points from the starting point:
    distances=zeros(length(dataOnTrack[:,1]))
    for i=1:length(dataOnTrack[:,1])
        distances[i]=rhumb_distance([start_end_points[1],start_end_points[2]], [dataOnTrack[i,1],dataOnTrack[i,2]])
    end

    # Lon-lat coordinates of the rectangle representing the Swath profile:
    lat_rect=Float64[];
    lon_rect=Float64[];
    push!(lon_rect,dataOnTrack[1,1]);
    push!(lat_rect,dataOnTrack[1,2]);
    for i=1:length(distances)

        lon_temp=dataOnTrack[i,1];
        lat_temp=dataOnTrack[i,2];

        # print([lon_temp,lat_temp])

        lon_rect_temp,lat_rect_temp=destination([lon_temp,lat_temp],azimuth_-90, HalfWidth);

        push!(lon_rect,lon_rect_temp);
        push!(lat_rect,lat_rect_temp);

    end
    for i=1:length(distances)

        lon_temp=dataOnTrack[end-i+1,1];
        lat_temp=dataOnTrack[end-i+1,2];

        lon_rect_temp,lat_rect_temp=destination([lon_temp,lat_temp],azimuth_+90, HalfWidth);

        push!(lon_rect,lon_rect_temp);
        push!(lat_rect,lat_rect_temp);

    end
    push!(lon_rect,dataOnTrack[1,1]);
    push!(lat_rect,dataOnTrack[1,2]);

    return total_distance,azimuth_,dataOnTrack,distances,lon_rect,lat_rect
end
### - ###
### - ###
### - ###
### - ###
"""
SegmentIntersection:

*Author:* 
Riccardo Nucci (riccardo.nucci9@gmail.com)

*Description:*
Given a file (ar a 2 columns DataFrame) with fault traces (as lon-lat rows; NaN-NaN separation between faults), it returns the intersection between fault traces and the profile (central track only).

*Notes:*
- num_points_track defines the sampling of the distances
"""
function SegmentIntersection(half_width::Real,start_end_points::Vector{<:Real},segments_path::String;num_points_track::Int64=500,delimiter::Char=' ')

    # Find intersection between track and segments:
    segments_df = CSV.read(segments_path, DataFrame; header=false, ignorerepeated=true, delim=delimiter)
    segments_df=segments_df[:,[1,2]] #Keep only lon-lat
    rename!(segments_df, [:lon, :lat])

    return SegmentIntersection(half_width,start_end_points,segments_df;num_points_track=num_points_track)

end
### - ###
### - ###
### - ###
### - ###
function SegmentIntersection(half_width::Real,start_end_points::Vector{<:Real},segments_df::DataFrame;num_points_track::Int64=500)

    _, _, dataOnTrack, distances, lon_rect, lat_rect = compute_tracks(half_width,start_end_points,num_points_track)
    LonLatCentralTrack=DataFrame(Distances=distances, Lon=dataOnTrack[:,1], Lat=dataOnTrack[:,2])
    LonLatRect=DataFrame(Lon=lon_rect, Lat=lat_rect)

    _,_,_,_,LonsLine,LatsLine = Intersect_track_blocksegments(segments_df, Matrix{Float64}(LonLatCentralTrack[:,[2,3]]),Matrix{Float64}(LonLatRect));

    # Only to compute the distance from the beginning of the profile:
    TempMatrix=[LonsLine LatsLine zeros(length(LonsLine)) zeros(length(LonsLine)) zeros(length(LonsLine)) ones(length(LonsLine)) ones(length(LonsLine)) ones(length(LonsLine))]
    _,_,_,temp=compute_GNSS_profiles(half_width*1.5, start_end_points, TempMatrix);
    # Results:
    DistancesLine=temp.Distance
    LonsLine=temp.Lon
    LatsLine=temp.Lat

    return DistancesLine,LonsLine,LatsLine
end
### - ###
### - ###
### - ###
### - ###
"""
IsolinesIntersection:

*Author:* 
Riccardo Nucci (riccardo.nucci9@gmail.com)

*Description:*
Gives the profile of structures defined through isolines. The input format is rigid: three columns separated as lon-lat-depth.
-  NaN rows separate different fault STRUCTURES (not isolines)
-  different isolines referring to the same structure are recognized on the basis of the depth field, no separation

*Notes:*
- num_points_track defines the sampling of the distances
"""
function IsolinesIntersection(half_width::Real,start_end_points::Vector{<:Real},isolines_df::DataFrame;num_points_track::Int64=500)

    _, _, dataOnTrack, distances, lon_rect, lat_rect = compute_tracks(half_width,start_end_points,num_points_track)
    LonLatCentralTrack=DataFrame(Distances=distances, Lon=dataOnTrack[:,1], Lat=dataOnTrack[:,2])
    LonLatRect=DataFrame(Lon=lon_rect, Lat=lat_rect)
    
    # Strategy is: separate the dataframe into blocks on the basis of the NaN values
    # for each block separate according for the different values of depths
    df_faults=divide_faults(isolines_df)

    IntersectionPoints = Matrix{Float64}[]

    for fault in df_faults
        isolines=divide_isolines(fault)
        IntersectionPointFault = Array{Float64}(undef, 0, 4)  # 0 rows, 3 columns
        for isoline in isolines
            _,_,_,_,LonsLine,LatsLine = Intersect_track_blocksegments(isoline, Matrix{Float64}(LonLatCentralTrack[:,[2,3]]),Matrix{Float64}(LonLatRect));

            if(!isempty(LonsLine))
            TempMatrix=[LonsLine LatsLine zeros(length(LonsLine)) zeros(length(LonsLine)) zeros(length(LonsLine)) ones(length(LonsLine)) ones(length(LonsLine)) ones(length(LonsLine))]
            _,_,_,temp=compute_GNSS_profiles(half_width*1.5, start_end_points, TempMatrix);
            # Results:
            DistancesLine=temp.Distance

            DepthsLine=[isoline[1,3] for i=1:length(LonsLine)]

            IntersectionPointFault= vcat(IntersectionPointFault, [DistancesLine LonsLine LatsLine DepthsLine])
            end
        end
        if(!isempty(IntersectionPointFault))
            push!(IntersectionPoints, IntersectionPointFault)
        end
    end
    
    return IntersectionPoints

end
### - ###
### - ###
### - ###
### - ###
"""
Intersect\\_track\\_blocksegments:

*Author:* 
Riccardo Nucci (riccardo.nucci9@gmail.com)

*Description:*
Given a dataframe of fault traces (as lon-lat rows; NaN-NaN separation between faults), it returns the geometric objects and coordinates of their intersections with the profile (both central track and rectangle)
"""
function Intersect_track_blocksegments(segments_df_input,track_ll_matrix,rectangle_ll_matrix)

my_segment_df=segments_df_input[:,[1,2]] #Suppose lon-lat for the first two columns

rename!(my_segment_df, [:lon, :lat])

segments_df=split_into_segments(my_segment_df) #Fault traces as segments

# Create a nx4 matrix
rows = Float64[]
buf = Float64[]
for r in eachrow(segments_df)

    x, y = Float64(r.lon), Float64(r.lat)
    if isnan(x) || isnan(y)
        empty!(buf)
        continue
    end
    push!(buf, x, y)
    if length(buf) == 4
        append!(rows, buf)
        empty!(buf)
    end
end

segments = reshape(rows, 4, :)'

println("Track-Segments intersection: ", size(segments,1), " segments found.")

GMTPolyLine=Compute_GMTPolyLine(track_ll_matrix)
GMTPoly=Compute_GMTPoly(rectangle_ll_matrix)

# Intersection points:
LonsLine=Float64[]
LatsLine=Float64[]
LonsPoly=Float64[]
LatsPoly=Float64[]

for idx=1:size(segments,1)

    segment_temp=segments[idx,:]

    GMTSegment=Compute_GMTPolyLine(Matrix{Float64}([segment_temp[1] segment_temp[2]; segment_temp[3] segment_temp[4]]))

    # Use gmtspatial to find intersection (spherical approximation with great circles (-jg) using the authalic radius) 
    inter_ds = GMT.gmtspatial(GMTSegment, GMTPolyLine, I="e")

    if(!isempty(inter_ds))
        push!(LonsLine,float(inter_ds[1,1]))
        push!(LatsLine,float(inter_ds[1,2]))
    end

    inter_ds = GMT.gmtspatial(GMTSegment, GMTPoly, I="e") 
    if(!isempty(inter_ds))
        push!(LonsPoly,float(inter_ds[1,1]))
        push!(LatsPoly,float(inter_ds[1,2]))
    end
end
	
println("Track-Segments intersection: ", length(LonsLine), " intersections for track.")
println("Track-Segments intersection: ", length(LonsPoly), " intersections for rectangle.")

println(LonsLine)
println(LatsLine)
println(LonsPoly)
println(LatsPoly)


return GMTPolyLine,GMTPoly,LonsPoly,LatsPoly,LonsLine,LatsLine

end
### - ###
### - ###
### - ###
### - ###
"""
Given a DataFrame with fault traces (lon-lat or lon-lat-depth) separated by NaN, split_into_segments atomize each fault into individual segments and put all toghether into a same DataFrame NaN-separated.
"""
function split_into_segments(input_df)

    df=copy(input_df)

    out = similar(df, 0) #crea un dataframe vuoto
    
    if(isnan(df[1,1]) || isnan(df[1,2]))
        df=df[2:end,:]
    end
    if(isnan(df[end,1]) || isnan(df[end,2]))
        df=df[1:end-1,:]
    end
    
    my_line_index = 1
    n = nrow(df)
    ncol_=ncol(df)
    
    while my_line_index <= n
    
        next_line_index=my_line_index
    
        while next_line_index <= n && ((!isnan(df[next_line_index,1])) && (!isnan(df[next_line_index,2])))
            next_line_index = next_line_index + 1
        end
    
        next_line_index=next_line_index-1
    
        if next_line_index - my_line_index >= 1
            for indx in my_line_index:(next_line_index-1)
                push!(out, df[indx,:])
                push!(out, df[indx+1,:])
                push!(out, fill(NaN, ncol_))
            end
        end
    
        my_line_index=next_line_index+2
    
    end
    
    return out
    
end
### - ###
### - ###
### - ###
### - ###
"""
Given a DataFrame with fault strucures (lon-lat or lon-lat-depth) separated by NaN, divide_faults return a vector of DataFrames each relative to a single fault
"""
function divide_faults(input_df)

    df=copy(input_df)

    outs= DataFrame[]   # vettore vuoto di DataFrame
    
    if(isnan(df[1,1]) || isnan(df[1,2]))
        df=df[2:end,:]
    end
    if(isnan(df[end,1]) || isnan(df[end,2]))
        df=df[1:end-1,:]
    end
    
    my_line_index = 1
    n = nrow(df)
    ncol_=ncol(df)

    while my_line_index <= n
    
        next_line_index=my_line_index
    
        while next_line_index <= n && ((!isnan(df.lon[next_line_index])) && (!isnan(df.lat[next_line_index])))
            next_line_index = next_line_index + 1
        end
    
        next_line_index=next_line_index-1
    
        if next_line_index - my_line_index >= 1
            indx_interval=my_line_index:next_line_index
            push!(outs, df[indx_interval,:])

        end
    
        my_line_index=next_line_index+2
    
    end
    
    return outs
    
end
### - ###
### - ###
### - ###
### - ###
"""
Given a DataFrame with fault structure defined by isolines (lon-lat-depth) WITHOUT NaN, divide_isolines return a vector of DataFrames, each relative to a single isoline
"""
function divide_isolines(input_df)

    df=copy(input_df)

    outs = DataFrame[]  # vettore vuoto di DataFrame
    
    if(isnan(df[1,1]) || isnan(df[1,2]))
        df=df[2:end,:]
    end
    if(isnan(df[end,1]) || isnan(df[end,2]))
        df=df[1:end-1,:]
    end
    
    my_line_index = 1
    n = nrow(df)
    ncol_=ncol(df)

    my_depth_value=df[1,3]

    while my_line_index <= n
    
        next_line_index=my_line_index
    
        while next_line_index <= n && (df[next_line_index,3]==my_depth_value)
            next_line_index = next_line_index + 1
        end

        if(next_line_index!=(n+1))
            my_depth_value=df[next_line_index,3]
        end
        
        next_line_index=next_line_index-1
    
        if next_line_index - my_line_index >= 1
            indx_interval=my_line_index:next_line_index
            push!(outs, input_df[indx_interval,:])
        end
    
        my_line_index=next_line_index+1
    
    end
    
    return outs

end
### - ###
### - ###
### - ###
### - ###
function Compute_GMTPolyLine(LonLatMatrix) #GMT line

	poly_coords = copy(LonLatMatrix)

	poly_file_temp = mktemp()[1]
	open(poly_file_temp, "w") do io
    
		for i in 1:size(poly_coords,1)
			println(io, "$(poly_coords[i,1]) $(poly_coords[i,2])")
		end
		
	end
	
	GMTPolyLine=GMT.gmtspatial(poly_file_temp,F="l")# F="l" per linea non chiusa

	return GMTPolyLine
end
### - ###
### - ###
### - ###
### - ###
function Compute_GMTPoly(LonLatMatrix) #closed GMT polygon

	poly_coords = copy(LonLatMatrix)

    if ((poly_coords[1,1] != poly_coords[end,1]) || (poly_coords[1,2] != poly_coords[end,2])) #Close it it is not
        poly_coords = vcat(poly_coords, poly_coords[1, :])
    end
    
    poly_file_temp = mktemp()[1]
    open(poly_file_temp, "w") do io
    
        for i in 1:size(poly_coords,1)
            println(io, "$(poly_coords[i,1]) $(poly_coords[i,2])")
        end
        
    end
    
    GMTPoly=GMT.gmtspatial(poly_file_temp,F="")

	return GMTPoly
end
### - ###
### - ###
### - ###
### - ###
function Parallel_and_Meridian_Intersection(half_width,start_end_points,PorM;num_points_track=500)

    _, _, dataOnTrack, distances, _, _ = compute_tracks(half_width,start_end_points,num_points_track)
    LonLatCentralTrack=DataFrame(Distances=distances, Lon=dataOnTrack[:,1], Lat=dataOnTrack[:,2])
    
    GMTPolyLine=Compute_GMTPolyLine(Matrix{Float64}(LonLatCentralTrack[:,[2,3]]))
    
    LonsParallel,LatsParallel,LonsMeridians,LatsMeridians = Intersect_track_parallel_and_meridians(GMTPolyLine,1,1);
    
    if(PorM==0)
        TempMatrix=[LonsParallel LatsParallel zeros(length(LonsParallel)) zeros(length(LonsParallel)) zeros(length(LonsParallel)) ones(length(LonsParallel)) ones(length(LonsParallel)) ones(length(LonsParallel))]
    else
        TempMatrix=[LonsMeridians LatsMeridians zeros(length(LonsMeridians)) zeros(length(LonsMeridians)) zeros(length(LonsMeridians)) ones(length(LonsMeridians)) ones(length(LonsMeridians)) ones(length(LonsMeridians))]
    end
    ~,~,~,temp=compute_GNSS_profiles(half_width*1.5, start_end_points, TempMatrix);
    
    # Results:
    DistancesPorM=temp.Distance
    LonsPorM=temp.Lon
    LatsPorM=temp.Lat;
    
    return DistancesPorM,LonsPorM,LatsPorM
    
end
### - ###
### - ###
### - ###
### - ###
function Intersect_track_parallel_and_meridians(GMTPolyLine,step_lon,step_lat)

	meridians = [((λ, -90.0), (λ, 90.0)) for λ in -180.0:step_lon:180.0]
	
	lats = collect(-90.0:step_lat:90.0)
	valid_lats = filter(φ -> abs(φ) < 90.0, lats)

	parallels1 = [((0.0, φ), (180.0, φ)) for φ in valid_lats]
	parallels2 = [((180.0, φ), (360.0, φ)) for φ in valid_lats]
	parallels=vcat(parallels1,parallels2)
	
	LonsMeridians=Float64[]
	LatsMeridians=Float64[]

	for (k,(p1,p2)) in enumerate(meridians)
		# Scrivi il segmento k in file temporaneo
		meridian_file = mktemp()[1]
		open(meridian_file, "w") do io
			println(io, "$(p1[1]) $(p1[2])")
			println(io, "$(p2[1]) $(p2[2])")
		end

		GMTMeridian=GMT.gmtspatial(meridian_file,F="l")

		inter_ds = GMT.gmtspatial(GMTMeridian, GMTPolyLine, I="e") 
		if(!isempty(inter_ds))
			#println("intersection found (Line)")
			#println(typeof(float(inter_ds[1])))
			push!(LonsMeridians,float(inter_ds[1,1]))
			push!(LatsMeridians,float(inter_ds[1,2]))
		end
    
	end
	
	LonsParallel=Float64[]
	LatsParallel=Float64[]

	for (k,(p1,p2)) in enumerate(parallels)
		# Scrivi il segmento k in file temporaneo
		parallel_file = mktemp()[1]
		open(parallel_file, "w") do io
			println(io, "$(p1[1]) $(p1[2])")
			println(io, "$(p2[1]) $(p2[2])")
		end

		GMTParallel=GMT.gmtspatial(parallel_file,F="l")

		inter_ds = GMT.gmtspatial(GMTParallel, GMTPolyLine, I="e") 
		if(!isempty(inter_ds))
			#println("intersection found (Line)")
			#println(typeof(float(inter_ds[1])))
			push!(LonsParallel,float(inter_ds[1,1]))
			push!(LatsParallel,float(inter_ds[1,2]))
		end
    
	end
	
	return LonsParallel,LatsParallel,LonsMeridians,LatsMeridians

end
### - ###
### - ###
### - ###
### - ###
function Grid_Vel_Profile_Percentiles(half_width,start_end_points,NetCDFModel,ParallelOrOrth;num_points_track=500,up_perc_vel=97.5,low_perc_vel=2.5)

    total_distance, azimuth_, _, _, _, _ = compute_tracks(half_width,start_end_points,num_points_track)
    Info=DataFrame(Lon_Start=[start_end_points[1]], Lat_Start=[start_end_points[2]], Lon_End=[start_end_points[3]], Lat_End=[start_end_points[4]], Total_Distance=[total_distance], Azimuth=[azimuth_], Width=[half_width*2])

    # Read the velo model
    vDat=Dataset(NetCDFModel)
    lon=Vector(vDat["lon"][:])
    lat=Vector(vDat["lat"][:])
    ve=Matrix(vDat["Ve"])
    vn=Matrix(vDat["Vn"])
    vp, vo = compute_parorth_grid(ve, vn, Info.Azimuth[1])

    if(ParallelOrOrth) #true is parallel
        veloModel=vp
    else
        veloModel=vo
    end

    temp_file = mktemp()[1]
    
    write_to_netcdf(temp_file,lon,lat,veloModel)

    _, _, _, ctv, cpv=compute_profile_from_grid(half_width,start_end_points,temp_file);
    
    GVel=percentile_from_cross_sections(up_perc_vel,CrossProfileValues=cpv);
    LVel=percentile_from_cross_sections(low_perc_vel,CrossProfileValues=cpv);
    MVel=percentile_from_cross_sections(50,CrossProfileValues=cpv);

    DistanceVel=ctv.Distances

    return DistanceVel,GVel,LVel,MVel
end