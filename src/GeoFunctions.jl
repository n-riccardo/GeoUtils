# Define the Earth's radius
# https://www.jpz.se/Html_filer/wgs_84.html
global const Earth_Radius = 6371.0087714 #Mean Radius (km)
global const Earth_SemimajorAxis = 6378.137 #km
global const Earth_SemiminorAxis = 6356.75231424518 #km
global const Earth_Eccentricity=0.081819190842622 #not squared


############################################################
### Rhumb Lines calculations (possibly they will be accurately implemented in julia packages)
### Only the Spherical case is implemented
############################################################
"""
Returns distance in km between two points along a rhumb line.
LonLat1: Start point, as Vector
LonLat2: End point, as Vector
Return the rhumbline distance in km
"""
function rhumb_distance(LonLat1::Vector{<:Real}, LonLat2::Vector{<:Real})

    lon_a = deg2rad(LonLat1[1])
    lat_a = deg2rad(LonLat1[2])

    lon_b = deg2rad(LonLat2[1])
    lat_b = deg2rad(LonLat2[2])
    
    Δλ = lon_b - lon_a
    ΔL = lat_b - lat_a

    ΔΣ = log((tan(pi / 4 + lat_b / 2))/(tan(pi / 4 + lat_a / 2)))
    
    #If the points lie on the same parallel:
    small_threshold=1e-11
    if abs(ΔΣ) > small_threshold
        q = ΔL / ΔΣ
    else
        q = cos(lat_a)
    end
    
    #Longitude difference must take the shortest path:
    if(abs(Δλ) > pi)
        if(Δλ > 0)
            Δλ=-(2 * pi - Δλ)
        else
            Δλ=2 * pi + Δλ
        end
    end
    
    dist = sqrt((ΔL^2) + (q^2) * (Δλ^2)) * Earth_Radius

    return dist
end
function rhumb_distance(lon1, lat1, lon2, lat2)

    LonLat1=[lon1,lat1]
    LonLat2=[lon2,lat2]
    return rhumb_distance(LonLat1, LonLat2)

end

"""
Returns azimuth between two points in degrees (azimuth is computed from North, positive clockwise).
LonLat1: Start point, as Vector
LonLat2: End point, as Vector
Return the azimuth in degrees
"""
function Eval_azimuth(LonLat1::Vector{<:Real}, LonLat2::Vector{<:Real})

    lon_a= deg2rad(LonLat1[1])
    lat_a = deg2rad(LonLat1[2])
    
    lon_b = deg2rad(LonLat2[1])
    lat_b = deg2rad(LonLat2[2])
    
    ΔΣ = log((tan(pi / 4 + lat_b / 2)) / (tan(pi / 4 + lat_a / 2)))
    Δλ = lon_b - lon_a
    
    if(abs(Δλ) > pi)
        if(Δλ > 0)
            Δλ=-(2 * pi - Δλ)
        else
            Δλ=2 * pi + Δλ
        end
    end
    
    return atan(Δλ, ΔΣ)*180/pi

end
function Eval_azimuth(lon1, lat1, lon2, lat2)

    LonLat1=[lon1,lat1]
    LonLat2=[lon2,lat2]
    return Eval_azimuth(LonLat1, LonLat2)
    
end

"""
Returns the point B, found at certain distance from point A and at constant azimuth.
LonLat1: Start point, as Vector
azimuth: The azimuth in degrees
distance: Distance in kilometers
Return point B coordinates as lon_b,lat_b
"""
function destination(LonLat1::Vector{<:Real}, azimuth::Real, distance::Real)

    lon_a = deg2rad(LonLat1[1])
    lat_a = deg2rad(LonLat1[2])

    theta = deg2rad(azimuth)

    d = distance
    delta = d / Earth_Radius
    ΔL = delta * cos(theta)
    lat_b = lat_a + ΔL
    ΔΣ = log((tan(pi / 4 + lat_b / 2)) / (tan(pi / 4 + lat_a / 2)))

    #if the points lie on the same parallel:
    if abs(ΔΣ) > 1e-11
        q = ΔL / ΔΣ
    else
        q = cos(lat_a)
    end
    
    Δλ = delta * sin(theta) / q
    lon_b = lon_a + Δλ
    
    if abs(lat_b) > pi / 2
        lat_b = if lat_b > 0
            pi - lat_b
        else
            -pi - lat_b
        end
    end
    
    lat_b = rad2deg(lat_b)
    lon_b = rad2deg(lon_b)
    lon_b = (540 + lon_b) % 360 - 180

    return [lon_b, lat_b]
end
function destination(lon1, lat1, azimuth, distance)

    LonLat1=[lon1,lat1]
    lon_lat_b= destination(LonLat1, azimuth, distance)
    return lon_lat_b[1], lon_lat_b[2]

end
"""
Given two points and two azimuths, this function returns the intersection point of the two rhumb lines.
long1: Longitude of the first point
lat1: Latitude of the first point
az1: Azimuth of the first rhumb line
long2: Longitude of the second point
lat2: Latitude of the second point
az2: Azimuth of the second rhumb line
Return the intersection point as long, lat
"""
function rhxrh(long1, lat1, az1, long2, lat2, az2)
    
    # Convert input angles to radians
    long1, lat1, az1 = deg2rad(long1), deg2rad(lat1), deg2rad(az1)
    long2, lat2, az2 = deg2rad(long2), deg2rad(lat2), deg2rad(az2)
    
    angle_from_east1_=pi/2-az1
    angle_from_east2_=pi/2-az2

    # put the angles in the range [-π, π]
    angle_from_east1 = mod(angle_from_east1_ + pi, 2*pi) - pi
    angle_from_east2 = mod(angle_from_east2_ + pi, 2*pi) - pi

    if((angle_from_east1==(pi/2)) || (angle_from_east1==(-pi/2)))
        slope1=Inf
    else
        slope1=tan(angle_from_east1)
    end

    if((angle_from_east2==(pi/2)) || (angle_from_east2==(-pi/2)))
        slope2=Inf
    else
        slope2=tan(angle_from_east2)
    end
    
    x1, y1 = PartialMercator(long1, lat1)
    x2, y2 = PartialMercator(long2, lat2)
    
    # Compute the intersections
    if(isinf(slope1) && !isinf(slope2))
        xhat=x1
        yhat=slope2*(xhat-x2)+y2
    elseif(!isinf(slope1) && isinf(slope2))
        xhat=x2
        yhat=slope1*(xhat-x1)+y1
    elseif(isinf(slope1) && isinf(slope2))
        xhat=NaN
        yhat=NaN
    else
        xhat = ((y2- y1) - ((slope2 * x2) - slope1 * x1)) / (slope1 - slope2)
        yhat = slope1 * (xhat - x1) + y1
    end
    
    newlong, newlat = PartialInverseMercator(xhat, yhat)
    
    newlong = mod(newlong + pi, 2*pi) - pi  # Wrap to [-π, π]
    newlongdeg, newlatdeg = rad2deg(newlong), rad2deg(newlat)
    
    return newlongdeg, newlatdeg
    
end

function PartialMercator(long, lat)
    x = long
    y = log.(tan.(pi/4 .+ lat/2))
    return x, y
end

function PartialInverseMercator(x, y)
    long = x
    lat = 2 .* atan(exp.(y)) .- pi/2
    return long, lat
end
function generate_path(lon1,lat1,lon2,lat2,incremental_distance)
    #calculation for loxodromes is spherical in GMT.jl, input space is adjusted!
    dataOnTrack=GMT.sample1d([lon1 lat1;lon2 lat2], resample="R+l", inc=string(incremental_distance)*"k"); 
    return dataOnTrack=Matrix(dataOnTrack)
end
############################################################
### Geodesic Lines calculations (using GeographicLib)
### The WGS84 case is used
############################################################
function geodesic_distance(lon1, lat1, lon2, lat2) #use WGS84

    _, _, dist_m, _=GeographicLib.inverse(GeographicLib.WGS84,lon1, lat1, lon2, lat2)
    dist_km=dist_m/1000;
    return dist_km

end

function haversine(LongLat1, LongLat2,  my_Earth_Radius)  #use the Spherical Earth approx.
    lon1, lat1 = LongLat1
    lon2, lat2 = LongLat2

    Δlat = deg2rad(lat2 - lat1)
    Δlon = deg2rad(lon2 - lon1)
    
    a = sin(Δlat / 2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(Δlon / 2)^2
    
    c = 2 * asin(min(1.0, sqrt(a))) 

    return my_Earth_Radius * c
end
function haversine(lon1, lat1, lon2, lat2) #already assume you want to use the const. Earth_Radius
    return haversine([lon1,lat1], [lon2, lat2],  Earth_Radius)
end

function geodesic_azimuth(lon1, lat1, lon2, lat2) #azimuth on the starting point of the great circle line

    azi, _, _, _=GeographicLib.inverse(GeographicLib.WGS84,lon1, lat1, lon2, lat2)
    return azi
end

function geodesic_destination(lon_1, lat_1, azimuth, dist_km)

    dist_m=dist_km*1000
    lon_2, lat_2, _, _, _ = GeographicLib.forward(GeographicLib.WGS84,lon_1, lat_1, azimuth, dist_m)
    return lon_2,lat_2

end

# The intersection between gedesics must be implemented

function geodesic_generate_path(lon1,lat1,lon2,lat2,incremental_distance)
    #calculation for geodesic (WGS84), input space is adjusted!
    dataOnTrack=GMT.sample1d([lon1 lat1;lon2 lat2], resample="R", inc=string(incremental_distance)*"k"); 
    return dataOnTrack=Matrix(dataOnTrack)
end

############################################################
### Geometries 
############################################################

function generate_circles(coords_and_radius; num_points=360)

    circles = Matrix{Float64}[]

    for row in eachrow(coords_and_radius)

        lon_center, lat_center, radius_km = row

        circle_points = zeros(num_points, 2)

        for (i, azimuth) in enumerate(range(0, 360, length=num_points))
            lon, lat = geodesic_destination(lon_center, lat_center, azimuth, radius_km)
            circle_points[i, 1] = lon
            circle_points[i, 2] = lat
        end

        push!(circles, circle_points)  # Store each circle as a tuple (lon_list, lat_list)
    end

    return circles
end
function generate_circles(lons, lats, radii; num_points=360)

    coords_and_radius=[lons lats radii]
    circles=generate_circles(coords_and_radius; num_points=num_points)
    return circles

end
function generate_circle(lon1, lat1, radius; num_points=360)

    coords_and_radius=[lon1 lat1 radius]
    circle=generate_circles(coords_and_radius; num_points=num_points)
    return circle[1]

end

function rectangleForGMT(region,llStep)

	lon1, lon2, lat1, lat2 = region

	# Lati
	lon_vec = lon1:llStep:lon2
	lat_vec = lat1:llStep:lat2

	# Costruzione del perimetro (in senso orario)
	lon = vcat(
		lon_vec,                      # lato basso (O → E)
		repeat([lon2], length(lat_vec)-1),  # lato destro (S → N)
		reverse(lon_vec)[2:end],      # lato alto (E → O)
		repeat([lon1], length(lat_vec)-1)   # lato sinistro (N → S)
	)

	lat = vcat(
		repeat([lat1], length(lon_vec)),     # lato basso
		lat_vec[2:end],                      # lato destro
		repeat([lat2], length(lon_vec)-1),   # lato alto
		reverse(lat_vec)[2:end]              # lato sinistro
	)

	# Chiudiamo il poligono
	regionLons = vcat(lon, lon[1])
	regionLats = vcat(lat, lat[1])
	
	return regionLons, regionLats

end

############################################################
### Grids
############################################################
function makeH3grid(order,lon1,lat1,nrings)

    base = latLngToCell(LatLng(deg2rad(lat1), deg2rad(lon1)), order)
    rings = gridDisk(base, nrings)

    # order 4 is approx 26 km
    # order 5 is approx 10 km

    my_cells=Matrix{Float64}[]

    for boundary in cellToBoundary.(rings)
        x=[]
        y=[]
        for geo in boundary
            push!(x, rad2deg(geo.lng))
            push!(y, rad2deg(geo.lat)) 
        end
        push!(x, rad2deg(boundary[1].lng)) #close the polygon
        push!(y, rad2deg(boundary[1].lat)) 

        push!(my_cells,[x y])
    end

    return my_cells

end
############################################################
### Utilities
############################################################

function mean_longitude(lon; w=ones(length(lon)))

    x = 0.0
    y = 0.0
    sw = sum(w)

    for (λ, wi) in zip(lon, w)
        r = deg2rad(λ)
        x += cos(r) * wi
        y += sin(r) * wi
    end

    mean_rad = atan(y/sw, x/sw)
    mean_deg = rad2deg(mean_rad)

    if(mean_deg == -180)
		mean_deg=180
	end

    return mean_deg
end
function mean_longitude(lon1, lon2; w1=1, w2=1)

    return mean_longitude([lon1,lon2], w=[w1,w2])
	
end

function wrap_longitude(lon::Real)

    lon_shifted = copy(lon)
	
    while lon_shifted <= -180
        lon_shifted += 360
    end
	
    while lon_shifted > 180
        lon_shifted -= 360
    end
	
    return lon_shifted == -180 ? 180.0 : lon_shifted

end

"""
    km2deg(distance_km::Float64) -> Float64

# Author:
Riccardo Nucci (riccardo.nucci9@gmail.com)

# Description:
Convert a distance from kilometers to degrees of latitude in a Spherical Earth.

"""
function km2deg(distance_km::Float64)::Float64

    degrees_per_km = 360.0 / (2 * π * Earth_Radius)

    return distance_km * degrees_per_km

end

"""
    km2deg_longitude(distance_km::Float64) -> Float64

# Author:
Riccardo Nucci (riccardo.nucci9@gmail.com)

# Description:
Convert a distance from kilometers to degrees of longitude in a Spherical Earth at a certain latitude.

"""
function km2deg_longitude(distance_km::Float64, latitude::Float64)::Float64
    
    lat_rad = deg2rad(latitude)
    
    # Length of one degree of longitude at the given latitude
    km_per_deg_longitude = Earth_Radius * cos(lat_rad) * 2 * π / 360
    
    # Convert distance in kilometers to degrees of longitude
    deg_longitude = distance_km / km_per_deg_longitude
    
    return deg_longitude

end

############################################################
### Euler Poles
############################################################

function euler_pole_from_rotation(omega; cov_omega=[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]) # input in deg/Myr

    omega_rad_Myr=deg2rad.(omega)

    omega_lon_rad=atan(omega_rad_Myr[2],omega_rad_Myr[1])
    omega_lat_rad=atan(omega_rad_Myr[3] / sqrt((omega_rad_Myr[1]^2)+(omega_rad_Myr[2]^2)) )
    magn_rad_Myr=sqrt((omega_rad_Myr[1]^2)+(omega_rad_Myr[2]^2)+(omega_rad_Myr[3]^2)) #rad/Myr
    
    wx=omega_rad_Myr[1]
    wy=omega_rad_Myr[2]
    wz=omega_rad_Myr[3]
    W_m=copy(magn_rad_Myr)
    h2 = wx^2 + wy^2

    factor=1/(W_m^2)

    JacobianMatrix=[ # lon,lat, mod
        -wy/h2 wx/h2 0.0;
        -factor*wx*wz/(sqrt(h2)) -factor*wy*wz/(sqrt(h2)) factor*sqrt(h2);
        wx/W_m wy/W_m wz/W_m
    ]

    cov_omega_rad=cov_omega .* ((pi/180)^2)
    covariance_result=JacobianMatrix*cov_omega_rad*(JacobianMatrix')

    omega_lon_deg=rad2deg.(omega_lon_rad)
    omega_lat_deg=rad2deg.(omega_lat_rad)
    magn_deg_Myr=rad2deg(magn_rad_Myr)
    covariance_result_deg=covariance_result .* ((180/pi)^2)

    sigma_lon = sqrt(covariance_result_deg[1,1])
    sigma_lat = sqrt(covariance_result_deg[2,2])
    sigma_magn = sqrt(covariance_result_deg[3,3])
    println("Euler Pole: \n lon: $omega_lon_deg pm $sigma_lon \n lat: $omega_lat_deg pm $sigma_lat \n rotation velocity: $magn_deg_Myr pm $sigma_magn deg/Myr \n")

    return omega_lon_deg, omega_lat_deg, magn_deg_Myr, covariance_result_deg

end


function rotation_from_euler_pole(omega_lon,omega_lat,omega_magn; cov_input)

    omega_lon_rad=deg2rad(omega_lon)
    omega_lat_rad=deg2rad(omega_lat)
    omega_magn_rad_Myr=deg2rad(omega_magn)


    omega_vec= omega_magn_rad_Myr .* [
        cos(omega_lat_rad)*cos(omega_lon_rad);
        cos(omega_lat_rad)*sin(omega_lon_rad);
        sin(omega_lat_rad)
    ]

    JacobianMatrix=[
        -omega_magn_rad_Myr*cos(omega_lat_rad)*sin(omega_lon_rad) -omega_magn_rad_Myr*sin(omega_lat_rad)*cos(omega_lon_rad) cos(omega_lat_rad)*cos(omega_lon_rad)
        omega_magn_rad_Myr*cos(omega_lat_rad)*cos(omega_lon_rad) -omega_magn_rad_Myr*sin(omega_lat_rad)*sin(omega_lon_rad) cos(omega_lat_rad)*sin(omega_lon_rad)
        0.0 omega_magn_rad_Myr*cos(omega_lat_rad) sin(omega_lat_rad)
    ]

    cov_input_rad=cov_input .* ((pi/180)^2)
    cov_output_rad=JacobianMatrix*cov_input_rad*(JacobianMatrix')

    omega_vec_deg=rad2deg.(omega_vec)
    cov_output_deg=cov_output_rad .* ((180/pi)^2)

    return omega_vec_deg, cov_output_deg # deg/Myr

end


function euler_velocity(lon, lat, omega; cov_input=[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])

    λ = deg2rad(lon)
    φ = deg2rad(lat)
    slat = sin(φ)
    clat = cos(φ)
    slon = sin(λ)
    clon = cos(λ)

    R_EARTH=copy(Earth_Radius)
    R_EARTH=R_EARTH*(10^6); #mm

    K_matrix = [
        -R_EARTH*clon*slat   -R_EARTH*slon*slat   R_EARTH*clat;
        R_EARTH*slon       -R_EARTH*clon         0.0
    ]

    omega_rad_yr=deg2rad.(omega) .* (10^(-6)) 
    cov_omega_rad_yr=cov_input .* ((pi/180)^2) .* (10^(-6))^2

    v = K_matrix * omega_rad_yr
    covariance = K_matrix * cov_omega_rad_yr * (K_matrix')

    return v[1], v[2], covariance # east, north [mm/yr]
end


"""
    estimate_euler_velocity(stations)

Return:
omega = [wx,wy,wz] #deg/Myr
"""
function estimate_euler_velocity(stations) #stations is a Matrix lon-lat-ve-vn-se-sn

    n = size(stations,1)

    G = zeros(2*n,3)
    d = zeros(2*n)
    W = zeros(2*n,2*n)

    R_EARTH=copy(Earth_Radius)
    R_EARTH=R_EARTH*(10^6); #mm

    for i in 1:n

        s=stations[i,:]
        λ = deg2rad(s[1])
        φ = deg2rad(s[2])

        slat = sin(φ)
        clat = cos(φ)
        slon = sin(λ)
        clon = cos(λ)

        G[2*i-1,:] = [
            -R_EARTH*clon*slat   -R_EARTH*slon*slat   R_EARTH*clat
        ]

        G[2*i,:] = [
            R_EARTH*slon       -R_EARTH*clon         0.0
        ]

        d[2*i-1] = s[3]
        d[2*i]   = s[4]

        W[2*i-1,2*i-1]= 1 / (s[5]^2)
        W[2*i,2*i]= 1/ (s[6]^2)

    end

    Wsqrt=sqrt.(W)

    G_W=Wsqrt*G
    d_W=Wsqrt*d
    
    # weighted lsq
    omega = G_W \ d_W

    N = (G_W')* G_W
    b=(G_W')*d_W
    I6 = Matrix(I, 3, 3)

    covariance = N \ I6

    chisq = (sum(abs2, d_W)) - ((b')*covariance*b)

    omega_deg_Myr = rad2deg.(omega) .* (10^6) 

    factor = (180/pi)*1e6
    covariance_deg_Myr = covariance .* (factor^2)

    d_pred=G * omega 
    residuals=d .- d_pred
    resE=residuals[1:2:end]
    resN=residuals[2:2:end]

    return omega_deg_Myr, covariance_deg_Myr, chisq, resE, resN

end

function estimate_euler_pole(stations)

    omega_deg_Myr, covariance_deg_Myr, chisq, resE, resN=estimate_euler_velocity(stations)
    omega_lon_deg, omega_lat_deg, magn_deg_Myr, covariance_result_deg=euler_pole_from_rotation(omega_deg_Myr; cov_omega=covariance_deg_Myr)

    return omega_lon_deg, omega_lat_deg, magn_deg_Myr, covariance_result_deg, chisq, resE, resN

end
