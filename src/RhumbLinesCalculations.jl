# Define the Earth's radius
# https://www.jpz.se/Html_filer/wgs_84.html
global const Earth_Radius = 6371.0087714 #km
global const Earth_SemimajorAxis = 6378.137 #km
global const Earth_SemiminorAxis = 6356.75231424518 #km
global const Earth_Eccentricity=0.081819190842622

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


function rhumb_distance_WGS84(LonLat1::Vector{<:Real}, LonLat2::Vector{<:Real}) # Not super-sure about it; check needed

    lon_a = deg2rad(LonLat1[1])
    lat_a = deg2rad(LonLat1[2])

    lon_b = deg2rad(LonLat2[1])
    lat_b = deg2rad(LonLat2[2])
    
    Δλ = lon_b - lon_a
    ΔL = lat_b - lat_a
    e_=Earth_Eccentricity

    Σ_b = log(sec(lat_b)+tan(lat_b))-(e_/2)*(log(1+e_*sin(lat_b))-log(1-e_*sin(lat_b)))
    Σ_a = log(sec(lat_a)+tan(lat_a))-(e_/2)*(log(1+e_*sin(lat_a))-log(1-e_*sin(lat_a)))

    ΔΣ=Σ_b-Σ_a

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