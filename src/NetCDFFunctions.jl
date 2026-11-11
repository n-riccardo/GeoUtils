using NCDatasets

"""
    write_to_netcdf(filename, lon, lat, values; lon_units="degrees_east", lat_units="degrees_north", 
                    value_name="z", value_units=" ", value_desc=" ")

Write longitude, latitude, and corresponding values to a NetCDF file (lon-lat order).
Arguments:
- filename: String, name of the output NetCDF file
- lon: Vector of unique longitude values
- lat: Vector of unique latitude values
- values: Matrix of values corresponding to (lon, lat) grid
Keyword Arguments:
- lon_units: Units for the longitude variable (default: "degrees_east")
- lat_units: Units for the latitude variable (default: "degrees_north")
- value_name: Name of the values variable (default: "z")
- value_units: Units of the values variable (default: " ")
- value_desc: Description of the values variable (default: " ")
"""
function write_to_netcdf(filename::String, lon::Vector{<:Real}, lat::Vector{<:Real}, values::Matrix{<:Real}; 
                         lon_units::String="degrees_east", lat_units::String="degrees_north", 
                         value_name::String="z", value_units::String=" ", value_desc::String=" ")

    ds = NCDataset(filename, "c")

        # Define dimensions
        defDim(ds, "lon", length(lon))
        defDim(ds, "lat", length(lat))
        
        # Define variables
        lon_var = defVar(ds, "lon", Float64, ("lon",))
        lat_var = defVar(ds, "lat", Float64, ("lat",))
        value_var = defVar(ds, value_name, Float64, ("lon", "lat"))  # (lon, lat) order

        # Add attributes to variables
        lon_var.attrib["units"] = lon_units
        lon_var.attrib["long_name"] = "longitude"
        lat_var.attrib["units"] = lat_units
        lat_var.attrib["long_name"] = "latitude"
        value_var.attrib["units"] = value_units
        value_var.attrib["long_name"] = value_desc
        
        # Write data to variables
        lon_var[:] = lon
        lat_var[:] = lat
        value_var[:, :] = values  # Write values in (lon, lat) order

        # Close the NetCDF file
        close(ds)

end
