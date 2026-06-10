#!/usr/bin/env julia

# Julia translation of GAMIT/GLOBK velrot.f for standard velocity files.
#
# The code follows the Fortran program's data flow:
#   read sys1 -> frame update -> read sys2 -> frame update -> build links
#   -> estimate Helmert rate parameters -> update sys1 -> write combined field.
#
# Some GAMIT helper routines are not present in this workspace.  The missing
# geodetic routines are implemented here with WGS84 ellipsoid constants.  Frame
# rotations are treated as zero unless explicit values are added to
# frame_rotation_vector().

using LinearAlgebra
using Printf
using Dates

const VELROT_VER = "1.01-julia"
const A_WGS84 = 6378137.0
# GAMIT const_param.h value.  A local geod_to_xyz.f in this folder uses
# 298.257223563, but velrot.f calls the GAMIT library GEOD_to_XYZ signature,
# which is tied to const_param.h.
const F_WGS84 = 1.0 / 298.257222101
const E2_WGS84 = F_WGS84 * (2.0 - F_WGS84)
const RAD_TO_MAS = (180.0 / pi) * 3600.0 * 1000.0

mutable struct Site
    lon::Float64
    lat::Float64
    xyz::Vector{Float64}
    vel_xyz::Vector{Float64}
    cov_neu::Matrix{Float64}
    name::String
end

mutable struct Options
    sys1_file::String
    sys1_frame::String
    sys2_file::String
    sys2_frame::String
    out_file::String
    out_frame::String
    fund_file::String
    height_weight::Float64
    param_opt::String
    num_parn::Int
    eq_dist::Float64
    cp_dist::Float64
    av_dist::Float64
end

mutable struct FitResult
    trans_parm::Vector{Float64}
    trans_parm_out::Vector{Float64}
    trans_sigma::Vector{Float64}
    cov_parm::Matrix{Float64}
    sum_prefit::Float64
    sum_postfit::Float64
    sum_weight::Float64
    num_data::Int
    chi_fit::Float64
    rms_fit_mm::Float64
end

function usage()
    println("""
    Usage:
      julia velrot.jl <sys1> <frame1> <sys2> <frame2> <outname> <out_frame> <fundamental sites> <height wght> <param_opt>

    Example:
      julia velrot.jl Slave.vel eura Master.vel eura output_julia.txt eura link.file 1 TR
    """)
end

function parse_options(args)::Options
    length(args) >= 1 || (usage(); error("missing sys1 file"))
    sys1_file = args[1]
    sys1_frame = uppercase(get(args, 2, "NONE"))
    sys2_file = get(args, 3, "")
    sys2_frame = uppercase(get(args, 4, isempty(sys2_file) ? "NONE" : ""))
    out_file = get(args, 5, "6")
    out_frame = uppercase(get(args, 6, sys1_frame))
    fund_file = get(args, 7, "")
    height_weight = length(args) >= 8 ? max(parse(Float64, args[8]), 1e-6) : 1.0
    param_opt = uppercase(get(args, 9, "TR"))

    num_parn = 6
    if occursin("T", param_opt)
        num_parn = 3
    end
    if occursin("R", param_opt)
        num_parn = 6
    end
    if occursin("S", param_opt)
        num_parn = 7
    end
    if occursin("L", param_opt)
        num_parn = 2
    end
    if num_parn > 7
        @warn "Too many options in $param_opt; setting options to TR"
        param_opt = "TR"
        num_parn = 6
    end

    return Options(sys1_file, sys1_frame, sys2_file, sys2_frame,
                   out_file, out_frame, fund_file, height_weight,
                   param_opt, num_parn, 0.0, 0.0, 0.0)
end

function geod_to_xyz(lon_deg::Float64, lat_deg::Float64, h::Float64=0.0)
    lon = deg2rad(lon_deg)
    lat = deg2rad(lat_deg)
    sinφ = sin(lat)
    cosφ = cos(lat)
    n = A_WGS84 / sqrt(1.0 - E2_WGS84 * sinφ^2)
    return [(n + h) * cosφ * cos(lon),
            (n + h) * cosφ * sin(lon),
            (n * (1.0 - E2_WGS84) + h) * sinφ]
end

function xyz_to_geod(xyz::AbstractVector{<:Real})
    x, y, z = xyz
    equ_rad = hypot(x, y)
    lat_p = atan(z, equ_rad)
    h_p = 0.0
    lat_i = lat_p
    h_i = h_p
    for _ in 1:50
        rad_curve = A_WGS84 / sqrt(1.0 - E2_WGS84 * sin(lat_p)^2)
        rad_lat = equ_rad * (1.0 - E2_WGS84 * rad_curve / (rad_curve + h_p))
        lat_i = atan(z, rad_lat)
        h_i = abs(lat_i) < pi / 4 ?
              equ_rad / cos(lat_i) - rad_curve :
              z / sin(lat_i) - (1.0 - E2_WGS84) * rad_curve
        if abs(h_i - h_p) < 1e-4 &&
           abs(lat_i - lat_p) * rad_curve < 1e-4
            break
        end
        h_p = h_i
        lat_p = lat_i
    end
    lon = atan(y, x)
    lon < 0 && (lon += 2pi)
    return rad2deg(lon), rad2deg(lat_i), h_i
end

function xyz_to_neu_rotation(pos_xyz::AbstractVector{<:Real})
    x, y, z = pos_xyz
    colat = atan(hypot(x, y), z)
    lon = atan(y, x)
    radius = norm(pos_xyz)
    return [-cos(colat)*cos(lon)  -cos(colat)*sin(lon)   sin(colat);
            -sin(lon)              cos(lon)              0.0;
             x/radius              y/radius              z/radius]
end

function rotate_geod(v::AbstractVector{<:Real}, direction::Symbol, pos_xyz::AbstractVector{<:Real})
    lon, lat, _ = xyz_to_geod(pos_xyz)
    r = xyz_to_neu_rotation(pos_xyz)
    if direction == :XYZ_to_NEU
        return r * collect(v), (lon, lat), r
    elseif direction == :NEU_to_XYZ
        return transpose(r) * collect(v), (lon, lat), r
    else
        error("unknown rotation direction $direction")
    end
end

function read_sys(path::String)
    sites = Site[]
    open(path, "r") do io
        for line in eachline(io)
            isempty(strip(line)) && continue
            startswith(line, " ") || continue
            fields = split(line)
            length(fields) >= 13 || continue
            parsed = tryparse.(Float64, fields[1:12])
            any(isnothing, parsed) && continue
            values = Float64[x::Float64 for x in parsed]
            name = fields[13][1:min(end, 8)]

            lon = values[1]
            lat = values[2]
            xyz = geod_to_xyz(lon, lat, 0.0)
            vel_neu = [values[4], values[3], values[10]] ./ 1000.0
            vel_xyz, _, _ = rotate_geod(vel_neu, :NEU_to_XYZ, xyz)

            cov = zeros(Float64, 3, 3)
            cov[1, 1] = values[8]^2 * 1e-6
            cov[2, 2] = values[7]^2 * 1e-6
            cov[1, 2] = values[9] * values[7] * values[8] * 1e-6
            cov[2, 1] = cov[1, 2]
            cov[3, 3] = values[12]^2 * 1e-6
            push!(sites, Site(lon, lat, xyz, vel_xyz, cov, name))
        end
    end
    @printf(" There are %5d sites in sys file %s\n", length(sites), path)
    return sites
end

function fortran_data_block(text::String, marker::String)
    lines = split(text, '\n')
    start_idx = findfirst(line -> occursin(marker, line), lines)
    start_idx === nothing && return ""
    parts = String[]
    started = false
    for raw in lines[start_idx:end]
        line = raw
        if !started
            slash = findfirst('/', line)
            slash === nothing && continue
            line = line[(slash + 1):end]
            started = true
        end
        stripped = strip(line)
        if isempty(stripped) || first(stripped) in ('*', 'c', 'C')
            continue
        end
        slash = findfirst('/', line)
        if slash !== nothing
            push!(parts, line[1:slash - 1])
            break
        end
        push!(parts, line)
    end
    return join(parts, "\n")
end

function strip_fortran_comments(block::String)
    kept = String[]
    for raw in split(block, '\n')
        line = strip(raw)
        isempty(line) && continue
        first(line) in ('*', 'c', 'C') && continue
        bang = findfirst('!', line)
        if bang !== nothing
            line = strip(line[1:bang - 1])
        end
        push!(kept, line)
    end
    return join(kept, "\n")
end

function builtin_frames()
    path = joinpath(@__DIR__, "frame_to_fra.f")
    isfile(path) || return Dict{String,Vector{Float64}}()
    text = read(path, String)
    names_block = fortran_data_block(text, "data frame_names")
    data_block = fortran_data_block(text, "data frame_data")
    names = [uppercase(strip(m.captures[1])) for m in eachmatch(r"'([^']+)'", names_block)]
    numeric_text = replace(strip_fortran_comments(data_block), "D" => "e", "d" => "e")
    values = [parse(Float64, m.match) for m in eachmatch(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?", numeric_text)]
    frames = Dict{String,Vector{Float64}}()
    n = min(length(names), length(values) ÷ 3)
    for i in 1:n
        frames[names[i]] = values[(3i - 2):(3i)]
    end
    return frames
end

function external_frames()
    frames = Dict{String,Vector{Float64}}()
    for path in ("frames.dat", joinpath(homedir(), "gg", "tables", "frames.dat"))
        isfile(path) || continue
        open(path, "r") do io
            for line in eachline(io)
                isempty(strip(line)) && continue
                startswith(line, " ") || continue
                fields = split(line)
                length(fields) >= 4 || continue
                vals = try
                    parse.(Float64, fields[2:4])
                catch
                    continue
                end
                unit = length(fields) >= 5 ? uppercase(fields[5]) : ""
                if startswith(unit, "RAD")
                    vals .*= 180.0 / pi
                end
                name = uppercase(fields[1][1:min(end, 8)])
                haskey(frames, name) || (frames[name] = vals)
            end
        end
    end
    return frames
end

function frame_rotation_vector(sys_frame::String, out_frame::String)
    sys = uppercase(strip(sys_frame))
    out = uppercase(strip(out_frame))
    sys == out && return zeros(3)
    startswith(sys, "NONE") && return zeros(3)

    frames = builtin_frames()
    merge!(frames, external_frames())
    sys_key = sys[1:min(end, 8)]
    out_key = out[1:min(end, 8)]
    if !haskey(frames, sys_key) || !haskey(frames, out_key)
        @warn "Could not find frame definition; treating $sys_frame -> $out_frame as zero rotation"
        return zeros(3)
    end

    scs = occursin(":O", sys) && !startswith(sys, "ITRF") ? 1.0 / 0.9562 : 1.0
    sco = occursin(":O", out) && !startswith(out, "ITRF") ? 1.0 / 0.9562 : 1.0
    return frames[sys_key] .* scs .* (pi / 180e6) .-
           frames[out_key] .* sco .* (pi / 180e6)
end

function frame_update!(sites::Vector{Site}, sys_frame::String, out_frame::String)
    rot_vec = frame_rotation_vector(sys_frame, out_frame)
    @printf(" Rotating from %-10s to %-10s using rotation vector %12.6f%12.6f%12.6f degs/Myrs\n",
            sys_frame, out_frame, rot_vec[1] * 180e6 / pi,
            rot_vec[2] * 180e6 / pi, rot_vec[3] * 180e6 / pi)
    sum(abs, rot_vec) < 1e-10 && return
    for site in sites
        site.vel_xyz .+= cross(rot_vec, site.xyz)
    end
end

function find_site(name::String, sites::Vector{Site})
    uname = uppercase(name)
    return findfirst(s -> uppercase(s.name) == uname, sites)
end

function distance(a::Site, b::Site)
    return norm(a.xyz - b.xyz)
end

function read_fund_sites!(opt::Options, sys1::Vector{Site}, sys2::Vector{Site})
    links = Tuple{Int,Int}[]
    if !isempty(strip(opt.fund_file))
        open(opt.fund_file, "r") do io
            for line in eachline(io)
                isempty(strip(line)) && continue
                startswith(line, " ") || continue
                fields = split(line)
                isempty(fields) && continue
                key = uppercase(fields[1])
                if key == "EQ_DIST"
                    opt.eq_dist = parse(Float64, fields[2])
                    opt.cp_dist = opt.eq_dist
                    empty!(links)
                    for (j, s1) in pairs(sys1), (k, s2) in pairs(sys2)
                        distance(s1, s2) < opt.eq_dist && push!(links, (j, k))
                    end
                elseif key == "NAMES"
                    empty!(links)
                    for (j, s1) in pairs(sys1)
                        k = find_site(s1.name, sys2)
                        k !== nothing && push!(links, (j, k))
                    end
                elseif key == "CP_DIST"
                    opt.cp_dist = parse(Float64, fields[2])
                elseif key == "AV_DIST"
                    opt.av_dist = parse(Float64, fields[2])
                else
                    name1 = fields[1]
                    name2 = length(fields) >= 2 ? fields[2] : name1
                    link1 = !startswith(name1, "-")
                    link2 = !startswith(name2, "-")
                    name1 = replace(name1, r"^[+-]" => "")
                    name2 = replace(name2, r"^[+-]" => "")
                    j = find_site(name1, sys1)
                    k = find_site(name2, sys2)
                    if j !== nothing && k !== nothing && link1 && link2
                        push!(links, (j, k))
                    end
                    if j !== nothing && !link1
                        filter!(lk -> lk[1] != j, links)
                    end
                    if k !== nothing && !link2
                        filter!(lk -> lk[2] != k, links)
                    end
                end
            end
        end
    else
        for (j, s1) in pairs(sys1)
            k = find_site(s1.name, sys2)
            k !== nothing && push!(links, (j, k))
        end
    end

    deduped = Tuple{Int,Int}[]
    seen = Set{Tuple{Int,Int}}()
    for lk in links
        if !(lk in seen)
            push!(deduped, lk)
            push!(seen, lk)
        end
    end
    @printf(" There are %5d matching sites in fundamental file %s\n", length(deduped), opt.fund_file)
    return deduped
end

function get_parts(coord::AbstractVector{<:Real}, rot_matrix::AbstractMatrix{<:Real}, num_parn::Int)
    xyz_part = zeros(Float64, 3, 7)
    xyz_part[1, 1] = 1.0
    xyz_part[2, 2] = 1.0
    xyz_part[3, 3] = 1.0
    xyz_part[2, 4] = -coord[3]
    xyz_part[3, 4] =  coord[2]
    xyz_part[1, 5] =  coord[3]
    xyz_part[3, 5] = -coord[1]
    xyz_part[1, 6] = -coord[2]
    xyz_part[2, 6] =  coord[1]
    if num_parn == 7
        xyz_part[:, 7] .= coord
    end
    neu_part = zeros(Float64, 3, 7)
    neu_part[:, 1:num_parn] .= rot_matrix * xyz_part[:, 1:num_parn]
    return neu_part, xyz_part
end

function estimate_transform(opt::Options, sys1::Vector{Site}, sys2::Vector{Site}, links)
    npar = opt.num_parn
    norm_eq = zeros(Float64, 7, 7)
    bvec = zeros(Float64, 7)
    if !occursin("T", opt.param_opt) && !occursin("L", opt.param_opt)
        for i in 1:3
            norm_eq[i, i] = 1e14
        end
    end
    if !occursin("S", opt.param_opt)
        norm_eq[7, 7] = 1e14
    end

    sum_prefit = 0.0
    sum_weight = 0.0
    num_data = 0
    for (ns1, ns2) in links
        s1 = sys1[ns1]
        s2 = sys2[ns2]
        dx = s2.vel_xyz - s1.vel_xyz
        dn, _, rot_matrix = rotate_geod(dx, :XYZ_to_NEU, s1.xyz)
        neu_part, _ = get_parts(s1.xyz, rot_matrix, npar)
        weights = [1.0 / (s1.cov_neu[1, 1] + s2.cov_neu[1, 1]),
                   1.0 / (s1.cov_neu[2, 2] + s2.cov_neu[2, 2]),
                   opt.height_weight / (s1.cov_neu[3, 3] + s2.cov_neu[3, 3])]
        num_use = abs(weights[3] / weights[2]) < 1e-5 ? 2 : 3
        for i in 1:num_use
            for j in 1:7
                bvec[j] += neu_part[i, j] * dn[i] * weights[i]
                for k in 1:7
                    norm_eq[j, k] += neu_part[i, j] * weights[i] * neu_part[i, k]
                end
            end
            sum_prefit += dn[i]^2 * weights[i]
            sum_weight += weights[i]
            num_data += 1
        end
    end

    cov_full = zeros(Float64, 7, 7)
    trans = zeros(Float64, 7)
    if num_data - npar > 0
        cov_active = inv(norm_eq[1:npar, 1:npar])
        trans[1:npar] .= cov_active * bvec[1:npar]
        cov_full[1:npar, 1:npar] .= cov_active
        dprefit = dot(bvec[1:npar], cov_active * bvec[1:npar])
        sum_postfit = sum_prefit - dprefit
        sum_postfit = max(sum_postfit, 0.0)
        chi_fit = sqrt(sum_postfit / (num_data - npar))
        rms_fit = sqrt(num_data / sum_weight) * chi_fit
        if chi_fit < 0.10
            @printf("# Chi of fit less than 0.10 (%6.4f). Resetting to 1.0\n", chi_fit)
            chi_fit = 1.0
        end
    else
        sum_postfit = 0.0
        chi_fit = 1.0
        rms_fit = 0.0
    end

    trans_out = zeros(Float64, 7)
    trans_out[1:3] .= trans[1:3] .* 1000.0
    trans_out[4:6] .= trans[4:6] .* RAD_TO_MAS
    trans_out[7] = trans[7] * 1e9

    sig = zeros(Float64, 7)
    sig[1:3] .= sqrt.(diag(cov_full)[1:3]) .* 1000.0 .* chi_fit
    sig[4:6] .= sqrt.(diag(cov_full)[4:6]) .* RAD_TO_MAS .* chi_fit
    sig[7] = sqrt(cov_full[7, 7]) * 1e9 * chi_fit

    return FitResult(trans, trans_out, sig, cov_full, sum_prefit,
                     sum_postfit, sum_weight, num_data, chi_fit,
                     rms_fit * 1000.0)
end

function transform_cov(neu_part::Matrix{Float64}, cov_parm::Matrix{Float64})
    return neu_part * cov_parm * transpose(neu_part)
end

function residual_for_link(fit::FitResult, opt::Options, s1::Site, s2::Site)
    dx = s2.vel_xyz - s1.vel_xyz
    _, _, rot_matrix = rotate_geod(dx, :XYZ_to_NEU, s1.xyz)
    neu_part, xyz_part = get_parts(s1.xyz, rot_matrix, opt.num_parn)
    for j in 1:3
        for k in 1:opt.num_parn
            dx[j] -= xyz_part[j, k] * fit.trans_parm[k]
        end
    end
    dn, _, _ = rotate_geod(dx, :XYZ_to_NEU, s1.xyz)
    cov_neu = transform_cov(neu_part, fit.cov_parm)
    return dn, cov_neu
end

function write_header(io, opt::Options, fit::FitResult, nfund::Int)
    now = Dates.now()
    @printf(io, "* VELROT Run on %4d/%2d/%2d %2d:%2d Version %s\n",
            year(now), month(now), day(now), hour(now), minute(now), VELROT_VER)
    @printf(io, "* SYSTEM 1 File    : %-10s %s\n", opt.sys1_frame, opt.sys1_file)
    @printf(io, "* SYSTEM 2 File    : %-10s %s\n", opt.sys2_frame, opt.sys2_file)
    @printf(io, "* FUNDAMENTAL File :            %s\n", opt.fund_file)
    @printf(io, "* OUTPUT FRAME     : %-10s PARAM_OPT %-8s\n", opt.out_frame, opt.param_opt)
    @printf(io, "* EQ_DIST          : %18.1f m, CP_DIST : %12.1f m*\n", opt.eq_dist, opt.cp_dist)
    @printf(io, "* AV_DIST          : %18.1f m* HEIGHT WEIGHT    : %10.6f\n", opt.av_dist, opt.height_weight)
    if fit.num_data > 0
        @printf(io, "* \n* RMS fit for %8d components from %4d stations was %10.2f mm/yr, NRMS %8.2f\n",
                fit.num_data, nfund, fit.rms_fit_mm, fit.chi_fit)
        println(io, "* Estimates of Transformation parameters are: ")
        labels = ["X-Offset", "Y-Offset", "Z-Offset", "X-Rot", "Y-Rot", "Z-Rot", "Scale"]
        units = ["(mm/yr)", "(mm/yr)", "(mm/yr)", "(mas/yr)", "(mas/yr)", "(mas/yr)", "(ppb/yr)"]
        for i in 1:opt.num_parn
            @printf(io, "*%5d %-8s   %10.4f %10.4f %s\n", i, labels[i], fit.trans_parm_out[i], fit.trans_sigma[i], units[i])
        end
    end
end

function output_sum(io, opt::Options, fit::FitResult, sys1::Vector{Site}, sys2::Vector{Site}, links)
    write_header(io, opt, fit, length(links))
    println(io, "*  Differences at the fundamental sites")
    println(io, "*   #  Name 1  Name Ref    dN (mm)    dE (mm)    dU (mm)    sN (mm)    sE (mm)    sU (mm)   sTN (mm)   sTE (mm)   sTU (mm)")
    summ = zeros(4)
    sumv = zeros(4)
    sumw = zeros(4)
    for (i, (ns1, ns2)) in enumerate(links)
        s1, s2 = sys1[ns1], sys2[ns2]
        dn, cov_neu = residual_for_link(fit, opt, s1, s2)
        sig_obs = sqrt.(diag(s1.cov_neu + s2.cov_neu)) .* 1000.0
        sig_tr = sqrt.(max.(diag(cov_neu), 0.0)) .* 1000.0
        @printf(io, "A%5d %-8s %-8s %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
                i, s1.name, s2.name, dn[1] * 1000.0, dn[2] * 1000.0, dn[3] * 1000.0,
                sig_obs[1], sig_obs[2], sig_obs[3], sig_tr[1], sig_tr[2], sig_tr[3])
        for j in 1:3
            v = s1.cov_neu[j, j] + s2.cov_neu[j, j]
            summ[j] += dn[j] / v
            sumv[j] += dn[j]^2 / v
            sumw[j] += 1.0 / v
        end
    end
    summ[4] = (summ[1] + summ[2]) / 2.0
    sumv[4] = (sumv[1] + sumv[2]) / 2.0
    sumw[4] = (sumw[1] + sumw[2]) / 2.0
    for (j, lab) in enumerate(["North", "East", "Up", "Horz"])
        wmean = summ[j] / sumw[j] * 1000.0
        nrms = sqrt(sumv[j] / length(links))
        wrms = sqrt((1.0 / sumw[j]) * length(links)) * nrms * 1000.0
        @printf("S Component %-5s # %5d WMean %6.2f WRMS %6.2f mm/yr, NRMS %7.3f\n",
                lab, length(links), wmean, wrms, nrms)
    end
end

function update_tran!(fit::FitResult, opt::Options, sys1::Vector{Site})
    for s1 in sys1
        _, _, rot_matrix = rotate_geod(s1.vel_xyz, :XYZ_to_NEU, s1.xyz)
        neu_part, xyz_part = get_parts(s1.xyz, rot_matrix, opt.num_parn)
        dx = xyz_part[:, 1:opt.num_parn] * fit.trans_parm[1:opt.num_parn]
        s1.vel_xyz .+= dx
        cov_neu = transform_cov(neu_part, fit.cov_parm)
        dsig_e = sqrt(max(s1.cov_neu[2, 2] + cov_neu[2, 2], 0.0)) * 1000.0
        dsig_n = sqrt(max(s1.cov_neu[1, 1] + cov_neu[1, 1], 0.0)) * 1000.0
        dsig_u = sqrt(max(s1.cov_neu[3, 3] + cov_neu[3, 3], 0.0)) * 1000.0
        rho = (s1.cov_neu[1, 2] + cov_neu[1, 2]) / (dsig_n * dsig_e / 1e6)
        s1.cov_neu[1, 1] = dsig_n^2 / 1e6
        s1.cov_neu[2, 2] = dsig_e^2 / 1e6
        s1.cov_neu[3, 3] = dsig_u^2 / 1e6
        s1.cov_neu[1, 2] = rho * sqrt(s1.cov_neu[1, 1] * s1.cov_neu[2, 2])
        s1.cov_neu[2, 1] = s1.cov_neu[1, 2]
    end
end

function average_cluster!(cluster::Vector{Site})
    isempty(cluster) && return
    sum_av = zeros(3)
    sum_var = zeros(3)
    sum_rho = 0.0
    for site in cluster
        dn, _, _ = rotate_geod(site.vel_xyz, :XYZ_to_NEU, site.xyz)
        for j in 1:3
            sum_av[j] += dn[j] / site.cov_neu[j, j]
            sum_var[j] += 1.0 / site.cov_neu[j, j]
        end
        sum_rho += site.cov_neu[1, 2] / sqrt(site.cov_neu[1, 1] * site.cov_neu[2, 2])
    end
    mean_neu = sum_av ./ sum_var
    mean_var = 1.0 ./ sum_var
    mean_rho = sum_rho / length(cluster)
    for site in cluster
        site.vel_xyz, _, _ = rotate_geod(mean_neu, :NEU_to_XYZ, site.xyz)
        site.cov_neu .= 0.0
        for j in 1:3
            site.cov_neu[j, j] = mean_var[j]
        end
        site.cov_neu[1, 2] = mean_rho * sqrt(mean_var[1] * mean_var[2])
        site.cov_neu[2, 1] = site.cov_neu[1, 2]
    end
end

function av_frame!(opt::Options, sys1::Vector{Site}, sys2::Vector{Site})
    opt.av_dist <= 0 && return
    used1 = falses(length(sys1))
    used2 = falses(length(sys2))

    for i in eachindex(sys1)
        cluster = Site[sys1[i]]
        idx1 = Int[i]
        idx2 = Int[]
        for k in (i + 1):length(sys1)
            if !used1[k] && distance(sys1[i], sys1[k]) < opt.av_dist
                push!(cluster, sys1[k])
                push!(idx1, k)
            end
        end
        for k in eachindex(sys2)
            if !used2[k] && distance(sys1[i], sys2[k]) < opt.av_dist
                push!(cluster, sys2[k])
                push!(idx2, k)
            end
        end
        if length(cluster) > 1
            average_cluster!(cluster)
            used1[idx1] .= true
            used2[idx2] .= true
        end
    end

    for i in eachindex(sys2)
        used2[i] && continue
        cluster = Site[sys2[i]]
        idx2 = Int[i]
        for k in (i + 1):length(sys2)
            if !used2[k] && distance(sys2[i], sys2[k]) < opt.av_dist
                push!(cluster, sys2[k])
                push!(idx2, k)
            end
        end
        if length(cluster) > 1
            average_cluster!(cluster)
            used2[idx2] .= true
        end
    end
end

function check_cp(site::Site, others::Vector{Site}, dist::Float64, symbol::Char)
    dist <= 0 && return ' '
    for other in others
        distance(site, other) < dist && return symbol
    end
    return ' '
end

function write_velocity_line(io, prefix::Char, site::Site, symbol::Char)
    dn, (lon, lat), _ = rotate_geod(site.vel_xyz, :XYZ_to_NEU, site.xyz)
    lon = lon < 0 ? lon + 360.0 : lon
    dsig_n = sqrt(site.cov_neu[1, 1]) * 1000.0
    dsig_e = sqrt(site.cov_neu[2, 2]) * 1000.0
    dsig_u = sqrt(site.cov_neu[3, 3]) * 1000.0
    rho = site.cov_neu[1, 2] / (dsig_n * dsig_e / 1e6)
    @printf(io, "%c%10.5f %10.5f %8.2f %7.2f %7.2f %7.2f %7.2f %7.2f %6.3f  %8.2f %7.2f %7.2f %-8s%c\n",
            prefix, lon, lat, dn[2] * 1000.0, dn[1] * 1000.0,
            dn[2] * 1000.0, dn[1] * 1000.0, dsig_e, dsig_n, rho,
            dn[3] * 1000.0, dn[3] * 1000.0, dsig_u, site.name, symbol)
end

function output_frame(io, opt::Options, sys1::Vector{Site}, sys2::Vector{Site})
    println(io, "\n* SYSTEM 1 Velocities transformed to SYSTEM 2 ")
    println(io, "*   Long.       Lat.         E & N Rate      E & N Adj.      E & N +-   RHO        H Rate   H adj.    +-  SITE")
    println(io, "*  (deg)      (deg)           (mm/yr)       (mm/yr)       (mm/yr)                 (mm/yr)")
    for s1 in sys1
        write_velocity_line(io, ' ', s1, check_cp(s1, sys2, opt.cp_dist, '*'))
    end

    println(io, "\n* SYSTEM 2 Velocities except those in SYSTEM 1 ")
    println(io, "*   Long.       Lat.         E & N Rate      E & N Adj.      E & N +-   RHO        H Rate   H adj.    +-  SITE")
    println(io, "*  (deg)      (deg)           (mm/yr)       (mm/yr)       (mm/yr)                 (mm/yr)")
    sys1_names = Set(s.name for s in sys1)
    for s2 in sys2
        prefix = s2.name in sys1_names ? '-' : ' '
        write_velocity_line(io, prefix, s2, check_cp(s2, sys1, opt.cp_dist, '+'))
    end
end

function run_velrot(args=ARGS)
    opt = parse_options(args)
    println("\n VELROT: Velocity field comparison and combination Version $VELROT_VER\n")
    sys1 = read_sys(opt.sys1_file)
    frame_update!(sys1, opt.sys1_frame, opt.out_frame)
    sys2 = read_sys(opt.sys2_file)
    frame_update!(sys2, opt.sys2_frame, opt.out_frame)
    links = read_fund_sites!(opt, sys1, sys2)
    fit = estimate_transform(opt, sys1, sys2, links)

    if opt.out_file == "6"
        output_sum(stdout, opt, fit, sys1, sys2, links)
        update_tran!(fit, opt, sys1)
        av_frame!(opt, sys1, sys2)
        output_frame(stdout, opt, sys1, sys2)
    else
        open(opt.out_file, "w") do io
            output_sum(io, opt, fit, sys1, sys2, links)
            update_tran!(fit, opt, sys1)
            av_frame!(opt, sys1, sys2)
            output_frame(io, opt, sys1, sys2)
        end
    end
    return fit
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_velrot()
end
