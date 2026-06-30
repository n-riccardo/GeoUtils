############################################################
### VISR R/W
############################################################

function write_visr_input(filename,
    Site, Long, Lat,
    E_Rate, sigma_E,
    N_Rate, sigma_N,
    rho)

n = length(Site)

open(filename, "w") do io
for i in 1:n

# limita nome a 8 caratteri (Fortran a8)
site = rpad(Site[i][1:min(end,8)], 8)

lon = Long[i]
lat = Lat[i]

ux = E_Rate[i]
sx = sigma_E[i]

uy = N_Rate[i]
sy = sigma_N[i]

cxy = rho[i]

@printf(io,
"%-8s%10.4f%10.4f%7.2f%5.2f%7.2f%5.2f%7.3f\n",
site, lon, lat, ux, sx, uy, sy, cxy)
end
end
end

"""
    read_strain_output_visr(filepath::String)
 
Legge un file di output del codice Fortran per la deformazione della crosta.
 
Il file contiene una riga di intestazione seguita da righe dati con il formato:
  longitude latitude vx dvx vy dvy cxy w dw exx dexx exy dexy eyy deyy emax demax emin demin shr dshr azi dazi dilat ddilat dis weight chisq nsite [index = N]
 
Le righe con `index = 1` contengono dati completi (29 valori numerici).
Le righe con `index = 2` (o altri valori > 1) contengono solo lon/lat;
tutti gli altri campi vengono impostati a `NaN`.
 
# Valori restituiti
Restituisce una `NamedTuple` con i seguenti campi:
 
## Vettori di coordinate univoci (ordinati in senso crescente)
- `lons`    :: Vector{Float64}  — longitudini univoche (N_lon,)
- `lats`    :: Vector{Float64}  — latitudini univoche  (N_lat,)
 
## Matrici 2D per ogni campo fisico
Dimensioni: (N_lon, N_lat), con lon crescente lungo le righe e lat crescente lungo le colonne.
- `vx`, `dvx`         :: Matrix{Float64}   — velocità est e sua incertezza [mm/yr]
- `vy`, `dvy`         :: Matrix{Float64}   — velocità nord e sua incertezza [mm/yr]
- `cxy`               :: Matrix{Float64}   — correlazione vx-vy
- `w`, `dw`           :: Matrix{Float64}   — rotazione e sua incertezza [nrad/yr]
- `exx`, `dexx`       :: Matrix{Float64}   — strain est-est e sua incertezza [nstrain/yr]
- `exy`, `dexy`       :: Matrix{Float64}   — strain est-nord e sua incertezza [nstrain/yr]
- `eyy`, `deyy`       :: Matrix{Float64}   — strain nord-nord e sua incertezza [nstrain/yr]
- `emax`, `demax`     :: Matrix{Float64}   — strain principale massimo e sua incertezza [nstrain/yr]
- `emin`, `demin`     :: Matrix{Float64}   — strain principale minimo e sua incertezza [nstrain/yr]
- `shr`, `dshr`       :: Matrix{Float64}   — shear massimo e sua incertezza [nstrain/yr]
- `azi`, `dazi`       :: Matrix{Float64}   — azimut dell'asse di strain max e sua incertezza [deg]
- `dilat`, `ddilat`   :: Matrix{Float64}   — dilatazione e sua incertezza [nstrain/yr]
- `dis`               :: Matrix{Float64}   — 
- `weight`            :: Matrix{Float64}   — peso della stima
- `chisq`             :: Matrix{Float64}   — chi-quadro della soluzione
- `nsite`             :: Matrix{Float64}   — numero di siti usati
 
## Maschere booleane per valori di `index`
Matrici (N_lon, N_lat) che sono `true` dove `index == N`:
- `mask_index2`  :: Matrix{Bool}   — punti con index = 2
- (aggiuntivi per ogni valore di index > 1 trovato nel file)
  Restituiti come dizionario `index_masks :: Dict{Int, Matrix{Bool}}`
"""
function read_strain_output_visr(filepath::String)
 
    # ------------------------------------------------------------------ #
    # 1.  Raccoglie tutte le righe nel vettore raw        #
    # ------------------------------------------------------------------ #
    #longitude latitude      vx+dvx       vy+dvy    cxy         w+dw         exx+dexx         exy+dexy         eyy+deyy         emax+demax       emin+demin        shr+dshr         azi+dazi        dilat+ddilat     dis.    weight    chisq   nsite
    # Nomi dei 27 campi dopo lon/lat (nell'ordine in cui appaiono nel file)
    field_names = [
        :vx, :dvx, :vy, :dvy, :cxy,
        :w, :dw,
        :exx, :dexx, :exy, :dexy, :eyy, :deyy,
        :emax, :demax, :emin, :demin,
        :shr, :dshr, :azi, :dazi,
        :dilat, :ddilat,
        :dis, :weight, :chisq, :nsite
    ]
    n_fields = length(field_names)   # 27
 
    # Struttura temporanea per accumulare i dati riga per riga
    lons_raw    = Float64[]
    lats_raw    = Float64[]
    data_raw    = Vector{Float64}[]  # ogni elemento è un vettore di n_fields valori
    index_raw   = Int[]
 
    open(filepath, "r") do io
        readline(io)   # salta la riga di intestazione
        readline(io)   # salta la seconda riga

        for line in eachline(io)
            isempty(strip(line)) && continue
 
            if occursin("index", line)
                # Formato: "lon  lat  index = N"
                parts = split(line)
                lon   = parse(Float64, parts[1])
                lat   = parse(Float64, parts[2])
                idx   = parse(Int,     parts[end])
                push!(lons_raw,  lon)
                push!(lats_raw,  lat)
                push!(index_raw, idx)
                push!(data_raw,  fill(NaN, n_fields))
            else
                # Formato: "lon  lat  v1  v2  ... v27"
                parts = split(line)
                lon   = parse(Float64, parts[1])
                lat   = parse(Float64, parts[2])
                vals  = [parse(Float64, parts[i+2]) for i in 1:n_fields]
                push!(lons_raw,  lon)
                push!(lats_raw,  lat)
                push!(index_raw, -1) # -1 valore dove non ci sono righe "index="
                push!(data_raw,  vals)
            end
        end
    end
 
    # ------------------------------------------------------------------ #
    # 2.  Determina la griglia                                             #
    # ------------------------------------------------------------------ #
    lons_uniq = sort(unique(lons_raw))   # N_lon valori crescenti
    lats_uniq = sort(unique(lats_raw))   # N_lat valori crescenti

    longrid=[lons_uniq[i] for i in axes(lons_uniq,1), _ in axes(lats_uniq,1)]
    latgrid=[lats_uniq[j] for _ in axes(lons_uniq,1), j in axes(lats_uniq,1)]

    N_lon = length(lons_uniq)
    N_lat = length(lats_uniq)
 
    # Indici rapidi per lon/lat → posizione nella griglia
    lon_idx = Dict(v => i for (i, v) in enumerate(lons_uniq))
    lat_idx = Dict(v => i for (i, v) in enumerate(lats_uniq))
 
    # ------------------------------------------------------------------ #
    # 3.  Riempie le matrici 2D  (N_lon × N_lat)                          #
    # ------------------------------------------------------------------ #
    # Inizializza tutte le matrici a NaN
    matrices = [fill(NaN, N_lon, N_lat) for _ in 1:n_fields]
 
    # Trova tutti i valori di index presenti (eccetto -1)
    all_indices = sort(unique(index_raw))
    extra_indices = filter(x -> x != -1, all_indices)
 
    # Matrici booleane per ogni index != -1
    index_masks = Dict{Int, Matrix{Bool}}(
        k => fill(false, N_lon, N_lat) for k in extra_indices
    )
 
    for k in 1:length(lons_raw)
        i = lon_idx[lons_raw[k]]
        j = lat_idx[lats_raw[k]]
        idx = index_raw[k]
 
        # Riempie i campi fisici (NaN già presenti per index != -1)
        if idx == -1
            for f in 1:n_fields
                matrices[f][i, j] = data_raw[k][f]
            end
        end
 
        # Aggiorna la maschera booleana per i valori di index != -1
        if idx != -1
            index_masks[idx][i, j] = true
        end
    end
 
    # ------------------------------------------------------------------ #
    # 4.  Costruisce la NamedTuple di output                             #
    # ------------------------------------------------------------------ #
    # Associa ogni matrice al suo nome simbolico
    field_dict = Dict(field_names[f] => matrices[f] for f in 1:n_fields)
 
    return (
        # Coordinate
        lons   = lons_uniq,
        lats   = lats_uniq,
        longrid= longrid,
        latgrid=    latgrid,
 
        # Velocità
        vx     = field_dict[:vx],
        dvx    = field_dict[:dvx],
        vy     = field_dict[:vy],
        dvy    = field_dict[:dvy],
        cxy    = field_dict[:cxy],
 
        # Rotazione
        w      = field_dict[:w],
        dw     = field_dict[:dw],
 
        # Tensore di deformazione
        exx    = field_dict[:exx],
        dexx   = field_dict[:dexx],
        exy    = field_dict[:exy],
        dexy   = field_dict[:dexy],
        eyy    = field_dict[:eyy],
        deyy   = field_dict[:deyy],
 
        # Strain principali
        emax   = field_dict[:emax],
        demax  = field_dict[:demax],
        emin   = field_dict[:emin],
        demin  = field_dict[:demin],
 
        # Shear e azimut
        shr    = field_dict[:shr],
        dshr   = field_dict[:dshr],
        azi    = field_dict[:azi],
        dazi   = field_dict[:dazi],
 
        # Dilatazione
        dilat  = field_dict[:dilat],
        ddilat = field_dict[:ddilat],
 
        # Statistiche/metadati
        dis    = field_dict[:dis],
        weight = field_dict[:weight],
        chisq  = field_dict[:chisq],
        nsite  = field_dict[:nsite],
 
        # Maschere booleane per index != 1
        index_masks = index_masks,
    )
end


############################################################
### Wavelets R/W
############################################################


function write_wavformat(filename, Lon, Lat, Ve, Vn, Se, Sn, Name; 
    Vu=zeros(length(Lon)), Su=ones(length(Lon)), 
    Ren=zeros(length(Lon)), Reu=zeros(length(Lon)), Run=zeros(length(Lon)), 
    t1= (zeros(length(Lon)) .* NaN), t2=(zeros(length(Lon)) .* NaN))

VeloTseri = DataFrame(
    lon = Lon,
    lat = Lat,
    Ve  = Ve,
    Vn  = Vn,
    Vu  = Vu,
    Se  = Se,
    Sn  = Sn,
    Su  = Su,
    Ren = Ren,
    Reu = Reu,
    Run = Run,
    t1  = t1,
    t2  = t2,
    name = Name
)

CSV.write(filename, VeloTseri; delim=' ')

end


function wavelets_read(strain_result_path,velo_path,qmax_result_path,qmaxMask)

veloEastResults=CSV.read(velo_path*"east_plot.dat", delim=' ', header=false, ignorerepeated=true, DataFrame)
veloSouthResults=CSV.read(velo_path*"south_plot.dat", delim=' ', header=false, ignorerepeated=true, DataFrame)
veloResiduals=CSV.read(velo_path*"residual.dat", delim=' ', header=false, ignorerepeated=true, DataFrame)

LonV=veloEastResults[:,1];
LatV=veloEastResults[:,2];
VeloEast=veloEastResults[:,3];
VeloNorth=(veloSouthResults[:,3]) .* (-1);

ResLon=veloResiduals[:,1]
ResLat=veloResiduals[:,2]
ResN=-veloResiduals[:,7]
ResE=veloResiduals[:,8]

Lon=unique(LonV);
Lat=unique(LatV);
LonGrid=reshape(LonV, length(Lat), length(Lon));
LatGrid=reshape(LatV, length(Lat), length(Lon));
VeloEastGrid=reshape(VeloEast, length(Lat), length(Lon));
VeloNorthGrid=reshape(VeloNorth, length(Lat), length(Lon));

#------------------------------------------------

qmaxResults=CSV.read(qmax_result_path, delim=' ', header=false, ignorerepeated=true, DataFrame)

qmaxV=qmaxResults[:,3];
qmaxGrid=reshape(qmaxV, length(Lat), length(Lon));
qmaxGrid=Matrix(qmaxGrid')
Mask=float(copy(qmaxGrid));
Mask[qmaxGrid .>= qmaxMask] .= 1.0;
Mask[qmaxGrid .< qmaxMask] .= NaN;

#VeloUpGrid=Matrix(VeloUpGrid');

#GU.write_to_netcdf("./VeloEast.nc",Lon,Lat,VeloEastGrid);
#GU.write_to_netcdf("./VeloNorth.nc",Lon,Lat,VeloNorthGrid);
#GU.write_to_netcdf("./NetCDFFiles/VeloUp.nc",Lon,Lat,VeloUpGrid);

#------------------------------------------------

strainRateResults=CSV.read(strain_result_path, delim=' ', header=true, ignorerepeated=true, DataFrame)

Exx=strainRateResults.Dphph;
Eyy=strainRateResults.Dthth;
Exy=-strainRateResults.Dthph; #!!! The minus sign is because in Tape notation the versor for latitude points towards south

ExxGrid=reshape(Exx, length(Lat), length(Lon));
EyyGrid=reshape(Eyy, length(Lat), length(Lon));
ExyGrid=reshape(Exy, length(Lat), length(Lon));


#-- Return everything

LonGrid=Matrix(LonGrid')
LatGrid=Matrix(LatGrid')
VeloEastGrid=Matrix(VeloEastGrid')
VeloNorthGrid=Matrix(VeloNorthGrid')
ExxGrid=Matrix(ExxGrid')
EyyGrid=Matrix(EyyGrid')
ExyGrid=Matrix(ExyGrid')


return (
    LonGrid = LonGrid,
    LatGrid = LatGrid,
    Lon = Lon,
    Lat= Lat,
    VeloEastGrid = VeloEastGrid,
    VeloNorthGrid = VeloNorthGrid,
    ExxGrid = ExxGrid,
    EyyGrid = EyyGrid,
    ExyGrid = ExyGrid,
    Mask = Mask,
    ResLon = ResLon,
    ResLat = ResLat,
    ResE = ResE,
    ResN = ResN,
    qmaxGrid=qmaxGrid
)

end
