function LabelFocalMechanisms(Mrr,Mtt,Mff,Mrt,Mrf,Mtf)

# In GMT they use r, t, f system that is up, south, east (Harvard/Global CMT convention)
# Same convention for GCMT for the six moment-tensor elements (p=f): Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, where r is up, t is south, and p is east. 

fthrust,fstrikeslip,fnormal=FocalMechanismsType(Mrr,Mtt,Mff,Mrt,Mrf,Mtf)

condThrust=(fthrust .> fstrikeslip) .& (fthrust .> fnormal);
condStrikeSlip=(fstrikeslip .> fthrust) .& (fstrikeslip .> fnormal);
condNormal=(fnormal .> fthrust) .& (fnormal .> fstrikeslip);

# it is possible (in principle) that the FM is perfectly in between more conditions!

return condThrust, condStrikeSlip, condNormal

end

"""
PrincipalAxesFM

Get the eigenvalues and eigenvector decomposition of the moment tensor.
Expected vectors are Mrr,Mtt,Mff,Mrt,Mrf,Mtf where the convention followed is r=up, t=south, and f=east; representing the components of a vector of moment tensors.
Do not froget that the exp.

The output reference system is ENU
Eigenvalues refer to the FULL moment tensor M as provided in input
P-N-T (Pressure-Null-Tension) vectors are the (orthonormal) eigenvectors corresponding to the min, intermediate and max eigenvector, respectively 
For each Moment Tensor the isotropic part is also returned: +(1/3)*trace(M)
"""
function PrincipalAxesFM(Mrr_in,Mtt_in,Mff_in,Mrt_in,Mrf_in,Mtf_in;exp=ones(length(Mrr_in))) 

    Mrr=Mrr_in .* (10.0 .^ (exp));
    Mtt=Mtt_in .* (10.0 .^ (exp));
    Mff=Mff_in .* (10.0 .^ (exp));
    Mrt=Mrt_in .* (10.0 .^ (exp));
    Mrf=Mrf_in .* (10.0 .^ (exp));
    Mtf=Mtf_in .* (10.0 .^ (exp));
    # Convention is r is up, t (s) is south, and f (e) is east.
    # Results are given in ENU (East-North-Up) convention

    MinEigenValues=[]
    InterEigenValues=[]
    MaxEigenValues=[]

    IsoPart=[]

    Pvectors=zeros(length(Mtt),3)
    Nvectors=zeros(length(Mtt),3)
    Tvectors=zeros(length(Mtt),3)

    NumberOfMomentTensors=length(Mtt)

    for i=1:NumberOfMomentTensors

        mee=Mff[i]
        mnn=Mtt[i]
        muu=Mrr[i]
        
        men=-Mtf[i]
        meu=Mrf[i]
        mnu=-Mrt[i]

        # remove the isotropic part (deal with double couple only)
        E=(mee + mnn + muu)/3

        Mmatrix=[mee men meu;
                 men mnn mnu;
                 meu mnu muu]
        
        # compute eigenvalues and eigenvectors
        F = eigen(Mmatrix)
        # F contains the eigenvalues in F.values and the orthonormal eigenvectors in the COLUMNS of the matrix F.vectors
        # (The kth eigenvector can be obtained from the slice F.vectors[:, k].)
        Indices=sortperm(F.values)
        EigenvaluesOrdered=F.values[Indices]
        EigenvectorsOrdered=F.vectors[:, Indices]
    
        MinEigenValue=EigenvaluesOrdered[1]
        InterEigenValue=EigenvaluesOrdered[2]
        MaxEigenValue=EigenvaluesOrdered[3]
    
        # Eigenvectors
        Pvector=EigenvectorsOrdered[:,1]
        Nvector=EigenvectorsOrdered[:,2]
        Tvector=EigenvectorsOrdered[:,3]

        push!(IsoPart,E)
        push!(MinEigenValues,MinEigenValue)
        push!(InterEigenValues,InterEigenValue)
        push!(MaxEigenValues,MaxEigenValue)
    
        Pvectors[i,:]=[Pvector[1] Pvector[2] Pvector[3]] #Insert the eigenvectors in the ROWS
        Nvectors[i,:]=[Nvector[1] Nvector[2] Nvector[3]]
        Tvectors[i,:]=[Tvector[1] Tvector[2] Tvector[3]]

    end


    return MinEigenValues,InterEigenValues,MaxEigenValues,IsoPart,Pvectors,Nvectors,Tvectors

end


"""
HorPrincipalAxesFM

Same of PrincipalAxesFM but using horizontal M tensor
"""
function HorPrincipalAxesFM(Mtt_in,Mff_in,Mtf_in;exp=ones(length(Mtt_in))) #Same of PrincipalAxesFM but 2D

    Mtt=Mtt_in .* (10.0 .^ (exp));
    Mff=Mff_in .* (10.0 .^ (exp));
    Mtf=Mtf_in .* (10.0 .^ (exp));

    MinEigenValues=[]
    MaxEigenValues=[]

    IsoPart=[]

    Pvectors=zeros(length(Mtt),2)
    Tvectors=zeros(length(Mtt),2)

    NumberOfMomentTensors=length(Mtt)

    for i=1:NumberOfMomentTensors

        mee=Mff[i]
        mnn=Mtt[i]  
        men=-Mtf[i]

        # remove the isotropic part (deal with double couple only)
        E=(mee + mnn)/2

        Mmatrix=[mee men;
                 men mnn]
        
        # compute eigenvalues and eigenvectors
        F = eigen(Mmatrix)
        Indices=sortperm(F.values)
        EigenvaluesOrdered=F.values[Indices]
        EigenvectorsOrdered=F.vectors[:, Indices]
    
        MinEigenValue=EigenvaluesOrdered[1]
        MaxEigenValue=EigenvaluesOrdered[2]
    
        # Eigenvectors
        Pvector=EigenvectorsOrdered[:,1]
        Tvector=EigenvectorsOrdered[:,2]

        push!(IsoPart,E)
        push!(MinEigenValues,MinEigenValue)
        push!(MaxEigenValues,MaxEigenValue)
    
        Pvectors[i,:]=[Pvector[1] Pvector[2]]
        Tvectors[i,:]=[Tvector[1] Tvector[2]]

    end


    return MinEigenValues,MaxEigenValues,IsoPart,Pvectors,Tvectors

end

"""
HorPrincipalAxesFM

Compute the projection of principal axes on hor. plane (FM principal dirs. as seen from above)
Return the vector used to plot the projected principal axes
"""

function ProjectedPrincipalAxes(Mrr_in,Mtt_in,Mff_in,Mrt_in,Mrf_in,Mtf_in;exp=ones(length(Mrr_in)),normalize=true)

    #Compute principal axes first:
    MinEigenValues,InterEigenValues,MaxEigenValues,_,Pvectors,Nvectors,Tvectors=PrincipalAxesFM(Mrr_in,Mtt_in,Mff_in,Mrt_in,Mrf_in,Mtf_in;exp=exp) 

    P_E=Pvectors[:,1] .* MinEigenValues
    P_N=Pvectors[:,2] .* MinEigenValues
    N_E=Nvectors[:,1] .* InterEigenValues
    N_N=Nvectors[:,2] .* InterEigenValues
    T_E=Tvectors[:,1] .* MaxEigenValues
    T_N=Tvectors[:,2] .* MaxEigenValues

    if(normalize)
        my_norms=sqrt.(MinEigenValues.^2 + InterEigenValues.^2 .+ MaxEigenValues.^2); #Frobenius norm
    else
        my_norms=ones(length(MinEigenValues))
    end

    P_E_nor=P_E ./ my_norms
    P_N_nor=P_N ./ my_norms
    N_E_nor=N_E ./ my_norms
    N_N_nor=N_N ./ my_norms
    T_E_nor=T_E ./ my_norms
    T_N_nor=T_N ./ my_norms

    return P_E_nor,P_N_nor,N_E_nor,N_N_nor,T_E_nor,T_N_nor


end

"""
FocalMechanismsType

Evaluate the type of focal mechanism in the sense proposed by Frohlich and Apperson (1992). No need for exp
"""
function FocalMechanismsType(Mrr,Mtt,Mff,Mrt,Mrf,Mtf)

    _,_,_,_,Pvectors,Nvectors,Tvectors=PrincipalAxesFM(Mrr,Mtt,Mff,Mrt,Mrf,Mtf) 

    fthrust=zeros(length(Mrr))
    fstrikeslip=zeros(length(Mrr))
    fnormal=zeros(length(Mrr))

    NumberOfMomentTensors=length(Mtt)

    for i=1:NumberOfMomentTensors

        # Eigenvectors
        Pvector=Pvectors[i,:]
        Nvector=Nvectors[i,:]
        Tvector=Tvectors[i,:]

        #compute the angle of each eigenvector w.r.t. the horizontal plane
        Pangle=atan(Pvector[3],sqrt(Pvector[1]^2+Pvector[2]^2))
        if(Pangle<0)
           Pangle=Pangle+pi
        end
        Interangle=atan(Nvector[3],sqrt(Nvector[1]^2+Nvector[2]^2))
        if(Interangle<0)
           Interangle=Interangle+pi
        end
        Tangle=atan(Tvector[3],sqrt(Tvector[1]^2+Tvector[2]^2))
        if(Tangle<0)
            Tangle=Tangle+pi
        end

        fthrust[i]=sin(Tangle)^2
        fstrikeslip[i]=sin(Interangle)^2
        fnormal[i]=sin(Pangle)^2

    end

    return fthrust,fstrikeslip,fnormal

end

function FocalMechanismsColor(Mrr,Mtt,Mff,Mrt,Mrf,Mtf; colors=[1 0 0; 0 1 0 ; 0 0 1]) #default colors: red,green,blue

    fthrust,fstrikeslip,fnormal=FocalMechanismsType(Mrr,Mtt,Mff,Mrt,Mrf,Mtf)

    x=sqrt.(fthrust)
    y=sqrt.(fstrikeslip)
    z=sqrt.(fnormal)

    rgbcolors, _, _, _=map_ternary_color_octant(x,y,z,colors[:,1],colors[:,2],colors[:,3])
    hex_colors=rgb_to_hex(rgbcolors)

    return hex_colors

end

function MomentTensorsSum(Mrr,Mtt,Mff,Mrt,Mrf,Mtf,exps)

    Mrrsum=sum(Mrr.*(10.0 .^ exps))
    Mttsum=sum(Mtt.*(10.0 .^ exps))
    Mffsum=sum(Mff.*(10.0 .^ exps))
    Mrtsum=sum(Mrt.*(10.0 .^ exps))
    Mrfsum=sum(Mrf.*(10.0 .^ exps))
    Mtfsum=sum(Mtf.*(10.0 .^ exps))

    MSumVec=[Mrrsum,Mttsum,Mffsum,Mrtsum,Mrfsum,Mtfsum]

    # Normalize in some way
    MaxValue=maximum(MSumVec)
    expFinal=floor(log10(MaxValue))

    MsumFinal=MSumVec ./ (10^(expFinal))

    return MsumFinal[1],MsumFinal[2],MsumFinal[3],MsumFinal[4],MsumFinal[5],MsumFinal[6],expFinal

end

function SumFMInCells(Lon,Lat,Depth,Mrr,Mtt,Mff,Mrt,Mrf,Mtf,exps,cells,region) #use region for speed-up calculations

    FMData=[Lon Lat Depth Mrr Mtt Mff Mrt Mrf Mtf exps] 
    SumFM = Matrix{Float64}(undef, 0, size(FMData,2))

    indices_mapped=[]
    my_count=0

    #cells is expected to be a vector of matrices
    for (cell_indx, cell) in enumerate(cells)
        
        my_cell=copy(cell)
    
        test_lon=mean(my_cell[:,1])
        test_lat=mean(my_cell[:,2])
    
        if((test_lon >= region[1]) & (test_lon<=region[2]) & (test_lat >= region[3]) & (test_lat <= region[4]))
        
            polycontour=[my_cell[k,:] for k in axes(my_cell,1)];
            myPolygon=LibGEOS.Polygon([polycontour]);
            withinPoly=Bool[]

            for j in eachindex(Lon)
                onePoint=LibGEOS.Point(Lon[j],Lat[j]);
                push!(withinPoly,LibGEOS.intersects(onePoint,myPolygon));
            end
    
            myFM=FMData[withinPoly,:]

            if(!isempty(myFM))
                my_count=my_count+1
                Mrr,Mtt,Mff,Mrt,Mrf,Mtf,exp=MomentTensorsSum(myFM[:,4],myFM[:,5],myFM[:,6],myFM[:,7],myFM[:,8],myFM[:,9],myFM[:,10])
                SumFM=vcat(SumFM,[mean(my_cell[:,1]) mean(my_cell[:,2]) mean(myFM[:,3]) Mrr Mtt Mff Mrt Mrf Mtf exp])    
                push!(indices_mapped,cell_indx)
            end
        end
    end

    return SumFM,indices_mapped
end