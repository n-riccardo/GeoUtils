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

function LabelFocalMechanisms(psmecafile)

    # Assuming -Sd format: lat_orig,lon_orig,depth,mrr,mss,mee,mrs,mre,mse,exp,x,y,ev_id


    FocalMecDataGMT=Matrix(FocalMecData[:,[2,1,3,4,5,6,7,8,9,10,11,12]]);
    cT, cS, cN =SAK.LabelFocalMechanisms(FocalMecDataGMT[:,4],FocalMecDataGMT[:,5],FocalMecDataGMT[:,6],
                    FocalMecDataGMT[:,7],FocalMecDataGMT[:,8],FocalMecDataGMT[:,9]);

    return condThrust, condStrikeSlip, condNormal

end

"""
PrincipalAxesFM

Get the eigenvalues and eigenvector decomposition of the moment tensor.
Expected vectors are Mrr,Mtt,Mff,Mrt,Mrf,Mtf where the convention followed is r=up, t=south, and f=east; representing the components of a vector of moment tensors

The output reference system is ENU
Eigenvalues refer to the FULL moment tensor M as provided in input
P-N-T (Pressure-Null-Tension) vectors are the (orthonormal) eigenvectors corresponding to the min, intermediate and max eigenvector, respectively 
For each Moment Tensor the isotropic part is also returned: +(1/3)*trace(M)
"""
function PrincipalAxesFM(Mrr,Mtt,Mff,Mrt,Mrf,Mtf) 

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
function HorPrincipalAxesFM(Mtt,Mff,Mtf) #Same of PrincipalAxesFM but 2D

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
FocalMechanismsType

Evaluate the type of focal mechanism in the sense proposed by Frohlich and Apperson (1992)
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