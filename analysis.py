import numpy as np
from DataParsing import CGT
from scipy.integrate import simpson

#Bond length
b = 1

#Get segmental end-to-end distance

def segRee(D: np, Blck: int, Atom_Type: str):
    """
    Parse the hydrophobic or hydrophilic data
    
    D: data from the dump file
    Blck: length of block (1,5,10,50)
    Atom_Type: either "Hydrophobic" or "Hydrophilic"

    return the end-to-end distance of each hydrophobic/hydrophilic segment
    """
    #Ree parameters
    M = 100 # Number of molecules per simulation
    N = 100 # Number of atoms per molecule

    Nf = D.shape[0]
    Mseg = N//Blck*M

    # Reshape by Nf x Mseg x Blck x 9       
    D = D[1:].reshape(Nf-1,Mseg,Blck,9) #exclude first frame

    Xseg = D[:,:,:,3:6]

    if(Atom_Type == 'Hydrophilic'):
        # Get just the odd segments (hydrophilic segments)
        Xhp =  Xseg[:,1::2,:,:]
    elif(Atom_Type == 'Hydrophobic'):
        # Get just the even segments (hydrophobic segments)
        Xhp =  Xseg[:,0::2,:,:]
    
    # End-to-End vectors of the segments
    dX = Xhp[:,:,-1,:]-Xhp[:,:,0,:]
    Rhp = np.linalg.norm(dX,axis=2)
    return Rhp

def ComputeTauC(mRee: np):
    """
    Compute correlation time from the mean configurational trajectory
    """
    dRee = mRee - mRee.mean()
    result = np.correlate(dRee, dRee, mode = 'full')
    result = result[result.size // 2:] # Keep the positive lags
    acf = result / result[0] #Normalized autocorrelation
    
    #Apply numerical integration
    positive_acf = acf[acf > 0]  # Ignore negative part to avoid integration errors
    tau_c = simpson(positive_acf, dx=1) #Correlation time
    return tau_c
# end def

def mbar(inputs, error = 1e-6, iteration = int(1e4)):

    """
    Reference: https://doi.org/10.1103/PhysRevLett.63.1195
    Compute the unbiased potential energy distribution using k umbrella simulations
        P0 = N(U)*exp(-beta*U), N(U) is the density of states
    f_n: 1 x k
        Free energy difference
    cTau: 1 x k
        Correlation energy
    eta: k x bin
        Biasing term
        Constant force and temperature: -f*R
        Umbrella Sampling that enforces R0 through a spring with strength k: k*(R-R0)
    N_n: 1 x bin
        Number of observation of a particular order parameter R
    n_m: 1 x k
        Number of observation in each umbrella simulation
    T: float
        Temperature of the simulation
    """

    #Start by assuming free energy difference is 0
    f_n = np.zeros(inputs[0].shape[0])

    for i in range(iteration):
        f_n_new = iterate(inputs, f_n)
        diff = np.abs(f_n-f_n_new)
        if diff.max() < error:
            print(f'Convergence is achieved after {i} iterations')
            return f_n-f_n[0]
        # end if
        f_n = f_n_new
    # end for

    f_n = f_n-f_n[0] #Reference with respect to the extended coil
    print(f'Did not converge with error < {error}')
    print(f'Difference from the last iteration:{diff}')
    return f_n
# end def


def iterate(inputs, f_n):
    cTau, eta, N_n, n_m, T = inputs

    # The bias factor
    beta = 1/T # Energy scaling 1/kT
    bias_factor = np.exp(-beta*eta)

    g_n = 1+2*cTau
    numerator = 1/g_n @ N_n
    exponent = -beta*(eta-f_n[:, np.newaxis])
    denominator = n_m*1/g_n@np.exp(exponent)
    P = (numerator / denominator)[np.newaxis, :]*bias_factor
    sum = np.sum(P, axis = 1)
    f_n = np.log(sum)/-beta

    return f_n
# end def