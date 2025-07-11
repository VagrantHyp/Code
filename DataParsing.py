import numpy as np
import matplotlib.pylab as plt
from scipy import stats
from scipy.integrate import simpson  # for numerical integration

#Modify this everytime new temp is tested for any configuration
#This is for assigning color to each temperature tested
temps = [0.1, 0.3, 0.4, 0.5, 0.75, 1, 2, 5]

#Create an object that represents simulation tested under a constant temperature
class CGT:

    def __init__(self, f, t, Nseg, ReePaths, RgPaths):
        
        self.f = np.array(f)
        self.t = t
        self.Nseg = Nseg
        self.color = 'C' + str(temps.index(t))

        self.RgPaths = RgPaths
        self.ReePaths = ReePaths
        
        #Check if the input force and data are consistent
        if len(f) != len(ReePaths):
            raise ValueError("Length of force and data doesn't match")
        
        #Global Ree parameter
        self.Nlpf = 10009
        self.M = 100; self.N = 100
        self.Ndata = int(self.M*self.N)
        self.Ncols = 9
        self.Mseg = self.N//self.Nseg * self.M 
        
        #Data in equilibrium
        self.ParsedRee = [] #Nf x M*N (10000) x col
        self.ParsedRg = [] #Nf x M (100)
        
        #Averaged data per time frame
        self.RgData = [] #Nf x mean Rg
        self.ReeData = [] #Nf x mean Ree
        
        #Mean extension over all timesteps
        self.aveRg = np.empty(0) 
        self.aveRee = np.empty(0)
        
        #Draw both Rg and Ree vs time
        fig, axes = plt.subplots(nrows = 1, ncols= 2, figsize = (20, 8))
        axes[0].set_title('Rg')
        axes[1].set_title('Ree')
        
        for i in range(len(self.RgPaths)):

            #Rg
            t, data = self.RgParse(self.RgPaths[i])
            Rg = data.mean(axis = 1) #Mean at each tf
            
            #Find equilibrium time
            eqTf = self.findEquil(t[1:],Rg[1:]) #Exclude first frame
            if eqTf == t[-1]:
                print(f'No equilirbium found: F{self.f[i]}')
            mask = (t>=eqTf)
            self.ParsedRg.append(data[mask])
            self.RgData.append([t[mask], Rg[mask]])
            self.aveRg = np.append(self.aveRg, Rg[mask].mean())
            
            axes[0].plot(t, Rg, label = f'F{self.f[i]}')
            axes[0].plot(t[mask], Rg[mask], color = 'black')
            
            #Ree
            t, data = self.ReeParse(self.ReePaths[i])
            Ree = self.ReeCalc(data)
            
            #Find equilibrium time
            eqTf = self.findEquil(t[1:],Ree[1:]) #Exclude first frame
            if eqTf == t[-1]:
                print(f'No equilirbium found: F{self.f[i]}')
            mask = (t>=eqTf)
            self.ParsedRee.append(data[mask])
            self.ReeData.append([t[mask], Ree[mask]])
            self.aveRee = np.append(self.aveRee, Ree[mask].mean())
            
            axes[1].plot(t, Ree, label = f'F{self.f[i]}')
            axes[1].plot(t[mask], Ree[mask], color = 'black')
        
        plt.tight_layout()
        plt.legend()
        plt.show()
        
    #Parse Rg data into workable matrix
    def RgParse(self, path):
        Nlpf = 101
        D = np.loadtxt(path)
        Nf = D.shape[0]//101
        D = D.reshape(Nf,Nlpf,2)
        
        #Get timeframe from first row
        t = D[:,0,0]/1e5
        t = t - t[0] #Start from 0
        #Get Rg from second column
        Rg = D[:,1:,1]
        return t,Rg
    
    #Parse Ree data into workable matrix
    def ReeParse(self, path): 
        with open(path, 'r') as file:
            
            LINES = file.readlines()
            
            Nf = len(LINES) // (10009)

            file.seek(0) # go back to start of file

            all_frames=[]

            while True:
                # Skip the non-numeric header lines (assumed to be 9 lines)
                for _ in range(9):
                    line = file.readline()
                    if not line:
                        break  # End of file
                else:
                    # Read the next block of data
                    frame_data = np.loadtxt(file, max_rows=self.Nlpf-9)
                    all_frames.append(frame_data)
                    continue
                break  # Exit the loop if EOF
        
        # Combine all frames into a single numpy array
        D = np.vstack(all_frames).reshape(Nf,self.Ndata,self.Ncols)
        
        #Sort each frame by atom type
        for f in range(Nf):
            ids = D[f,:,0]
            key = np.argsort(ids)
            D[f,:,:] = D[f,key,:]

        t = np.linspace(0, (Nf-1)*3, Nf)
        
        return t, D
        
    #Insert parsed Ree data and calculate Ree
    def ReeCalc(self, D):
        Nf = D.shape[0]
        
        #Extract all atoms that are type 3
        Type3 = D[D[:,:,2] == 3]
        Type3 = Type3.reshape(Nf, 100, 9)

        #Extract all atoms that are type 4
        Type4 = D[D[:,:,2] == 4]
        Type4 = Type4.reshape(Nf, 100, 9)

        #Extract Position vector
        Type3 = Type3[:,:,3:6]
        Type4 = Type4[:,:,3:6]

        #Extract End-to-end vector
        dX = Type3-Type4
        Ree = np.linalg.norm(dX, axis = 2).mean(axis=1)
        
        return Ree
    
    # Equilibrium criterion:
    # Linear regression of average RG data less than cutoff = 0.01
    # Timescale normalized about the correlation time
    # Source: https://chem.libretexts.org/Bookshelves/Biological_Chemistry/Concepts_in_Biophysical_Chemistry_(Tokmakoff)/06%3A_Dynamics_and_Kinetics/22%3A_Biophysical_Reaction_Dynamics/22.05%3A_Time-Correlation_Functions
    def findEquil(self, t, Rg):
        cutoff = 0.05

        #Compute autocorrelation
        dRg = Rg - Rg.mean()
        result = np.correlate(dRg, dRg, mode = 'full')
        result = result[result.size // 2:] # Keep the positive lags
        acf = result / result[0] #Normalized autocorrelation
        
        #Apply numerical integration
        positive_acf = acf[acf > 0]  # Ignore negative part to avoid integration errors
        tau_c = simpson(positive_acf, dx=1) #Correlation time
        t_scaled = t/tau_c
        
        for i in range(len(t)-int(tau_c)):
            res = stats.linregress(t_scaled[i:], Rg[i:])
            if res.pvalue >= cutoff:
                return t[i]
        
        return t[-1] #If no equilibrium found, return last frame
