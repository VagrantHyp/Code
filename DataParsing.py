import numpy as np
import matplotlib.pylab as plt

#Modify this everytime new temp is tested for any configuration
#This is for assigning color to each temperature tested
temps = [0.1, 0.3, 0.4, 0.5, 0.75, 1, 2, 5]

#Create an object that represents simulation tested under a constant temperature
class CGT:

    def __init__(self, f, t, Nseg, ReePaths, RgPaths=None):
        
        self.f = f
        self.t = t
        self.Nseg = Nseg
        self.color = 'C' + str(temps.index(t))

        self.RgPaths = RgPaths if RgPaths is not None else ""
        self.ReePaths = ReePaths
        
        #Global Ree parameter
        self.Nlpf = 10009
        self.M = 100; self.N = 100
        self.Ndata = int(self.M*self.N)
        self.Ncols = 9
        self.Mseg = self.N//self.Nseg * self.M 
        
        #Raw data
        self.ParsedRee = [] #Nf x M*N (10000) x col
        self.ParsedRg = [] #Nf x M (100)
        
        #Averaged data per time frame
        self.RgData = [] #Nf x mean Rg
        self.ReeData = [] #Nf x mean Ree
        
        for p in self.RgPaths:
            t, data = self.RgParse(p)
            self.ParsedRg.append(data)
            self.RgData.append([t, data.mean(axis = 1)])
        
        for p in self.ReePaths:
            data = self.ReeParse(p)
            self.ParsedRee.append(data)
            self.ReeData.append(self.ReeCalc(data))
        
        #Averaged data per force
        self.aveRg = np.array(self.mean('Rg')) if RgPaths is not None else []
        self.aveRee = np.array(self.mean('Ree'))
        
    #Parse Rg data into workable matrix
    def RgParse(self, path):
        Nlpf = 101
        D = np.loadtxt(path)
        Nf = D.shape[0]//101
        D = D.reshape(Nf,Nlpf,2)
        
        #Get time from first row
        t = D[:,0,0]
        #Get Rg from second column
        Rg = D[:,1:,1]
        return t,Rg
    
    #Parse Ree data into workable matrix
    def ReeParse(self, path): 
        with open(path, 'r') as file:

            LINES = file.readlines()

            Nf = len(LINES) // self.Ndata

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

        return D

    #Insert parsed Ree data and calculate Ree
    def ReeCalc(self, D):
        Nf = D.shape[0]
        
        #Extract all atoms that are type 3
        Type3 = D[D[:,:,2] == 3]
        Type3 = Type3.reshape(Nf, 100, 9)

        #Extract all atoms that are type 4
        Type4 = D[D[:,:,2] == 4]
        Type4 = Type4.reshape(Nf, 100, 9)

        # Sort by molecular id
        for t in range(Nf):
            ids3 = Type3[t, :, 1]
            key3 = np.argsort(ids3)
            Type3[t,:,:] = Type3[t, key3, :]
            
            ids4 = Type4[t, :, 1]
            key4 = np.argsort(ids4)
            Type4[t,:,:] = Type4[t, key4, :]

        #Extract Position vector
        Type3 = Type3[:,:,3:6]
        Type4 = Type4[:,:,3:6]

        #Extract End-to-end vector
        dX = Type3-Type4
        Ree = np.linalg.norm(dX, axis = 2).mean(axis=1)
        
        #Create a timeframe array
        nf = np.linspace(0,(Nf-1)*300000,Nf)
        
        return [nf,Ree]
    
    #Plot tf vs extension under constant temperature
    #DataType must be either "Rg" or "Ree"
    def Plot(self, DataType):
        default = self.RgData
        if(DataType == 'Ree'):
            default = self.ReeData
        
        for i in range(len(self.f)):
                data = default[i]
                plt.plot(data[0],data[1], label = 'f=' + str(self.f[i]))
                
                #Plot Equilibirum Section
                equil = self.Equil(data[1], 'F = ' + str(self.f[i]))
                plt.plot(data[0][equil:], data[1][equil:], color = 'k')

        plt.legend()
        plt.xlabel("Timeframe")
        plt.ylabel('<'+ DataType + '> ($\sigma$)')
        plt.title('T = ' + str(self.t))

    #Find the timeframe at which equilibirum is achieved
    def Equil(self, data, str):
        #Start searching for equilibirum halfway
        half = len(data) // 2
        for i in range(len(data))[half:-1]:
                if np.abs(data[i] - data[-1]) < 0.25:
                    return i
        #If no i is returned, that means no equilibirum found. Calculate equilibirum from the second half of data
        print("Use temporary mean starting halfway:" + str)
        return half

    #Calculate ave Rg/Ree at each f
    def mean(self, Datatype):
        Mean = []
        data = self.RgData
        if(Datatype == 'Ree'):
            data = self.ReeData

        for i in range(len(self.f)):
            equil = self.Equil(data[i][1], 'F = ' + str(self.f[i]))
            Mean.append(data[i][1][equil:].mean())

        return Mean    

#Plot f if it is tested at at least two different t
def lsForces(lscgt):
    #Find the largest set of forces
    ls = []
    lsCGT = lscgt[:]
    force = lsCGT[0].f
    MTcgt = lsCGT[0] #Most Tested CGT
    for cgt in lsCGT[1:]:
        if(len(cgt.f) > len(force)):
            force = cgt.f
            MTcgt = cgt
    #Go through each force and see if it is tested at another temperature
    lsCGT.remove(MTcgt)
    for f in force:
        for cgt in lsCGT:
            if f in cgt.f:
                ls.append(f)
                break
    return ls
    

#Plot all temperatures at specific force
def Plot(lsCGT, f, DataType):
    for cgt in lsCGT:
        
        #Find index of f
        try:
            i = cgt.f.index(f)
        except ValueError:
            continue

        if(DataType == 'Rg'):
            data = cgt.RgData[i]
        else:
            data = cgt.ReeData[i]
        plt.plot(data[0], data[1], label = 'T = ' + str(cgt.t) + "$\epsilon$/$k_B$")
        
        # Plot equilibirum section
        equil = cgt.Equil(data[1], 'T = ' + str(cgt.t) + "$\epsilon$/$k_B$")
        plt.plot(data[0][equil:], data[1][equil:], color = 'k')

        
    plt.legend()
    plt.xlabel('Tf')
    plt.ylabel('<'+ DataType +'> ($\sigma$)')
    plt.title('F = ' + str(f) + '$\epsilon$/$\sigma$')