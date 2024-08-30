import numpy as np
import matplotlib.pylab as plt

#Modify this everytime new temp is tested for any configuration
#This is for assigning color to each temperature tested
temps = [0.1, 0.3, 0.4, 0.5, 0.75, 1, 2, 5]

#Class for creating an object that represents simulation tested under a constant temperature
class CGT:

    def __init__(self, f, t, RgPaths, ReePaths):
        
        self.f = np.array(f)
        self.t = t
        self.color = 'C' + str(temps.index(t))

        self.RgPaths = RgPaths
        self.ReePaths = ReePaths

        #Organize data into matrix
        self.RgData = []
        self.ReeData = []

        for p in self.RgPaths:
            self.RgData.append(self.Data(p))
        
        for p in self.ReePaths:
            self.ReeData.append(self.ReeCalc(p))
        
        #Average extension at each force
        self.aveRg = np.array(self.mean('Rg'))
        self.aveRee = np.array(self.mean('Ree'))
        
    #Parse Rg data into workable 3d matrix
    def Data(self, path):
        Nlpf = 101
        D = np.loadtxt(path)
        Nf = D.shape[0]//101
        D = D.reshape(Nf,Nlpf,2)
        #Get vector of frame times
        t = D[:,0,0]
        RG_All = D[:,1:,1]
        #Get Rg data for each frame and molecule
        Rg = RG_All.mean(axis=1)
        return [t,Rg]
    
    #Parse Ree data into workable 3d matrix
    def ReeCalc(self,path): 
        Nlpf = 10009
        M = 100; N = 100
        Ndata = int(M*N)
        Ncols = 9

        with open(path, 'r') as file:

            LINES = file.readlines()

            Nf = len(LINES) // Ndata
            Nfi = int(LINES[1])

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
                    frame_data = np.loadtxt(file, max_rows=Nlpf-9)
                    all_frames.append(frame_data)
                    continue
                break  # Exit the loop if EOF

        # Combine all frames into a single numpy array
        D = np.vstack(all_frames).reshape(Nf,Nlpf-9,Ncols)

        #Calculate end-to-end
        ri = np.zeros((Nf,M,Ncols))
        rf = np.zeros((Nf,M,Ncols))

        for t in range(Nf):
            frame = D[t]

            # Vectorized filtering
            type3_data = frame[frame[:Ndata, 2] == 3]  # Select first NData rows with type 3
            type4_data = frame[frame[:Ndata, 2] == 4]  # Select first NData rows with type 4

            # Sort by id (assuming id is in the first column, index 0)
            type3_data = type3_data[np.argsort(type3_data[:, 0])]
            type4_data = type4_data[np.argsort(type4_data[:, 0])]

            # Extract positions directly
            ri[t] = type3_data
            rf[t] = type4_data


        #Extact position
        ri = ri[:,:,3:6]
        rf = rf[:,:,3:6]

        #Calculate Ree and Distances 
        Ree = np.sqrt(np.square(rf-ri).sum(axis=-1)).mean(axis=1) 
        
        #Return the Nf vector and the average
        nf = np.linspace(Nfi,(Nfi//300000+Nf-1)*300000,Nf)

        return [nf,Ree]
    
    #Plot Data under constant temperature
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
        plt.ylabel('<'+ DataType + '>')
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
        plt.plot(data[0], data[1], label = 'T = ' + str(cgt.t))
        
        # Plot equilibirum section
        equil = cgt.Equil(data[1], 'T = ' + str(cgt.t))
        plt.plot(data[0][equil:], data[1][equil:], color = 'k')

        
    plt.legend()
    plt.xlabel('Tf')
    plt.ylabel('<'+ DataType +'>')
    plt.title('F = ' + str(f))