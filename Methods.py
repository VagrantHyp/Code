import numpy as np
import matplotlib.pylab as plt

#Modify this everytime new temp is tested for any configuration
temps = [0.1, 0.3, 0.5, 0.75, 1, 2, 5]

class CGT:

    def __init__(self, f, t, RgPaths, ReePaths):
        
        self.f = f
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
        
        #Average data at each force
        self.aveRg = self.mean('Rg')
        self.aveRee = self.mean('Ree')
        
    #Rg Methods
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

    #Ree Methods
    def ReeCalc(self, path): 
        Nlpf = 10009
        M = 100; N = 100
        Ndata = int(M*N)
        Ncols = 9

        with open(path, 'r') as file:
            LINES = file.readlines()

        D = np.loadtxt(path, skiprows=9, max_rows=10000)
        Nf = len(LINES) // Ndata
        Nfi = int(LINES[1])
        for i in range(Nf)[1:]:
            D = np.append(D,np.loadtxt(path, skiprows = 9+10009*i, max_rows = 10000), axis = 0)
        D = D.reshape(Nf,Nlpf-9,Ncols)
        
        #Calculate end-to-end
        ri = np.zeros((Nf,M,Ncols))
        rf = np.zeros((Nf,M,Ncols))
        Ree = np.zeros((Nf,M))

        #Find all data with type 3 or 4 molecule (first and last molecule in a chain)
        for t in range(Nf):
            ji = 0
            jf = 0
            for i in range(Ndata):
                if D[t,i,2] == 3:
                    ri[t,ji] = D[t,i,:]
                    ji += 1
                elif D[t,i,2] == 4:
                    rf[t,jf] = D[t,i,:]
                    jf += 1
            #Sort according to id
            ri[t] = ri[t][ri[t][:,1].argsort()]
            rf[t] = rf[t][rf[t][:,1].argsort()]


        #Extact position
        ri = ri[:,:,3:6]
        rf = rf[:,:,3:6]

        #Calculation of distance
        for t in range(Nf):
            for i in range(M):
                dx = rf[t,i,0] - ri[t,i,0]
                dy = rf[t,i,1] - ri[t,i,1]
                dz = rf[t,i,2] - ri[t,i,2]
                Ree[t,i] = np.sqrt(dx**2+dy**2+dz**2)

        #Return the Nf vector and the average
        nf = np.linspace(Nfi,(Nfi//300000+Nf-1)*300000,Nf)
        Ree = Ree.mean(axis = 1)
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
                equil = self.Equil(data[1])
                plt.plot(data[0][equil:], data[1][equil:], color = 'k')

        plt.legend()
        plt.xlabel("Timeframe")
        plt.ylabel('<'+ DataType + '>')
        plt.title('T = ' + str(self.t))

    def Equil(self, data):
        #Start searching for equilibirum halfway
        half = len(data) // 2
        for i in range(len(data))[half:]:
                if np.abs(data[i] - data[i-1]) < 0.1:
                    return i
        print("Use temporary mean starting halfway")
        return half

    #Calculate ave Rg/Ree at each f
    def mean(self, Datatype):
        Mean = []
        data = self.RgData
        if(Datatype == 'Ree'):
            data = self.ReeData

        for i in range(len(self.f)):
            equil = self.Equil(data[i][1])
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
    

#Plot at specific force
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
        equil = cgt.Equil(data[1])
        plt.plot(data[0][equil:], data[1][equil:], color = 'k')

        
    plt.legend()
    plt.xlabel('Tf')
    plt.ylabel('<'+ DataType +'>')
    plt.title('F = ' + str(f))