import matplotlib.image as img
image = img.imread('C:\Users\Haris\Desktop\OR - Seminarski\Lenna.png') #input path of desired image
from matplotlib import pyplot as plt
import numpy as np
import random

m=0 #global coordinates
n=0
   
def MatrixConvolution(A,B):
    if len(A[:,1])!=len(B[:,1]) or len(A[1,:])!=len(B[1,:]):
        raise Exception("Two non-compatible matrix")
    
    R=len(A[:,1])
    C=len(B[1,:])
    sum1=0
    for j in range(0,C):
        for i in range(0,R):
            sum1+=A[i,j]*B[R-i-1,C-j-1]
              
    return sum1
                
    
def SobelOperator(image):
        
    R=len(image[:,1])
    C=len(image[1,:])
    np.ndarray(shape=(3,3), dtype=int, order='F')
    Ahorizontal=np.array([[ -1, 0, 1], [  -2, 0, 2], [-1, 0, 1],])
    Avertical=np.array([[ -1, -2, -1], [  0, 0, 0], [1, 2, 1],])
    Acopy=np.copy(image)
    for i in range(1,R-1):
        for j in range(1,C-1):
            B=np.array([[ image[i-1,j-1], image[i-1,j], image[i-1,j+1]], [  image[i,j-1], image[i,j], image[i,j+1]], [image[i+1,j-1], image[i+1,j], image[i+1,j+1]]])
            Gy=MatrixConvolution(Avertical,B)
            Gx=MatrixConvolution(Ahorizontal,B)
            Acopy[i,j]=np.sqrt(Gy*Gy+Gx*Gx)
                
    return Acopy

class ACO:
    tau_init=float
    N=int #number of iterations
    L=int #number of construction steps
    K=int #number of ants
    alpha=float
    beta=float
    phi=float #pheromone decay coefficient
    rho=float #pheromone evaporation coefficient
    treshold=float
    PheromoneMap=[] #map of pheromones
    InformationMatrix=[] #map of eta values (heuristic information)
    R=int #row length of Pheromone Map
    C=int #column length of Pheromone Map
    delta_tau_Matrix=[] #matrix that holds elements for delta_tau(i,j) calculation
    RoutteMatrix=[] #matrix of ants previous positions
    
    def DecayPheromones(self): #before global pheromone update we need to decay all of the pheromones
        for i in range(0, self.R):
            for j in range(0, self.C):
               self.PheromoneMap[i,j]=(1-self.rho)*self.PheromoneMap[i,j]
    
    def GlobalUpdate(self):
        for i in range(0, self.R):
            for j in range(0, self.C):
                self.PheromoneMap[i,j] += self.delta_tau_Matrix[i][j]
    
    def LocalUpdate(self, pixel , k, l):
        
        i0=pixel[0]
        j0=pixel[1]
        xMinLim=i0-1
        xMaxLim=i0+1
        yMinLim=j0-1
        yMaxLim=j0+1
        
        if i0==0: #if ant reaches final edge of the image, neighbourhood mustn't consist of those coordinates
            xMinLim==0
        if j0==0:
            yMinLim==0
        if i0>=self.R-1:
            xMaxLim=self.R-1
        if j0>=self.C-1:
            yMaxLim=self.C-1
        
        neighbourhood=[]
       
        for i in range(xMinLim, xMaxLim+1): #creating the neighbourhood
            for j in range(yMinLim, yMaxLim+1):
                if (i!=i0 or j!=j0):
                    neighbourhood.append([i,j])
        
        for node in neighbourhood: #removing those already visited
            for positions in self.RoutteMatrix[k]:
                if positions[0]==node[0] and positions[1]==node[1]:
                    neighbourhood.remove(node)
                        
                    break
            else:
                continue
            break
        
        u=random.random()
        p=0
        
        if not neighbourhood: #if neighbourhood is empty use the same pixel
            m=i0
            n=j0
            self.RoutteMatrix[k].append([m,n])
            
        else:
                  
            j=0
            brojac=0           
            while u>p: #testing the neighbourhood and using Roulette Wheel system to pick the new pixel
                
                p+=float(pow(self.PheromoneMap[neighbourhood[j][0],neighbourhood[j][1]][0],self.alpha))*float(pow(self.InformationMatrix[neighbourhood[j][0]][neighbourhood[j][1]], self.beta))
                
                j+=1
                brojac+=1
                if j==len(neighbourhood):
                    j=0
                    
                if brojac >15: #if p is increasing too slow we can randomly choose next pixel
                    j=random.randint(0, len(neighbourhood))
                    break
                
            self.RoutteMatrix[k].append(neighbourhood[j-1])
            m=neighbourhood[j-1][0]
            n=neighbourhood[j-1][1]
               
        self.PheromoneMap[m,n]=(1-self.phi)*self.PheromoneMap[m,n]+self.phi*self.tau_init #local pheromone update
        self.delta_tau_Matrix[m][n]+=self.InformationMatrix[m][n]/float(l+1) #we used heuristic information eta to eliminate the noise which is caused by initial random positioning of ants
        
        
    def Run(self):
        for n in range(0, self.N): #number of iterations
            for l in range(0, self.L): #number of ant steps
                for k in range(0, self.K): #number of ants
                        
                    self.LocalUpdate(self.RoutteMatrix[k][l], k, l)
                                           
            self.DecayPheromones()
            self.GlobalUpdate()
            
    def ShowImage(self):
        FinalImage=[]
        
        for i in range(0, self.R):
            lista=[]
            for j in range(0, self.C):
                if self.PheromoneMap[i,j][0]<self.tau_init:
                    lista.append(0)
                    
                else:
                    lista.append(1)
            FinalImage.append(lista)
                        
        plt.imshow(FinalImage, interpolation='nearest')
        plt.show()
    
    
    def __init__(self, image, tau_init=0.1, N=2, L=50, K=5000, alpha=1.0, beta=2.0, phi=0.05, rho=0.1, treshold=0.6):
        if type(tau_init)!=float:
            raise Exception ("Tau_init must be float!")
        self.tau_init=tau_init
        if type(N)!=int:
            raise Exception ("Number of iterations N must be int!")
        self.N=N
        if type(L)!=int:
            raise Exception ("Number of ant steps L must be int!")
        self.L=L
        if type(L)!=int:
            raise Exception ("Number of ants K must be int!")
        self.K=K
        if type(alpha)!=float:
            raise Exception ("Alpha must be float!")        
        self.alpha=alpha
        if type(beta)!=float:
            raise Exception ("Beta must be float!")
        self.beta=beta
        if type(phi)!=float:
            raise Exception ("Phi must be float!")
        self.phi=phi
        if type(rho)!=float:
            raise Exception ("Rho must be float!")
        self.rho=rho
        if type(treshold)!=float:
            raise Exception ("Treshold must be float!")
        self.treshold=treshold
        self.R=len(image[:,1])
        self.C=len(image[1,:])
        if len(image[0,0])<3:
            raise Exception ("Image format must be png!")
        self.PheromoneMap=np.copy(image)
                
        for i in range(0,self.R):
            for j in range(0,self.C):
                self.PheromoneMap[i, j]=self.tau_init    
                
                        
        for i in range(0,self.R):
            lista1=[]
            for j in range(0,self.C):
                lista1.append(0)
            
            self.delta_tau_Matrix.append(lista1)
        
        SobelMatrix=SobelOperator(image)         
        for i in range(0,self.R):
            lista=[]
            for j in range(0,self.C):
                if SobelMatrix[i,j][0]*np.sqrt(len(SobelMatrix[i,j]))>self.treshold: #treshold
                    lista.append(2)
                    
                else:
                    lista.append(0)
                    
            self.InformationMatrix.append(lista)
        
        for i in range(0, self.K): #creating a matrix of all ant positions
            
            self.RoutteMatrix.append([])
            self.RoutteMatrix[i].append([random.randint(1, self.R-1), random.randint(1, self.C-1)])
      
slika=ACO(image)  
slika.Run()          
slika.ShowImage()



