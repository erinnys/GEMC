# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import matplotlib.pyplot as plt

#Minimum Imae Convention
# sigma  angstrom
# epsilon  kJ/mol
#Boltzmann constant  kJ/mol⋅K

Ar={'sigma':3.401,'epsilon':0.978638}
He={'sigma':2.556,'epsilon':0.08368}
bolt=8.314462618e-3  

numbsit=[[269],[269]]
lab=['Ar']
lac=np.array([[[215.4434],[215.4434],[215.4434]],[[215.4434],[215.4434],[215.4434]]])
rad=2.5
b=8.314462618e-3*273
step=10000
sigma=[3.401]
epsilon=[0.978638]


class model:
    def __init__(self,latvec,num,rad,b,config,lab,sigma,epsilon):
        self.latvec=latvec
        self.num=num
        self.rad=rad
        self.b=b
        self.config=config
        self.tot=[0 for i in range(len(self.latvec))]
        for i in range(len(self.num)):
            for j in range(len(self.num[i])):
                self.tot[i]+=self.num[i][j]
        self.lab=lab
        self.perio='off'
        self.sigma=sigma
        self.epsilon=epsilon
        
    def dist(self,i,j,k,z,y):
        return round(np.linalg.norm(self.config[y][i].T[j]-self.config[y][k].T[z]),5)
    
    def initconfig(self):
        con=[[] for i in self.latvec]
        for i in range(len(self.latvec)):
            self.config.append([])
            con[i]=np.random.rand(len(self.latvec[i]),self.tot[i])*self.latvec[i]-self.latvec[i]/2
        for k in range(len(self.latvac)):
            ac=0
            for i in self.num[k]:
                self.config[k].append(con[k][:,ac:ac+i])
                ac=ac+i
        while self.checkconfig()==False:
            con=[[] for i in self.latvec]
            for i in range(len(self.latvec)):
                self.config.append([])
                con[i]=np.random.rand(len(self.latvec[i]),self.tot[i])*self.latvec[i]-self.latvec[i]/2
            for k in range(len(self.latvac)):
                ac=0
                for i in self.num[k]:
                    self.config[k].append(con[k][:,ac:ac+i])
                    ac=ac+i

            
    def checkconfig(self,index):
        chek=0
        for i in range(len(self.config[index])):
            for j in range(len(self.config[index][i].T)):
                for k in range(len(self.config[index])):
                    for z in range(len(self.config[index][k].T)):
                        if i==k and j==z:
                            pass
                        else:
                            if self.dist(index,i,j,k,z)<self.rad:
                                chek=chek+1
                            else:
                                pass
        if chek==0:
            return True
        else:
            return False
        
    def configE(self,index,index2,cutoff):
        if self.perio=='off':
            if index==1:
                utot=0
                for i in range(len(self.config[index2])):
                    for j in range(len(self.config[index2][i].T)):
                        for k in range(len(self.config[index2])):
                            for z in range(len(self.config[index2][k].T)):
                                if i==k and j==z:
                                    pass
                                elif self.dist(index2,i,j,k,z)>cutoff:
                                    pass
                                else:
                                    utot=utot+4*((epsilon[i]*epsilon[k])**0.5)*((((self.sigma[i]+self.sigma[k])/2)/self.dist(index2,i,j,k,z))**12-(((self.sigma[i]+self.sigma[k])/2)/self.dist(i,j,k,z))**6)
                return utot/2 
        else:
            if index==1:
                utot=0
                for i in range(len(self.config[index2])):
                    for j in range(len(self.config[index2][i].T)):
                        for k in range(len(self.config[index2])):
                            for z in range(len(self.config[index2][k].T)):
                                for a in range(3):
                                    for b in range(3):
                                        for c in range(3):
                                            if i==k and j==z and a==b==c==0:
                                                pass
                                            elif round(np.linalg.norm(self.config[index2][i].T[j]-self.config[index2][k].T[z])+np.array([a,b,c])*self.latvec,5)>cutoff:
                                                pass
                                            else:
                                                utot=utot+4*((epsilon[i]*epsilon[k])**0.5)*((((self.sigma[i]+self.sigma[k])/2)/round(np.linalg.norm(self.config[index2][i].T[j]-self.config[index2][k].T[z])+np.array([a,b,c])*self.latvec,5))**12-(((self.sigma[i]+self.sigma[k])/2)/round(np.linalg.norm(self.config[index2][i].T[j]-self.config[index2][k].T[z])+np.array([a,b,c])*self.latvec,5))**6)
                return utot/2 
            return 0
    #計算config能量
    # 1 : Lennard-Jones
    def volchange(self,length):
        cho=np.random.choice([1,-1])*length*round(np.random.rand(),5)
        oldE=self.configE(1,0)+self.configE(1,1)
        self.config[0]=self.config[0]*(1+cho)
        self.config[1]=self.config[1]*(1-cho)
        if 1<((self.latvec[0][0]*cho)**(3*self.tot[0])*(self.latvec[1][0]/cho)**(3*self.tot[1])/((self.latvec[0][0])**(3*self.tot[0])*(self.latvec[1][0])**(3*self.tot[1])))*np.exp(-self.b*((self.configE(1,0)+(self.configE(1,1))-oldE))):
            pass
        else:
            if np.random.random()<((self.latvec[0][0]*cho)**(3*self.tot[0])*(self.latvec[1][0]/cho)**(3*self.tot[1])/((self.latvec[0][0])**(3*self.tot[0])*(self.latvec[1][0])**(3*self.tot[1])))*np.exp(-self.b*((self.configE(1,0)+(self.configE(1,1))-oldE))):
                pass
            else:
                self.config[0]=self.config[0]/(1+cho)
                self.config[1]=self.config[1]/(1-cho)
    
    def parexchange(self):
        cho=np.random.choice([0,1])
        cho1=np.random.randint(0,len(self.config[cho]))
        cho2=np.random.randint(0,len(self.config[cho][cho1].T))
        temp1=self.config[0]
        temp2=self.config[1]
        oldE=self.configE(1,0)+self.configE(1,1)
        temp3=self.config[cho][cho1].T
        temp4=np.delete(temp3,[cho2],0)
        temp5=np.append(temp3[cho2],self.config[abs(1-cho)][cho1].T)
        self.config[cho][cho1]=temp4.T
        self.config[abs(1-cho)][cho1]=temp5.T
        if np.exp(-self.b*(self.configE(1,0)+(self.configE(1,1))-oldE))*self.tot[cho]*self.latvec[abs(1-cho)][0]**3/((self.latvec[cho][0]**3)*(self.tot[abs(1-cho)]+1))>1:
            self.tot[cho]=self.tot[cho]-1
            self.tot[abs(1-cho)]=self.tot[abs(1-cho)]+1
            self.num[cho][cho1]=self.num[cho][cho1]-1
            self.num[abs(1-cho)][cho1]=self.num[abs(1-cho)][cho1]+1
        else:
            if np.random.random()<np.exp(-self.b*(self.configE(1,0)+(self.configE(1,1))-oldE))*self.tot[cho]*self.latvec[abs(1-cho)][0]**3/((self.latvec[cho][0]**3)*(self.tot[abs(1-cho)]+1)):
                self.tot[cho]=self.tot[cho]-1
                self.tot[abs(1-cho)]=self.tot[abs(1-cho)]+1
                self.num[cho][cho1]=self.num[cho][cho1]-1
                self.num[abs(1-cho)][cho1]=self.num[abs(1-cho)][cho1]+1
            else:
                self.config[0]=temp1
                self.config[1]=temp2
        return 0
    
    
    def ranchoiceandmove(self,steplength):
        cho=np.random.choice([0,1])
        cho0=np.random.choice([1,-1])
        cho1=np.random.randint(0,len(self.config[cho]))
        cho2=np.random.randint(0,len(self.config[cho][cho1].T))
        cho3=np.random.randint(0,len(self.latvec[cho]))
        cho4=steplength*round(np.random.rand(),5)
        oldE=self.configE(1,cho)
        self.config[cho][cho1][cho3,cho2]=self.config[cho][cho1][cho3,cho2]+(cho0*cho4*self.latve[cho][cho3])/100
        if np.exp(-self.b*(self.configE(1,cho)-oldE))>=1:
            pass
        else:
            if np.random.random()<np.exp(-self.b*(self.configE(1,cho)-oldE)):
                pass
            else:
                self.config[cho1][cho3,cho2]=self.config[cho1][cho3,cho2]-(cho0*cho4*self.latvec[cho][cho3])/100
    
    
    
    def density(self,ind):
        return self.tot[ind]/(self.latvec[ind][0]**3)
    
    #隨機取樣與移動
    def rdf(self,lim):
        
        step=0.1
        for i in range(len(self.lab)):
            for k in range(len(self.lab)):
                data=np.zeros(10*lim)
                for z in range(len(self.config[i].T)):
                    dist=[]
                    for j in range(len(self.config[k].T)):
                        dist.append(self.dist(i,z,k,j))
                    for a in range(len(data)):
                        for b in dist:
                            if b>a*step and b<a*step+step:
                                data[a]=data[a]+1
                for c in range(len(data)):
                    data[c]=data[c]/((4*np.pi*(((c+1)*step))**3-(c*step)**3)/3)
                data=data/len(self.config[i].T)
                data[0]=0
                plt.figure()
                plt.plot(np.arange(0,10,0.1),data)
                
                
            
    
mod=model(lac,numbsit,rad,b,[],lab,sigma,epsilon)
mod.initconfig()
mod.picture(0)
eig=[]
eig.append(mod.configE(1))
for i in range(step):
    mod.ranchoiceandmove(10)
    eig.append(mod.configE(1))

plt.figure()
plt.plot(range(step+1),eig)
mod.picture(step)