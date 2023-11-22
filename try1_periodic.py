# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import matplotlib.pyplot as plt
import copy
#Minimum Imae Convention
# sigma  angstrom
# epsilon  kJ/mol
#Boltzmann constant  kJ/mol⋅K

Ar={'sigma':3.401,'epsilon':0.978638}
He={'sigma':2.556,'epsilon':0.08368}
bolt=8.314462618e-3  

numbsit=[[26],[26]]
lab=['Ar']
lac=np.array([[[100],[100],[100]],[[100],[100],[100]]])
rad=2.5
b=8.314462618e-3*87.3
step=10000
sigma=[3.401]
epsilon=[0.978638]


class model:
    def __init__(self,latvec,num,rad,b,config,lab,sigma,epsilon,cutoff):
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
        self.cutoff=cutoff
        
    def dist(self,y,i,j,k,z):
        return round(np.linalg.norm(self.config[y][i].T[j]-self.config[y][k].T[z]),5)
    
    def initconfig(self):
        self.config=[]
        con=[[] for i in self.latvec]
        for i in range(len(self.latvec)):
            self.config.append([])
            con[i]=np.random.rand(len(self.latvec[i]),self.tot[i])*self.latvec[i]-self.latvec[i]/2
        for k in range(len(self.latvec)):
            ac=0
            for i in self.num[k]:
                self.config[k].append(con[k][:,ac:ac+i])
                ac=ac+i
        while self.checkconfig(0)==False:
            con=np.random.rand(len(self.latvec[0]),self.tot[0])*self.latvec[0]-self.latvec[0]/2
            ac=0
            for i in self.num[0]:
                self.config[0].append(con[:,ac:ac+i])
                ac=ac+i
        while self.checkconfig(1)==False:
            con=np.random.rand(len(self.latvec[1]),self.tot[1])*self.latvec[1]-self.latvec[1]/2
            ac=0
            for i in self.num[1]:
                self.config[1].append(con[:,ac:ac+i])
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
        
    def configE(self,index,index2):
        if self.perio=='off':
            if index==1:
                utot=0
                for i in range(len(self.config[index2])):
                    for j in range(len(self.config[index2][i].T)):
                        for k in range(len(self.config[index2])):
                            for z in range(len(self.config[index2][k].T)):
                                if i==k and j==z:
                                    pass
                                elif self.dist(index2,i,j,k,z)>self.cutoff:
                                    pass
                                else:
                                    utot=utot+4*((epsilon[i]*epsilon[k])**0.5)*((((self.sigma[i]+self.sigma[k])/2)/self.dist(index2,i,j,k,z))**12-(((self.sigma[i]+self.sigma[k])/2)/self.dist(index2,i,j,k,z))**6)
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
                                            elif round(np.linalg.norm(self.config[index2][i].T[j]-self.config[index2][k].T[z]+np.array([a,b,c])*self.latvec[index2][0][0]),5)>self.cutoff:
                                                pass
                                            else:
                                                utot=utot+4*((epsilon[i]*epsilon[k])**0.5)*((((self.sigma[i]+self.sigma[k])/2)/round(np.linalg.norm(self.config[index2][i].T[j]-self.config[index2][k].T[z]+np.array([a,b,c])*self.latvec[index2][0][0]),5))**12-(((self.sigma[i]+self.sigma[k])/2)/round(np.linalg.norm(self.config[index2][i].T[j]-self.config[index2][k].T[z]+np.array([a,b,c])*self.latvec[index2][0][0]),5))**6)
                return utot/2 
    #計算config能量
    # 1 : Lennard-Jones
    def volchange(self,length):
        cho=np.random.choice([1,-1])*(length/100)*round(np.random.rand(),5)
        oldE=self.configE(1,0)+self.configE(1,1)
        mod.latvec[0]*(1+cho)
        mod.latvec[1]*(1-cho)
        for i in range(len(self.config[0])):
            self.config[0][i]=self.config[0][i]*(1+cho)
        for i in range(len(self.config[1])):
            self.config[1][i]=self.config[1][i]*(1-cho)
        if 1<round(((self.latvec[0][0][0])**(3*self.tot[0])*(self.latvec[1][0][0])**(3*self.tot[1])/((self.latvec[0][0][0]/(1+cho))**(3*self.tot[0])*(self.latvec[1][0][0]/(1+cho))**(3*self.tot[1])))*np.exp(-self.b*(self.configE(1,0)+(self.configE(1,1))-oldE)),5):
            pass
        else:
            if np.random.random()<round(((self.latvec[0][0][0])**(3*self.tot[0])*(self.latvec[1][0][0])**(3*self.tot[1])/((self.latvec[0][0][0]/(1+cho))**(3*self.tot[0])*(self.latvec[1][0][0]/(1+cho))**(3*self.tot[1])))*np.exp(-self.b*(self.configE(1,0)+(self.configE(1,1))-oldE)),5):
                pass
            else:
                print('VNO')
                mod.latvec[0]/(1+cho)
                mod.latvec[1]/(1-cho)
                for i in range(len(self.config[0])):
                    self.config[0][i]=self.config[0][i]/(1+cho)
                for i in range(len(self.config[1])):
                    self.config[1][i]=self.config[1][i]/(1-cho)
        
    
    def parexchange(self):
        cho=np.random.choice([0,1])
        cho1=np.random.randint(0,len(self.config[cho]))
        cho2=np.random.randint(0,len(self.config[cho][cho1].T))
        temp1=copy.deepcopy(self.config[0])
        temp2=copy.deepcopy(self.config[1])
        oldE=self.configE(1,0)+self.configE(1,1)
        temp3=self.config[cho][cho1].T
        temp4=np.delete(temp3,[cho2],0)
        temp5=np.append(self.config[abs(1-cho)][cho1].T,temp3[cho2].reshape(1,3),axis=0)
        self.config[cho][cho1]=temp4.T
        self.config[abs(1-cho)][cho1]=temp5.T
        if np.exp(-self.b*(self.configE(1,0)+self.configE(1,1)-oldE))*self.tot[cho]*self.latvec[abs(1-cho)][0][0]**3/((self.latvec[cho][0][0]**3)*(self.tot[abs(1-cho)]+1))>1:
            self.tot[cho]=self.tot[cho]-1
            self.tot[abs(1-cho)]=self.tot[abs(1-cho)]+1
            self.num[cho][cho1]=self.num[cho][cho1]-1
            self.num[abs(1-cho)][cho1]=self.num[abs(1-cho)][cho1]+1
        else:
            if np.random.random()<np.exp(-self.b*(self.configE(1,0)+self.configE(1,1)-oldE))*self.tot[cho]*self.latvec[abs(1-cho)][0][0]**3/((self.latvec[cho][0][0]**3)*(self.tot[abs(1-cho)]+1)):
                self.tot[cho]=self.tot[cho]-1
                self.tot[abs(1-cho)]=self.tot[abs(1-cho)]+1
                self.num[cho][cho1]=self.num[cho][cho1]-1
                self.num[abs(1-cho)][cho1]=self.num[abs(1-cho)][cho1]+1
            else:
                print('exNo')
                self.config[0]=temp1
                self.config[1]=temp2

    
    
    def ranchoiceandmove(self,steplength):
        cho=np.random.choice([0,1])
        cho0=np.random.choice([1,-1])
        cho1=np.random.randint(0,len(self.config[cho]))
        cho2=np.random.randint(0,len(self.config[cho][cho1].T))
        cho3=np.random.randint(0,len(self.latvec[cho]))
        cho4=steplength*round(np.random.rand(),5)
        oldE=self.configE(1,cho)
        self.config[cho][cho1][cho3,cho2]=self.config[cho][cho1][cho3,cho2]+(cho0*cho4*self.latvec[cho][cho3][0])/100
        if np.exp(-self.b*(self.configE(1,cho)-oldE))>=1:
            pass

        else:
            if np.random.random()<np.exp(-self.b*(self.configE(1,cho)-oldE)):
                pass
            else:
                print('movNO')
                self.config[cho1][cho3,cho2]=self.config[cho1][cho3,cho2]-(cho0*cho4*self.latvec[cho][cho3][0])/100

        
    
    
    
    def density(self,ind):
        return self.tot[ind]/(self.latvec[ind][0][0]**3)
    
    #隨機取樣與移動

                
                
            
    
mod=model(lac,numbsit,rad,b,[],lab,sigma,epsilon,10)
mod.initconfig()
for i in range(10):
    mod.ranchoiceandmove(10)
    mod.parexchange()
    mod.volchange(10)