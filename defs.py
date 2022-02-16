#!/usr/bin/env python
# coding: utf-8

from random import random
from math import sqrt,cos,sin,acos,asin,log,pi
from ROOT import gRandom
gRandom.SetSeed(0)

class muon:
    def __init__(self, x ,y, z, vx, vy, vz):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx 
        self.vy = vy
        self.vz = vz 
        self.hit_phs = 0
        self.absorbed = 0
        self.escaped = 0
        self.released = 0
    
    def __str__(self):
        return f"[{self.x}, {self.y}, {self.z}, {self.vx}, {self.vy}, {self.vz}]" 
            
    def propagate(self, sc):
        self.x = self.x  + self.vx
        self.y = self.y  + self.vy
        self.z = self.z  + self.vz
        
        if self.x > sc.xmin and self.x < sc.xmax         and self.y > sc.ymin and self.y < sc.ymax        and self.z > sc.zmin and self.z <= sc.zmax:
            return True
        
        else:
            return False
        
    def n_photons(self):
        n_ph = int(self.delE()/100)
        self.released += n_ph
        return n_ph
                
    
    def gen_photon(self):
        c = 1
        th = acos(1- 2*random())
        phi = 2*pi*random()
        vz = c*cos(th)
        vx = c*sin(th)*cos(phi)
        vy = c*sin(th)*sin(phi)
        return photon(self.x, self.y, self.z, vx, vy, vz)
    
    def delE(self):
        beta2 = 0.99934
        Z = 6
        A =12
        S = 1030*0.01
        eps = (0.1536/beta2)*(Z/A)*S
        Euler_c = 0.577
        lamda = gRandom.Landau()/eps
        dE = eps*(lamda + log((eps*(5.597*10**9)*beta2)/(1-beta2)*Z**2) + 1 - beta2 - Euler_c)
        #print(eps)

        return dE*10**3

            
        
        


# In[4]:


def gen_muon(scintillator):
    sc = scintillator
    
    dz = 0.01
    th = acos((1- .9375*random())**.25)
    phi = 2*pi*random()
    vz = -dz*cos(th)
    vx = dz*sin(th)*cos(phi)
    vy = dz*sin(th)*sin(phi)
    x = sc.xmin + (sc.xmax-sc.xmin)*random()
    y = sc.ymin + (sc.ymax-sc.ymin)*random()
    z = sc.zmax
    return muon(x, y, z, vx, vy, vz)


# In[5]:


class photon:
    def __init__(self, x , y, z, vx, vy, vz):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        
    def __str__(self):
        return f"[{self.x}, {self.y}, {self.z}, {self.vx}, {self.vy}, {self.vz}]" 
    
    def proptime(self, surface):
        sur = surface
        t = (sur.d - (sur.a*self.x + sur.b*self.y + sur.c*self.z))/(sur.a*self.vx + sur.b*self.vy + sur.c*self.vz + 10**-30)
        
        return t
        
    def propagate(self, scintillator):
        sc = scintillator
        min_t = 10**20
        for surf in sc.surfaces:
            t = self.proptime(surf)
            if t > 0 and t < min_t:
                min_t = t
        
        d = self.propdistance(min_t)
        #print(min_t, d)
        if self.absorption(d, sc.l_sc):
            return False
        else:
            self.x = self.x + self.vx*min_t
            self.y = self.y + self.vy*min_t
            self.z = self.z + self.vz*min_t
            
            if abs(self.x - sc.xmin) < 10**-8:
                self.x = sc.xmin
            elif abs(self.x - sc.xmax) < 10**-8:
                self.x = sc.xmax
            if abs(self.y - sc.ymin) < 10**-8:
                self.y = sc.ymin
            elif abs(self.y - sc.ymax) < 10**-8:
                self.y = sc.ymax
            if abs(self.z - sc.zmin) < 10**-8:
                self.z = sc.zmin
            elif abs(self.z - sc.zmax) < 10**-8:
                self.z = sc.zmax
            
            
            if not sc.in_scintillator(self):
                return False
        return True
        
        
    
        
    def propdistance(self, t): 
        x = self.x + self.vx*t
        y = self.y + self.vy*t
        z = self.z + self.vz*t
            
        return sqrt(x**2 + y**2 + z**2)
        
    def absorption(self,d, l_sc):
        if d > (-l_sc*log(1 - random())):
            return True
        else:
            return False
        
    def reflection(self, scintillator):
        sc = scintillator
        
        if abs(self.z - sc.zmin) < 10**(-8) or abs(self.z - sc.zmax) < 10**(-8):
            if sqrt((self.vx**2 + self.vy**2)/(self.vx**2 + self.vy**2 + self.vz**2 )) >  sc.eta:
                self.vz = - self.vz
                return True
            elif random() < sc.refl:
                self.vz = - self.vz
                return True
            else:
                return False
        
        elif abs(self.y - sc.ymin) < 10**(-8) or abs(self.y - sc.ymax) < 10**(-8):
            if sqrt(((self.vx**2 + self.vz**2)/self.vx**2 + self.vy**2 + self.vz**2 ) ) >  sc.eta:
                self.vy = - self.vy
                return True
            elif random() < sc.refl:
                self.vy = - self.vy
                return True
            else:
                return False         
        
        elif abs(self.x - sc.xmin) < 10**(-8) or abs(self.x - sc.xmax) < 10**(-8):
            if sqrt((self.vz**2 + self.vy**2)/(self.vx**2 + self.vy**2 + self.vz**2 )) >  sc.eta:
                self.vx = - self.vx
                return True
            elif random() < sc.refl:
                self.vx = - self.vx
                return True
            else:
                return False         
    
        
        
        


# In[6]:


class Detector:
    def __init__(self,xmin,xmax, ymin, ymax, zmin, zmax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax 
        self.zmin = zmin
        self.zmax = zmax 
    
    def check_hit(self,photon):
        ph = photon
        if ph.x>self.xmin and ph.x<self.xmax and ph.y>self.ymin and ph.y<self.ymax and ph.z>self.zmin and ph.z< self.zmax:
            return True
        else:
            return False


# In[7]:


class surface():
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        
    
            


# In[8]:


class Scintillator:
    def __init__(self,xmin,xmax, ymin, ymax, zmin, zmax, l_sc, eta, refl):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax 
        self.zmin = zmin
        self.zmax = zmax 
        self.l_sc = l_sc
        self.eta = eta
        self.refl = refl
        self.surfaces = [surface(0,0,1,zmin),surface(0,0,1,zmax),surface(1,0,0,xmin),surface(1,0,0,xmax),surface(0,1,0,ymin),surface(0,1,0,ymax)]

    def in_scintillator(self, photon):
        ph = photon
        if ph.x>=self.xmin and ph.x<=self.xmax and ph.y>=self.ymin and ph.y<=self.ymax and ph.z>=self.zmin and ph.z<= self.zmax:
            return True
        else:
            return False
        
        
def Event(i):
    dt = Detector( 12.5 - 0.00001, 12.5+ 0.00001 , -0.5,0.5, -0.5, 0.5)
    sc = Scintillator(-12.5, 12.5, -12.5, 12.5, -0.5, 0.5,100,1.59, 0.95 )
    
    #Events = 1000
    released = 0
    detected = 0
    absorbed = 0
    escaped =  0
    #for i in range(Events):
    mu = gen_muon(sc)
    while(mu.propagate(sc)):
        nph = mu.n_photons()
        for i in range(nph):
            ph = mu.gen_photon()
            for j in range(100):
                if ph.propagate(sc):
                    if dt.check_hit(ph):
                        mu.hit_phs += 1
                        break
                    elif not ph.reflection(sc):
                        mu.escaped += 1
                        break
                else:
                    mu.absorbed += 1
                    break
        detected = mu.hit_phs 
        escaped = mu.escaped
        released = mu.released
        absorbed = mu.absorbed
        
    #hist = create_hist(detected, Events)
    #tree = create_tree(detected,escaped,released,absorbed, Events)
    #Write_file(hist, tree)
    out = [detected,escaped,released,absorbed]
    return out





