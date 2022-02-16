#!/usr/bin/env python
# coding: utf-8

# In[1]:


from multiprocessing import Pool
from ROOT import TTree, TBranch, TFile, TH1F, TCanvas
import numpy as np
from defs import*
P = Pool(8)


# In[2]:


def create_tree(data, events):
    det = np.int_([0])
    esc = np.int_([0])
    rel = np.int_([0])
    abso = np.int_([0])
    t = TTree('Ph_stat', 'PhotonStatistics_tree')
    t.Branch('detected',det, 'detected/I')
    t.Branch('escaped',esc, 'escaped/I')
    t.Branch('released',rel, 'released/I')
    t.Branch('absorbed',abso, 'absorbed/I')
    for i in range(events):
        det[0] = int(data[i][0])
        esc[0] = int(data[i][1])
        rel[0] = int(data[i][2])
        abso[0] = int(data[i][3])

        t.Fill()
    
    return t


# In[3]:


def create_hist(data, events, r_min, r_max, bins):
    histph    = TH1F( 'histph', 'Photon Detected vs Muon Events',bins, r_min, r_max )
    for i in range(events):
        histph.Fill(data[i][0])
    
    histph.Fit("landau")
    return histph


# In[4]:


def Write_file(hist, tree):
    f = TFile('PhotonStatisticsdata.root', 'RECREATE')
    hist.Write();
    tree.Write();
    f.Close()
    


# In[5]:


Events = 10000
data = P.map(Event, range(Events))  # Generates muon events Parallel


# In[6]:


c2 = TCanvas()
hist = create_hist(data, Events, r_min = 0, r_max = 1000, bins = 100)
hist.Draw()
c2.Draw()


# In[7]:


tree = create_tree(data,Events)
Write_file(hist, tree)

