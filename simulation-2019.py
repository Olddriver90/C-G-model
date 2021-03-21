# This program was designed to simulate the process of evaporative isotope enrichment  on water surface.
# V------the volume of water mass
# E------evaporation loss from water surface
# δC:  isotopic value of canal water
# δE:  isotopic value of water vapour
# δA:  isotopic value of ambient moisture
# T:  temperature of water-vapour interface
# H:  relative humidity,range from 0.0 to 1.0
# α:  equilibrium fractionation factor,depending on temperature
# ε:  the total fractionation factor,ε=(ε*)+(εk), ε*=(α-1)/α; εk=Ck*(1-h)
T_C=16
T_K=T_C+273.15
h=0.65
δC_O=-7.39
δC_H=-50.63
δA_O=-14.90
δA_H=-103.74

import math
α_O=math.exp(-0.007685+6.7123/T_K-1666.4/T_K**2+350410/T_K**3)
ε1_O=(α_O-1)*1000
εK_O=14.2*(1-h) # Ck=14.2 for O
ε_O=εK_O+ε1_O/α_O

α_H=math.exp(1158.8*T_K**3/10**12-1620.1*T_K**2/10**9+794.84*T_K/10**6-161.04/10**3+2999200/T_K**3)
ε1_H=(α_H-1)*1000
εK_H=12.5*(1-h) # Ck=12.5 for H
ε_H=εK_H+ε1_H/α_H

f=1
dV=0.00001
x=1

result=[]
for i in range(1,1000):
    f=f*(1-dV)
    δE_O=(δC_O/α_O-h*δA_O-ε_O)/(1-h+εK_O*0.001)
    δE_H=(δC_H/α_H-h*δA_H-ε_H)/(1-h+εK_H*0.001)
    δC_O=(δC_O-x*dV*δE_O-(1-x)*dV*δC_O)/(1-dV)
    δC_H=(δC_H-x*dV*δE_H-(1-x)*dV*δC_H)/(1-dV)
    result.append([i,(1-f)*1,δE_O,δE_H,δC_O,δC_H])

import pandas as pd
names=['order','E/V','vapor_O','vapor_H','18O','2H']
data=pd.DataFrame(result, columns=names)
data.to_csv('C:/Users/11488/Desktop/manuscript/南水北调中线渠水蒸发研究/蒸发模型/Python/2019py/simulation-2019.csv', index=0)
