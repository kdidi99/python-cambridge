import pandas as pd
import numpy as np
import math
import os
import glob
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D     


def extract(path):
    #f = open(path, "r")
    r=0
    t=0
    e=0
    df = pd.DataFrame([], columns=["bondlength", "angle", "energy"])
    for infile in glob.glob( os.path.join(path, '*.out')):
        f = open(infile, "r")

        for line in f:
            #find r value
            if "Output" in line:
                mystr = str(line)
                char1 = 'r'
                char2 = 'h'
                r = float(mystr[mystr.find(char1)+1 : mystr.find(char2)-1])
                
                #find theta value
                char1 = 'e'
                char2 = '.out'
                t = float(mystr[mystr.find(char1)+3 : mystr.find(char2)])
                
            #find energy
            if "SCF Done:" in line:
                l = line.split()
                e = float(l[-5])
    
        df_tmp = pd.DataFrame([(r, t, e)], columns=["bondlength", "angle", "energy"])
        df = df.append(df_tmp)
        f.close()
    return df
#f= open("C:/Users/Kieran/Documents/Programming/Python_Cambridge/ex2/H2Ooutfiles\H2O.r0.70theta100.0.out", "r")


def plot(df):
    fig = plt.figure()
    ax = Axes3D(fig)
    trisurf = ax.plot_trisurf(df.bondlength, df.angle, df.energy, cmap=cm.jet, linewidth=0.2, edgecolor = 'grey')
    fig.colorbar(trisurf, ax = ax, shrink = 0.5, aspect = 5) 
    ax.set_title('Energy Surface')  
    # Adding labels 
    ax.set_xlabel('bondlength [A]', fontweight ='bold')  
    ax.set_ylabel('angle [degrees]', fontweight ='bold')  
    ax.set_zlabel('energy [kJ/mol]', fontweight ='bold') 
    plt.show()
    
def vibration(df):
 
    eqe2 = df['energy'].min()
    df_min = df[df.energy == df.energy.min()]
    eqr= float(df_min.bondlength)
    eqt= float(df_min.angle)
    eqe= float(df_min.energy)
    print("eqr:" +str(eqr))
    print("eqt:" +str(eqt))
    print("eqe:" +str(eqe))
    #subsetting dataframe to only fit near the equilibrium geometry
    #important to use .copy() to really create new df 
    rfitdata = df.loc[(df["angle"] == eqt) & (df["bondlength"] >= (eqr-0.2)) & (df["bondlength"] <= (eqr+0.2))].copy()
    print(rfitdata.head())
    print(rfitdata.info())
    #convert units for fit to m^2 (bondlength)and J(energy)
    rfitdata.bondlength = ((rfitdata.bondlength-eqr)*1e-10)**2
    rfitdata.energy*=4.35974e-18
    rfit = np.polyfit(rfitdata.bondlength, rfitdata.energy,2)
    #extract spring constant from linear fit
    kr = rfit[1]*2
    print("kr:" + str(kr))
    #calculate reduced mass and frequencies from data
    u1 = 2 * 1.66053886e-27
    u2 = 0.5 * 1.66053886e-27
    prefactor = 1/(2*(math.pi))
    print(prefactor)
    #calculate frequency of symmetric stretch from u1 and kr
    v1 = prefactor*math.sqrt(kr/u1)*3.3356e-11
    print("v1:" + str(v1))

    #calculate frequency of bend from u2 and kt
    tfitdata = df.loc[(df["bondlength"] == eqr) & (df["angle"] >= (eqt-5)) & (df["angle"] <= (eqt+5))].copy()
    print(tfitdata.head())
    print(tfitdata.info())
    #convert units for fit to rad(bond angle)and J(energy)
    tfitdata.angle = ((tfitdata.angle-eqt)*math.pi/180)**2
    tfitdata.energy*=4.35974e-18
    #fit data to get kt
    tfit = np.polyfit(tfitdata.angle, tfitdata.energy,2)
    print(tfit)
    kt= tfit[1]*2
    print("kt:" + str(kt))
    #calculate frequency for bend from u2 and kt
    v2 = prefactor*math.sqrt(kt/(u2*(eqr*1e-10)**2))*3.3356e-11
    print("v2:" + str(v2))
    
    
def main():
    path = input("Please enter the directory pathname for your input files:\n")
    df = extract(path)
    plot(df)
    vibration(df)


#C:\Users\Kieran\Documents\Programming\Python_Cambridge\ex2\H2Soutfiles
    
#C:/Users/Kieran/Documents/Programming/Python_Cambridge/ex2/H2Ooutfiles
