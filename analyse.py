# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 19:48:55 2014

@author: martin
"""

import csv
from numpy import *
import matplotlib.pyplot as plt
from uncertainties import  ufloat
from scipy.optimize import curve_fit

def load_scope(path, scale = True):
    data = open(path,'rb')
    time = []
    voltage = []
    reader = csv.reader(data, delimiter = ',')
    for line in reader:
        time.append(float(line[3]))
        voltage.append(float(line[4]))
    return array(time), array(voltage)
    
def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = False
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
 
    return array(maxtab), array(mintab)
    

def make_LaTeX_table(data,header, flip= 'false', onedim = 'false'):
    output = '\\begin{tabular}{'
    #Get dimensions
    if(onedim == 'true'):
        if(flip == 'false'):
        
            data = array([[i] for i in data])
        
        else:
            data = array([data])
        

    row_cnt, col_cnt = data.shape
    header_cnt = len(header)
    
    if(header_cnt == col_cnt and flip== 'false'):
        #Make Format
        output += '|'
        for i in range(col_cnt):
            output += 'c|'
        output += '}\n\\hline\n'+ header[0]
        for i in range (1,col_cnt):
            output += ' & ' + header[i]
        output += ' \\\\\n\\hline\n'
        for i in data:
            output += str(i[0])
            for j in range(1,col_cnt):
                output += ' & ' + str( i[j])
            output += '\\\\\n'
        output += '\\hline\n\\end{tabular}\n'
                            
        return output
    else:
        if(row_cnt == header_cnt):
            output += '|c|' + (col_cnt)*'c' + '|}\n\\hline\n'
            for i in range(row_cnt):
                output += header[i]
                for j in range(col_cnt):
                    output += ' & ' + str(data[i][j])
                output += '\\\\\n\\hline\n'
                
            output += '\\end{tabular}\n'
            return output
        else:
            return 'ERROR'

    
def lin_reg(x,y):
    N = len(x)
    sumx = x.sum()
    sumy = y.sum()
    sumxx = (x*x).sum()
    sumxy = (x*y).sum()
    m = (sumxy -  sumx*sumy/N)/(sumxx- sumx**2/N)
    b = sumy/N - m*sumx/N
    
    sy = sqrt(((y - m*x - b)**2).sum()/(N-1))
    m_err = sy *sqrt(N/(N*sumxx - sumx**2))
    b_err= m_err * sqrt(sumxx/N)
    return ufloat(m,m_err), ufloat(b,b_err)
    
    
###### ANFANG DER ANALYSE ##################################################
time, voltage = load_scope('./Messwerte/ALL0000/F0000CH1.CSV')

voltage = voltage -10.2 #offset aus CSV File
max_ , min_ = peakdet(voltage, 1, time )

max_ = max_.T

# Erste Grafik erzeugen
plt.close()
plt.plot(time*1000, voltage , 'k-', label = "Messwerte")
plt.plot(max_[0]*1e3,max_[1] ,'ro', markersize=5, label="Lokale Maxima")
plt.xlabel("t [ms]")
plt.ylabel("$U_C$ [V]",rotation='horizontal')
plt.xlim(0,time[-1]*1000)
plt.legend()
plt.savefig("./Abb/abb1.png")
# extrema zu einer liste hinzufügen und minima spigeln



regression_times = max_[0]
regression_voltages =max_[1]
print(make_LaTeX_table(array([regression_times*1e3,regression_voltages,[round(log(i),3) for i in regression_voltages]]).T,[r"$\frac{t}{ms}$",r"$\frac{U_C}{V}$",r"$ln(\frac{U_C}{V})$"],flip= 'false'))

plt.close()
plt.plot(1e3*regression_times,log(regression_voltages), 'x')

m,b = lin_reg(regression_times,log(regression_voltages))
print('m=%s'%str(m))
T_ex = -1/m
L= ufloat(10.11,0.03)*1e-3
C = ufloat(2.093,0.003)*1e-9
print("T_ex = 1/m = %s" % str(T_ex))
R_eff = 2*L/T_ex
print("R = 2L/T_ex = %s" % str(R_eff))

f_1res =  1/(2*pi)* (1/(L*C)- R_eff**2/(2*L**2))**0.5
q = 1/R_eff * (L/C)**0.5
breit = R_eff/ L 
f_1 = 1/(2*pi)*(-R_eff/(2*L)+(R_eff**2/(4*L**2)+1/(L*C))**0.5)
f_2 = 1/(2*pi)*(R_eff/(2*L)+(R_eff**2/(4*L**2)+1/(L*C))**0.5)
print('Guetfaktor: %s' % q)
print('Resonanzfrequez 1: %s' % f_1res)
print('Kurvenbreite : %s' % breit)
print('f_1: %s und f_2: %s'  % (f_1,f_2))
t = linspace(0,0.0005)

plt.plot(1e3*t,m.n * t +b.n )
plt.xlabel('t [ms]')
plt.ylabel(r'$ln(\frac{U_C}{V})$',rotation='horizontal')
plt.savefig('.\Abb\abb2.png')

print('#####################')
##### Berechnung des Aperiodischen Grenzfalls

R_apexp = ufloat(3220,20)

R_ap = (4*L/C)**0.5

print('Gemessener Wert: %s und berechneter Wert: %s' % (R_apexp,R_ap))


# Aufgabe c Resonanz

#Daten laden

data = loadtxt("./Messwerte/resonanz", unpack='true')

amplitude = data[2]/data[1]
freq = data[0]* 1000

max_res, min_res = peakdet(amplitude,1.0,freq)
bar = max_res[0][1]/sqrt(2)
x = linspace(0,90)
print('Resonanzfrequenz: %s und Guete des Schwingkreises: %s' % (str(max_res[0][0]),str(max_res[0][1])))
print()
plt.close()
plt.plot(freq*1e-3,amplitude, 'x')
plt.plot(x,bar+0*x)
plt.xlabel('f [kHz]')
plt.ylabel(r'$\frac{U_C}{U_0}$', rotation='horizontal')
plt.show()
plt.savefig('.\Abb\abb3.png')



### Phasenverschiebung  bestimmen 
file_id = ['03','04','05','06','07','08','09','10','11','12','13','14','15','16','17']
cnt_1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]
cnt_2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]
phase = []
freq = array([ 5000*i for i in range(2,17)])
def cos_fit(t, A, w , phi):
    return A* cos(t*w+phi)
    
for i in range(0,15):
    time , ch1 = load_scope('./Messwerte/ALL00%s/F00%sCH1.CSV' % (file_id[i],file_id[i]))
    time , ch2 = load_scope('./Messwerte/ALL00%s/F00%sCH2.CSV'% (file_id[i],file_id[i]))
    max_1,min_1 = peakdet(ch1,0.3,time)
    max_2, min_2 = peakdet(ch2,0.1, time)

    delta_t = max_1[cnt_1[0]][0] - max_2[cnt_2[0]][0]
    delta_phase = delta_t*freq[i]*360
    phase.append(delta_phase)
print phase


plt.close()
plt.plot(freq/1000,phase, 'x' , label= "Messwerte")
plt.xlabel('Frequenz [kHz] ')
plt.ylabel(r'$\Delta \varphi [\circ]$')
plt.legend()
plt.show()



    
