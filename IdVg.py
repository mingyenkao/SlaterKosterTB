import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from math import sqrt
from math import exp
import xlrd
import csv
import os
import sys
def IdVg_Lg(Leff, Cg, Cs, Cd, lamda, Vg_offset,R):
	W=1E-6
	Vds=0.5
	Vinj= (84271.1660839 + 84209.8107604) / 2
	m0= 9.11E-31
	mx= 0.40638722 * m0
	my= 0.41535855 * m0
	m_eff= 2*sqrt(mx*my)
	e = 1.6E-19
	kB= 1.38064852E-23
	T = 300
	pi= 3.1415926
	hbar = 1.0545718E-34
	Vgs=0
	x_plot=[]
	y_plot=[]
	
	for i in range(1,200):
		err_I=1
		ID=1E-12
		while(abs(err_I) > 0.000001):
			U = 0.5*e
			del_n_0=0
			err_U=1
			Vgs=i *0.01 -1 -Vg_offset - R*ID
			Vds=0.5-2*R*ID
			counter=0
			while(abs(err_U) > 0.000001):
				del_n= m_eff*kB*T/(2*pi* (hbar**2)) * (exp(-U/(kB*T))) * (1+exp(-e*Vds/(kB*T)))
				#print (U)
				U_new=Cg/(Cg+Cs+Cd) * (-e*Vgs)+Cd/(Cg+Cs+Cd) * (-e)*Vds + e**2 * del_n /(Cg+Cs+Cd)
				err_U=(U_new-U)/U
				U=U+0.001*(U_new-U)
				del_n_0=del_n
				counter=counter+1

			
			ID_new=W*del_n_0*e*Vinj*(1/(1+(2*Leff/lamda)))*(1-exp(-e*Vds/(kB*T)))/(1+(exp(-e*Vds/(kB*T))/(1+(2*Leff/lamda))))
			err_I=(ID_new-ID)/ID
			ID=ID+(ID_new-ID)
			#print(ID)
			#print(ID_new)
			#print(err_I)
		x_plot.append(Vgs+Vg_offset)
		y_plot.append(ID)
		#print(counter)
		print(ID)
		#print(Vgs)
		
		#print(del_n_0)
	plt.semilogy(x_plot,y_plot)

def print_data_point(filename,style):	
	x_plot=[]
	y_plot=[]
	get_data=[]
	file = open(filename, 'r')
	csvCursor=csv.reader(file)
	for row in csvCursor:
		get_data.append(row)
	for i in range (len(get_data)):
		x_plot.append(get_data[i][0])
		y_plot.append(get_data[i][1])
	plt.semilogy(x_plot,y_plot,style)
IdVg_Lg(60E-9,0.01046404,0.009590497,0.002201557,1E-9,0.3,10000)
IdVg_Lg(30E-9,0.01046404,0.031204614,0.004400697,1E-9,0.6,10000)
IdVg_Lg(15E-9,0.01046404,0.024394252,0.007421493,1E-9,0.7,10000)
IdVg_Lg(60E-9,0.01046404,0.009590497,0.002201557,1E-9,0.3,0)
IdVg_Lg(30E-9,0.01046404,0.031204614,0.004400697,1E-9,0.6,0)
IdVg_Lg(15E-9,0.01046404,0.024394252,0.007421493,1E-9,0.7,0)
print_data_point('15nm_Lg_0.5Vds.csv','g^')
print_data_point('30nm_Lg_0.5Vds.csv','ro')
print_data_point('60nm_Lg_0.5Vds.csv','bs')
#0.029601808
plt.show()


