# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 12:55:52 2018

@author: Martin
"""

#importar librerias 
import matplotlib.pyplot as plt
import numpy as np
import random
import pylab
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.datasets import load_svmlight_file
from sklearn.datasets.samples_generator import make_regression 
from scipy import stats

#def gradiente():

# asignacion de valores a las variables
x_tem=load_svmlight_file("C:\\Users\\Martin\\Downloads\\ex2x.dat")
y_tem=load_svmlight_file("C:\\Users\\Martin\\Downloads\\ex2y.dat")  

#asignacion de 
x=x_tem[1]
y=y_tem[1]

plt.scatter(x, y, color='black')
plt.xlabel('Edad en años')
plt.ylabel('Altura en metros')


m=len(y)#número de muestras 
#x=np.transpose(np.vstack((np.ones(m),x)))
alpha=0.01 #float(input("Alpha=  "))
#thetas
theta0=np.linspace(-3,3,1)#np.random.random()#np.linspace(-3,3,100))
theta1=np.linspace(-1,1,1)#random.random()#linspace(-1,1,100)
listo= False 
while not listo:
#error
        J = sum([(theta0 + theta1*x[i] - y[i])**2 for i in range(m)])

#gradiente 
        grad0= 1/m * sum([(theta0 + theta1*x[i] - y[i]) for i in range(m)])
        grad1= 1/m * sum([(theta0 + theta1*x[i] - y[i])*x[i] for i in range(m)])

#actualizacion 
        t0 = theta0 - alpha * grad0
        t1 = theta1 - alpha * grad1

        theta0=t0 
        theta1=t1

        e = sum( [ (theta0 + theta1*x[i] - y[i])**2 for i in range(m)] ) 
        if abs(J-e) <= 0.00001:
            listo = True 
            print("Ya sta")
            
for i in range (0,50):
    y_p= theta0 + theta1*x
    pylab.plot(x,y_p,'-')
    
   # return theta0,theta1
#if __name__ == '__main__':
#theta0, theta1 = gradiente(x,y)
#slope, intercept, r_value, p_value, slope_std_error = stats.linregress(x[:,1], y)
#for i in range(0,50):
#        y_predict = theta0 + theta1
#pylab.plot(x,y_predict,'-')
#    
