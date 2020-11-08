# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 13:32:01 2020

@author: 19787146
"""
from dataclasses import dataclass
import numpy as np
from astropy.constants import G
from astropy import constants as const
from sgp4.api import Satrec
from sgp4.api import jday
import matplotlib.pyplot as plt


class Sattelite:
    def __init__(self,t,s):
        self.Sat_num =0
        self.Inclination =0
        self.RAAN=0
        self.Eccentricity=0
        self.Argument_of_perigee=0
        self.Mean_anom=0
        self.Mean_motion=0
        self.gm =0
        self.T =0
        self.a =0
        self.a =0
        self.A = 0
        self.E_anoma =0
        self.n = 0
        temp = ''
        i=0
        self.TrueAnom = 0
        self.Mnew =0
        self.u = const.GM_earth.value*7.464959999242148368082
        #self.u = const.GM_earth.value
        self.toll =0.0000001
        self.dell = 0
        self.alpha=0
        self.r = np.zeros(3)
        self.v = [0,0,0]
        
        count = 0
        count2 = 0
        while i<len(t):
            if(t[i]==' '):
                if t[i]==' ' and t[i+1]==' ':
                    i+=1
                count+=1
                if count==2:
                    self.Sat_num=int(temp)
                if count==3:
                    self.Inclination=float(temp)
                if count==4:
                    self.RAAN=float(temp)
                if count==5:
                    dot = '0.'
                    new= dot+temp.strip()
                    self.Eccentricity=float(new)
                if count==6:
                    self.Argument_of_perigee=float(temp)
                if count==7:
                    self.Mean_anom=float(temp)
                if count==8:
                    self.Mean_motion=float(temp)
                temp=''
            temp += t[i]
            i+=1
        self.Mean_motion=float(temp)
        i=0
        
        temp=''
        self.Inclination = self.Inclination*(np.pi/180)
        self.RAAN = self.RAAN*(np.pi/180) -(78*(np.pi/180))
        self.Argument_of_perigee = self.Argument_of_perigee*(np.pi/180)
        self.Mean_anom = self.Mean_anom*(np.pi/180)
        constt = 7.2722052166431
        self.n = self.Mean_motion*7.2722052166431*0.00001
        
        
    def printt(self):
#         Prints the second line elements of the Two Line Paramators.
        #print(' Sat_num:             ',self.Sat_num)
      #  print(' Inclination:         ',self.Inclination)
      #  print(' RAAN:                ',self.RAAN)
      #  print(' Eccentricity:        ',self.Eccentricity)
       # print(' Argument_of_perigee: ',self.Argument_of_perigee)
       # print(' Mean_anom:           ',self.Mean_anom)
      #  print(' Mean_motion:         ',self.Mean_motion)
      #  print(' n:                   ',self.n)
     #   print('\n')
        Sat_info = ' Sat_num:             ' + str(self.Sat_num) + '\n' + ' Inclination:         ' + str(self.Inclination) +'\n' + ' RAAN:                '+ str(self.RAAN)+'\n'+ ' Eccentricity:        '+str(self.Eccentricity)+'\n'+' Argument_of_perigee: '+str(self.Argument_of_perigee) +'\n'+' Mean_anom:           '+str(self.Mean_anom)  +'\n'+ ' Mean_motion:         '+ str(self.Mean_motion) +'\n\n'
        return Sat_info
    def calc(self,tt):
        self.T = 1/self.Mean_motion
       # 1 04793U 70106A   20300.41028892 -.00000052  00000-0 -50271-4 0  9994       
#         a is in km
       # date = datetime.datetime(2020,28,10)
        self.a = np.power((np.power(((self.T)/(2*np.pi)),2)*self.u),1/3)
        self.Mnew = self.Mean_anom + self.n*(tt)
        vtemp = self.E_anom()
        self.TrueAnom = vtemp
        ut = self.Argument_of_perigee+vtemp
        pp = self.a*(1-np.power(self.Eccentricity,2))
        r = pp/(1+self.Eccentricity*np.cos(self.TrueAnom))
        self.r[0] = r*(np.cos(self.RAAN)*np.cos(self.TrueAnom+self.Argument_of_perigee)-np.sin(self.RAAN)*np.sin(self.TrueAnom+self.Argument_of_perigee)*np.cos(self.Inclination))
        self.r[1] = r*(np.sin(self.RAAN)*np.cos(self.TrueAnom+self.Argument_of_perigee)+np.cos(self.RAAN)*np.sin(self.TrueAnom+self.Argument_of_perigee)*np.cos(self.Inclination))
        self.r[2] = r*(np.sin(self.Argument_of_perigee+self.TrueAnom)*np.sin(self.Inclination))
        
        return self.r

    
    def E_anom(self): 

#         print(temp)
        E0=self.Mnew
        En = self.Mnew
        if(En<np.pi):
            En+=En+(self.Eccentricity/2)
        else:
            En += En-(self.Eccentricity/2)
        
        for i in range(20):
            f = En - self.Eccentricity*np.sin(En)-self.Mnew
            df = 1-self.Eccentricity*np.cos(En)
            E0 = En
            En = En -(f/df)
#         vfinal = np.arccos((np.cos(En)-self.Eccentricity)/(1-self.Eccentricity*np.cos(En)))
        vfinal = 2*np.arctan(np.power((1+self.Eccentricity)/(1-self.Eccentricity),0.5)*np.tan(En/2))
#         print("En: ",En)
#         En = En -f/fd
#         f = En-self.Eccentricity*np.sin(En)
#         fd = 1-self.Eccentricity*np.cos(En)
#         print('\n ratio: ',f/fd)
        
#         True_A = 2*np.arctan(np.power((1+self.Eccentricity)/(1-self.Eccentricity),0.5)*np.tan(En/2))
#         f = np.arctan((self.Eccentricity*np.sin(True_A))/(1+self.Eccentricity*np.cos(True_A)))
#         print('\n True_A:',True_A)
# #         print('\n E_naught: ',En)
        return vfinal    
    #def gen_data_time():
        
    def gen_data(self,rev_per_day,points,offset,t):
        
        R_vec = np.zeros((3,points))
        temp = self.Mean_motion/(24*60*60)
        period = 1/(temp)
        #print(period)
        interval = period/points
#         print(interval)

        for i in range(points):
            #r = tmat.dot(temp.calc(tot))
            tmat = ECI2ECEF(t+i*interval)
            r = tmat.dot(self.calc(interval*(i)+offset))
            R_vec[0,i]=r[0]
            R_vec[1,i]=r[1]
            R_vec[2,i]=r[2]
  
        return R_vec
    
def ECI2ECEF(t):
    global ECEF_list
    #GMST
    we = 7.29211510*0.00001
    a = 0.7790572732640
   # b = 1.00273781191135448
    #thetag = we*t - 2*np.pi*(a+b*JD)
    thetag = we*t - 2*np.pi*(a)
    tmat=    np.zeros((3,3))
    tmat[0,0] = np.cos(thetag)
    tmat[0,1] = np.sin(thetag)
    tmat[1,0] =-np.sin(thetag)
    tmat[1,1] =np.cos(thetag)
    tmat[2,2] =1
    return tmat