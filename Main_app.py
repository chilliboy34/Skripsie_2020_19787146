# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 14:29:21 2020

@author: 19787146
"""
global Sat_num
Sat_num = 0
global sat_list
global launch_time
launch_time = 1597360978.725441
sat_list =[]



import ISS_spyder as sat

from tkinter import *
import time
from datetime import datetime
import plotly
import plotly.graph_objs as go
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from sgp4.api import Satrec
from sgp4.api import jday
import pandas as pd
import geopandas as gpd
import datetime, time
import pyproj as proj
import mayavi
import math
import matplotlib.patches as mpatches
import pytz
root = Tk()

def remove_sat():
    print("remove sat")
    temp = int(entry3.get())
    sat_list.pop(temp)
    
def load_sat():
    print("Loading sattelite")
    t = entry1.get()
    s = entry2.get()
    sat1 = sat.Sattelite(t,s)
    sat_list.append(sat1)
r_list = []
def track():
    global r_list
    
    r = np.zeros((3,200))
    points = 100
    SPG4_R_vec = np.zeros((3,points))
    per = (1/(15.49141851239476/(24*60*60)))
    interval = per/points
   
    ss1 = "1 25544U 98067A   20227.03014521  .00002523  00000-0  53627-4 0  9998"
    tt1 =  "2 00001  51.6461  62.1977 0001435  30.4470 104.2862 15.49164713240985"
    ss= '1 42982U 98067NE  20226.63852685  .00012705  00000-0  10205-3 0  9996'
    tt='2 42982  51.6351 349.1177 0003261 317.3782  42.6962 15.71399752160064'
    
    t_test ="2 00001  90  0 0000000000000001 30.4470 0 15.49164713240985"
    s_test = "1 25544U 98067A   20227.03014521  .00002523  00000-0  53627-4 0  9998"
    
    ISS = Satrec.twoline2rv(s_test,t_test)
# print('interval: ',interval)
# print(interval)

    for i in range(points):
        news = interval*i
        seconds = news % (24 * 3600) 
        hour = seconds // 3600
        seconds %= 3600
        minutes = seconds // 60
        seconds %= 60
        jd,fr = jday(2020,8,7,8+hour,0+minutes,0+seconds)
        e,r,v = ISS.sgp4(jd,fr)
        SPG4_R_vec[0,i]=r[0]
        SPG4_R_vec[1,i]=r[1]
        SPG4_R_vec[2,i]=r[2]
    x1 = SPG4_R_vec[0,:]
    y1 = SPG4_R_vec[1,:]
    z1 = SPG4_R_vec[2,:]
    c = ["blue","red","green","black","purple","brown","orange","cyan","grey","lime","teal","indigo","pink"]
    c_new = []
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #handles = []
    for i in range(len(sat_list)):
        temp = sat_list[i]
        r = temp.gen_data(temp.Mean_motion,100,0)
        r_list.append(r)
        x = r[0,:]
        y = r[1,:]
        z = r[2,:]
        #handles.append(str(sat_list[i].Sat_num))
        #c_new.append(c[i])
        ax.scatter(x,y,z,color=c[i],s=1,label=str(sat_list[i].Sat_num))
    ax.scatter(x1,y1,z1,color="indigo")
   # print(handles)
    ax.legend(title="Legend")
    s = 12742
    ax.scatter(0,0,0,s=s)
    
    #handles, labels1 = scatter.legend_elements(prop="sizes", alpha=0.6)
    #legend2 = ax.legend(handles, labels1, loc="upper right", title="Sizes")
    #ax.legend()
    ax.set_xlabel("X axis (km)")
    ax.set_ylabel("Y axis (km)")
    ax.set_zlabel("Z axis (km)")
    #ax.scatter(x1,y1,z1,color='red',s=1)
       # ax.legend()
    plt.show()
TOImat =0   
JD_global=0
timeload=0
def location():
    global JD_global 
    global TOImat
    global timeload
    throws=0
    print("Calculating location")    
    tempp = entry4.get()
    year, month, day, hour, minn, sec = map(int, tempp.split('-'))
   # year= int(tempp[0]+tempp[1]+tempp[2]+tempp[3])
   # month= int(tempp[4]+tempp[5])
  #  day=int(tempp[6]+tempp[7])
  #  hour=int(tempp[8]+tempp[9])
 #   minn=int(tempp[10]+tempp[11])
 #   sec=int(tempp[12]+tempp[13])
    
    JD_global,throws = jday(year,month,day,hour,minn,sec)
    d1 = datetime.datetime(2020,8,1,12,0,0)
    d2 = datetime.datetime(year, month, day,hour, minn, sec) 
    delta = d2-d1
    tot = delta.days*24*3600+delta.seconds
    timeload=tot
    print("Seconds_location",tot)
    tmat = ECI2ECEF(tot)
    TOImat = tmat
    c = ["pink","red","green","black","purple","brown","orange","cyan","grey","lime","teal","indigo","blue"]
    fig = plt.figure(figsize=(10,10)) #Adjusts the aspect ratio and enlarges the figure (text does not enlarge)
    ax = fig.gca(projection='3d')
    ax.set_zlim(-6900,6900)
    ax.set_xlim(-8000,8000)
    ax.set_ylim(-8000,8000)
    ECEF_list.clear() 
    ecef = proj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = proj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    lla_list =[]
    dicts = {'Sattelite':[],'Latitude':[],'Longitude':[]}
    dicts_orbit_map = {'Sattelite':[],'Latitude':[],'Longitude':[]}
    starts = np.zeros(3)
    starts[0] =3070.134092
    starts[1] =-4112.927314
    starts[2] =-3774.618733
    for i in range(len(sat_list)):
        lla_mat = np.zeros((3,len(sat_list)))
        temp = sat_list[i]
        r = tmat.dot(temp.calc(tot))
        #print("Seconds_single",tot)
        ECEF_coord = r.dot(tmat)
        
        print("Single location:",ECEF_coord)
        lon, lat, alt = proj.transform(ecef, lla, ECEF_coord[0]*1000, ECEF_coord[1]*1000, ECEF_coord[2]*1000, radians=True)
        lon = lon*100
        lat = lat*100
        alt = alt
        if (lon > 180):
            lon = 360-lon
        if (lon<-180):
            lon = -360-lon  
        if (lat > 90):
            lat = 180-lat
        if (lat<-90):
            lat = -180-lat     
        dicts['Sattelite'].append(str(i+1))
        dicts['Latitude'].append(lat)
        dicts['Longitude'].append(lon)
        
        lla_mat[0,i]=lon
        lla_mat[1,i]=lat
        lla_mat[2,i]=alt
        ECEF_list.append(ECEF_coord)
        xf = np.zeros(2)
        yf = np.zeros(2)
        zf = np.zeros(2)
        xf[0] =starts[0]
        xf[1]=r[0]
        yf[0]=starts[1]
        yf[1]=r[1]
        zf[0]=starts[2]
        zf[1]=r[2]
        ax.scatter(0,0,0,color="red",s=20)
       # ax.plot3D(xf,yf,zf,alpha=0.3)
       # ax.scatter(starts[0],starts[1],starts[2],s=100)
        ax.scatter(r[0],r[1],r[2],color=c[i],s=100,label=str(temp.Sat_num),alpha=1)
        #print("Vectore single value:",r)
    plt.show()
    ax.legend(title="Legend",loc="upper left")  
    
    if Checkk==1:
        for i in range(len(sat_list)):
            temp1 = sat_list[i]
            r1 = temp1.gen_data(temp1.Mean_motion,100,timeload,tot)
            x = r1[0,:]
            y = r1[1,:]
            z = r1[2,:]
            ax.scatter(x,y,z,color=c[i],s=5,label=str(sat_list[i].Sat_num))
        ax.scatter(x[0],y[0],z[0],color="red",s=15)
        ax.scatter(x[1],y[1],z[1],color="red",s=15)
        ax.scatter(x[2],y[2],z[2],color="red",s=15)
            
    for i in range(len(sat_list)): 
        lla_mat_map = np.zeros((3,len(sat_list)))
        r = temp.calc(tot)
        ECEF_coord = tmat.dot(r)
        lon, lat, alt = proj.transform(ecef, lla, ECEF_coord[0]*1000, ECEF_coord[1]*1000, ECEF_coord[2]*1000, radians=True)
        lon = lon*100
        lat = lat*100
        alt = alt
        if (lon > 180):
            lon = 360-lon
        if (lon<-180):
            lon = -360-lon
        
        if (lat > 90):
            lat = 180-lat
        if (lat<-90):
            lat = -180-lat
            
        dicts_orbit_map['Sattelite'].append(str(i+1))
        dicts_orbit_map['Latitude'].append(lat)
        dicts_orbit_map['Longitude'].append(lon)
                
    df = pd.DataFrame(dicts) 
    gdf = gpd.GeoDataFrame(
    df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude))
    world = gpd.read_file(
    gpd.datasets.get_path("naturalearth_lowres"))
    fig5 = plt.figure(figsize=(10,10))
    ax5 = fig5.add_subplot(111)     
    world.plot(ax=ax5, color="lightgray")
    ax5.set(xlabel="Longitude(Degrees)", ylabel="Latitude(Degrees)",title="WGS84 Datum (Degrees)")
    categories = df.Sattelite
    cat = str(np.linspace(0,len(sat_list),1))
    gdf.plot(ax=ax5, color=c, categorical=True,markersize=20)
    plt.show()   
         
    uee = np.linspace(0, np.pi, 40)
    vee = np.linspace(0, 2 * np.pi, 40)
    xee = np.outer(79.819*np.sin(uee), 79.819*np.sin(vee))
    yee = np.outer(79.819*np.sin(uee), 79.819*np.cos(vee))
    zee = np.outer(79.819*np.cos(uee), 79.819*np.ones_like(vee))
    ax.scatter(xee, yee, zee,s=2,alpha=0.4)
    ax.set_xlabel("X axis (km)")
    ax.set_ylabel("Y axis (km)")
    ax.set_zlabel("Z axis (km)")
    plt.show()   
    df = pd.DataFrame(dicts) 
    gdf = gpd.GeoDataFrame(
    df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude))
    
    
def Mission_Control():
    sim_time_start = datetime.datetime.now()
    start = entry5.get()
    end = entry6.get()
    launchd = datetime.datetime(2020,8,1,12,0,0)

    year, month, day, hour, minn, sec = map(int, start.split('-'))
    date1 = datetime.datetime(year, month, day,hour,minn,sec)
    delta1 = date1-launchd
    
    year1, month1, day1, hour1, minn1, sec1 = map(int, end.split('-'))
    date2 = datetime.datetime(year1, month1, day1,hour1,minn1,sec1)
    text.delete(1.0,END)
    sums = datetime.datetime(1970,1,1,0,0)
    delta45 = date1-sums
    launchsecs = (delta45.days*24*3600+delta45.seconds)
    
    tot1 = delta1.days*24*3600+delta1.seconds
    print("tot1_mission:",tot1)
    text.insert(END,'From:'+str(date1)+' To:'+str(date2)+'\nDuration:'+str(tot1)+'s')
    delta2 = date2-date1
    tot2 = delta2.days*24*3600+delta2.seconds
    print("sim duration:",tot2)
    angle2=0
    anglevec = np.zeros((tot2,2,len(sat_list)))
    anglevec1 = np.zeros((tot2,len(sat_list)))
    anglevec2 = np.zeros((tot2,len(sat_list)))
    south = np.zeros(3)
    start21 = datetime.datetime(2020,8,1,12,0,0)
    end21 = datetime.datetime(2020,8,1,12,0,0)
    beta = 30*np.pi/180
    thresh = 15*(np.pi/180) 
    Re = 6378
    alpha = 18*np.pi/180
    home_base_ECEF = {}
    Rearth = np.zeros(3)
    Rearth[0] =3070.134092
    Rearth[1] = -4112.927314
    Rearth[2] =-3774.618733
    south2satv = np.zeros(3)
    r22 = np.sqrt(np.power(Rearth[0],2)+np.power(Rearth[1],2)+np.power(Rearth[2],2))
    Rearth = Rearth/r22
    Rearthf = [3070.134092,-4112.927314,-3774.618733]
    l3 = np.zeros(3)
    x=0
    y=0
    z=0
    testss=0
    testss2 = 0
    checkpoint=0
    #need to calculate ground base ECEF coordinate
    rmain = np.zeros((tot2,3))
    text.delete(1.0,END)
    

    

    for i in range(len(sat_list)):
        temp = sat_list[i]
        for p in range(tot2):
            
            
            #datehold = datetime.datetime.fromtimestamp(tot1+p)
            #timeload = datehold-launchd
            #date1
            #timeload = date1.days*24*3600+timeload.seconds+1
            
            tmat = ECI2ECEF(tot1+p)
            rmain[p,:]=tmat.dot(temp.calc(tot1+p))
            hold = rmain[p,:]
            l3 = hold -Rearth
            r11 = np.sqrt(np.power(hold[0],2)+np.power(hold[1],2)+np.power(hold[2],2))
            r33 = np.sqrt(np.power(l3[0],2)+np.power(l3[1],2)+np.power(l3[2],2))
            hold = hold/r11
            l3 = l3/r33
            testss2 = l3.dot(Rearth)
            #angle =np.pi/2 - np.arccos(l3.dot(Rearth))
            angle =np.arccos(l3.dot(Rearth))
            anglevec[p,0,i] = angle
            if testss2>=np.cos(20*(np.pi/180)):
                anglevec[p,1,i] = 1   
                
                #print("Value of index:",p)
               # print("Seconds_mission:",tot1+p)
              #  print("vec_mission:",(temp.calc(tot1+p)).dot(tmat))
               # print("Vec: ",rmain[p,:],"Number: ",p)
              #  print("Seconds_mision",p+tot1)
            else:
                anglevec[p,1,i]=0
        text.insert(END,'\n'+'Sat_num: '+str(i)+'\n')     
            
            
        for op1 in range(tot2-1):
            if (anglevec[op1+1,1,i]==1) and (anglevec[op1,1,i]==0) and checkpoint ==0:
                d = datetime.timedelta(seconds=op1)
                new_timestamp = date1+d
                start21 = new_timestamp
                print("start")
                #print("doublecheck_vec:",(temp.calc(tot1+op1).dot(ECI2ECEF(tot1+op1))))
                #print("ofset_mision:",tot1)
               # print(rmain[op1,:])
               # print("Index 2:",op1)##datetime.datetime.utcfromtimestamp(launchsecs+op1)
                #print("total seconds",(op1+tot1))
               ## start21 = datetime.datetime.fromtimestamp(launchsecs+op1-(2*60*60))
                checkpoint=1
                
            if (anglevec[op1+1,1,i]==0) and (anglevec[op1,1,i]==1) and checkpoint==1:
                d = datetime.timedelta(seconds=op1)
                
                new_timestamp = date1+d
                end21 = new_timestamp ##datetime.datetime.utcfromtimestamp(launchsecs+op1)
               ## end21 = datetime.datetime.fromtimestamp(launchsecs+op1-(2*60*60))
                noday = (end21-start21).total_seconds()
                if (noday>0):
                    text.insert(END,' From: '+str(start21)+' To: '+str(end21)+'\n')
                    text.insert(END,' Duration: '+str(datetime.timedelta(seconds=(end21-start21).total_seconds()))+'\n')
                checkpoint=0
    sim_time_end = datetime.datetime.now()
    text.insert(END,'\n'+'Simulation runtime: '+str((sim_time_end-sim_time_start))+'\n')
    text.insert(END,'Simulation duration: '+str((date2-date1))+'\n')      
        
def listt_sat():
    text.delete(1.0,END)
    for i in range(len(sat_list)):
        text.insert(END,sat_list[i].printt())
        
def load():
    s=0
   # 2 04793 101.6432   6.2370 0031769  46.6969 324.9382 12.53996379282508 
   # sat_test = sat.Sattelite("2 20439  98.2033 269.8085 0010685 326.6893  33.3618 14.31689905815666",s)
 #   sat_list.append(sat_test)
    sat1 = sat.Sattelite("2 00001  90  10 00008 30.4470 0 15.3",s)
    sat_list.append(sat1)
    sat2 = sat.Sattelite("2 00002  90  10 00008 30.4470 45 15.3",s)
    sat_list.append(sat2)
    sat3 = sat.Sattelite("2 00003  90  10 00008 30.4470 90 15.3",s)
    sat_list.append(sat3)
    sat4 = sat.Sattelite("2 00004  90  10 00008 30.4470 135 15.3",s)
    sat_list.append(sat4)
    sat5 = sat.Sattelite("2 00005  90  10 00008 30.4470 180 15.3",s)
    sat_list.append(sat5)
    sat6 = sat.Sattelite("2 00006  90  10 00008 30.4470 225 15.3",s)
    sat_list.append(sat6)
    sat7 = sat.Sattelite("2 00007  90  10 00008 30.4470 270 15.3",s)
    sat_list.append(sat7)
    sat8 = sat.Sattelite("2 00008  90  10 00008 30.4470 315 15.3",s)
    sat_list.append(sat8)
Checkk=0

def change():
    global Checkk
    
    if Checkk==1:
        Checkk=0
    else:    
        Checkk=1
ECEF_list = []

LLA_list = []
def ECI2ECEF(t):
    global ECEF_list
    #GMST
    we = 7.29211510*0.00001
    a = 0.7790572732640
   # b = 1.00273781191135448
    #thetag = we*t - 2*np.pi*(a+b*JD)
    #thetag = we*(t+(32*60*60)) - 2*np.pi*(a)
    thetag = we*(t) - 2*np.pi*(a)
    tmat=    np.zeros((3,3))
    tmat[0,0] = np.cos(thetag)
    tmat[0,1] = np.sin(thetag)
    tmat[1,0] =-np.sin(thetag)
    tmat[1,1] =np.cos(thetag)
    tmat[2,2] =1
    return tmat
def ECEF2LLA():
    global LLA_list
    a = 6378137
    asq = np.power(a,2)
    f = 1/(298.2572223563)
    e = 0.08181919084
    es = np.power(e,2)
    b = 6356752.314  
    A = 42697.673
    B = 42841.312  
    LLA_list.clear()
    LLA=np.zeros(3)
    for i in range(len(ECEF_list)):
        r = ECEF_list[i]
        x = r[0]
        y = r[1]
        z = r[2]
        wsqrd = np.power(x,2)+np.power(y,2)
        l = es*(1/2)
        ls = np.power(l,2)
        m = wsqrd*(1/asq)
        n = np.power((1-es)*(z*(1/b)),2)
        p = (m+n-4*ls)*(1/6)
        G = m*n*ls
        H = 2*np.power(p,3)+G
        #Need to check that H<Hmin
        #print("test:",(m+n-4*ls))
        #print("H:",H)
        C = np.cbrt((H+G+2*math.sqrt(H*G)))*1/np.cbrt(2)
        
        
        i = -(2*ls+m+n)*(1/2)
        P = np.power(p,2)
        beta =i*(1/3)-C-(P*(1/C))
        betas = np.power(beta,2)
        k = ls*(ls-m-n)
        t = np.power((np.power(betas-k,1/2))-(beta+i)/2,1/2) -1*np.sign(m-n)*(np.power(np.abs((beta-i)/2),1/2))
        
        
        F = np.power(t,4)+2*i*np.power(t,2)+2*l*(m-n)*t+k
        df = 4*np.power(t,3)+4*i*t+2*l*(m-n)
        dellt = -F*(1/df)
        u = t+dellt+l
        v = t+dellt-l
        w = np.power(wsqrd,1/2)
        #print("\n",z,u,w,v,"\n")
        phi = np.arctan(z*u/w*v)
        dellw = w*(1-1/u)
        dellz = z*(1-(1-es)/v)
        h = np.sign(u-1)*np.power((np.power(dellw,2))+(np.power(dellz,2)),1/2)
        lamda = np.arctan(y/x)
        LLA[0] = phi
        LLA[1] = lamda
        LLA[2] = h
        LLA_list.append(LLA)
    
    

lab1 = Label(root,text = "t:")
lab2 = Label(root,text = "s:")
lab3 = Label(root,text ="RNumber:")
lab4 = Label(root,text = "Enter date:")
lab5 = Label(root,text = "Enter start date:")
lab6 = Label(root,text = "Enter end date:")
#lab5 = Label(root,text = "Load:")
entry1 =  Entry(root)
entry2 = Entry(root)
entry3 = Entry(root)
entry4 = Entry(root)
entry5 = Entry(root)
entry6 = Entry(root)


lab1.grid(row=0,sticky = E)
lab2.grid(row=1,sticky =E)
lab3.grid(row=2,sticky =E)
lab4.grid(row=3,sticky=E)
lab5.grid(row=4,sticky=E)
lab6.grid(row=5,sticky=E)

entry1.grid(row=0,column=1)
entry2.grid(row=1,column=1)
entry3.grid(row=2,column =1)
entry4.grid(row=3,column=1)
entry5.grid(row=4,column=1)
entry6.grid(row=5,column=1)

Button(root,text ="Load new sattelite", command = load_sat).grid(row=6,sticky=W,pady=4)
Button(root,text = "List loaded sattelites",command =listt_sat).grid(row=7,sticky=W,pady=4)
#Button(root,text = "Plot sattelite(s) Orbit(s)",command = track).grid(row=8,sticky=W,pady=4)
Button(root,text = "Remove Loaded sattelite",command=remove_sat).grid(row=9,sticky=W,pady=4)
Button(root,text = "Find location:",command=location).grid(row=10,sticky=W,pady=4)
Button(root,text = "Calculate connection schedule:",command=Mission_Control).grid(row=11,sticky=W,pady=4)
Button(root,text = "Load:",command=load).grid(row=12,sticky=W,pady=4)

text = Text(root,height=8, width=50)
text.grid(row=13, column=0, sticky="w",pady=4)
scroll_y = Scrollbar(root, orient="vertical", command=text.yview)
text.config(yscrollcommand=scroll_y.set)
scroll_y.grid(row=13, column=2, sticky="w",pady=4)
vart = IntVar()
c = Checkbutton(root, text = "Show Orbit",var=vart,onvalue=1,offvalue=0, command=change)
c.grid(column=0,row=12)
root.mainloop()
