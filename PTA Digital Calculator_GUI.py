# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 22:48:10 2021

@author: drewy
"""
#imports
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import time
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
import tkinter.font as tkfont
from tkinter import ttk
from matplotlib.widgets import SpanSelector
#Main window setup 
root = tk.Tk()
root.title('PTA Digital Calculator')
root.geometry('1400x900')
root.iconbitmap(r'C:\Users\drewy\OneDrive\TTU\Petroleum Engineering TTU\Reservoir Engineering Fundamentals\Python\TTU.ico')
#fonts
headm = tkfont.Font(family='Century Gothic',size=15,weight='bold')
headl2 = tkfont.Font(family='Arial',size=12)
mymenu = tk.Menu(root,tearoff=0)
filemenu = tk.Menu(mymenu,tearoff=0)
filemenu.add_command(label='New')
mymenu.add_cascade(label='File',menu=filemenu)
editmenu = tk.Menu(mymenu,tearoff=0)
editmenu.add_command(label='Undo')
mymenu.add_cascade(label='Edit',menu=editmenu)
tools = tk.Menu(mymenu,tearoff=0)
algo_menu = tk.Menu(tools,tearoff=0)
tools.add_cascade(label='Algorithm',menu=algo_menu)
algo = tk.IntVar()
algo_menu.add_radiobutton(label='Automatic',variable=algo,value=1)
algo_menu.add_radiobutton(label='Manual',variable=algo,value=2)
mymenu.add_cascade(label='Tools',menu=tools)
helpmenu = tk.Menu(root,tearoff=0)
mymenu.add_cascade(label='Help',menu=helpmenu)
root.config(menu=mymenu)
#Separator1
sep2 = ttk.Separator(root,orient='horizontal')
sep2.pack(fill='x',side='top')
#Resframe
resframe = tk.Frame(root)
resframe.pack(side='left',fill='both')
l1 = tk.Label(resframe, text='Test Data',font=headm)
l1.grid(row=0,column=1,padx=5,pady=5)
l2=tk.Label(resframe, text='Filename:',font=headl2)
l2.grid(row=1,column=0,padx=5,pady=5)
filevar = tk.StringVar()
tp = None
def open_file():
    filename = filedialog.askopenfilename(defaultextension='.csv',filetypes=\
    [('All Files','*.*'),('Comma Separated Values','*.*')])
    if filename == '':
        filename = None 
    else:
        filevar.set(filename)
e1 = tk.Entry(resframe,textvariable=filevar).grid(row=1,column=1,padx=5,pady=15,ipadx=5,ipady=5)
b1 = tk.Button(resframe, text='Browse',relief='raised',command=open_file)
b1.grid(row=1,column=2,padx=7,pady=7)
l3 = tk.Label(resframe,text='Reservoir Data',font=headm).grid(row=2,column=1,padx=5,pady=25)
#input reservoir variables
params = ['Viscosity','Formation Volume Factor','Height','Porosity',\
          'Compressibility Factor','Rate','Wellbore Radius']
u = tk.StringVar(); FVF = tk.StringVar(); h = tk.StringVar(); phi=tk.StringVar()
ct = tk.StringVar(); q=tk.StringVar(); rw=tk.StringVar()
invar = [u,FVF,h,phi,ct,q,rw]
k=3
for i in range(len(params)):
    tk.Label(resframe,text=params[i],font=headl2).grid(row=k,column=i%2,padx=5,pady=5,ipadx=5,ipady=5)
    tk.Entry(resframe,textvariable=invar[i]).grid(row=k+1,column=i%2,padx=5,pady=5,ipadx=5,ipady=5)
    if (i%2>0):
        k+=2
var_r1 = tk.IntVar()
r1 = tk.Radiobutton(resframe,text='Build up',font=headl2,variable=var_r1,value=1)
r1.grid(row=11,column=0,padx=5,pady=5,ipadx=5,ipady=5)
r2 = tk.Radiobutton(resframe,text='Drawdown',font=headl2,variable=var_r1,value=2)
r2.grid(row=12,column=0,padx=5,pady=5,ipadx=5,ipady=5)
k = tk.StringVar(); s = tk.StringVar(); PD = tk.StringVar(); FE = tk.StringVar();path_im = tk.StringVar();IP=tk.StringVar()
#Separator2
sep2 = ttk.Separator(root,orient='vertical')
sep2.pack(fill='y',side='left',expand=True)
#Graph Frame
frame2 = tk.Frame(root)
frame2.pack(side='left',fill='both',expand='True')
fig = plt.Figure()
fig.suptitle('Plotting Area',fontsize=16)
c1 = FigureCanvasTkAgg(fig,master=frame2)
c1.draw()
c1.get_tk_widget().pack(side='top',fill='both',expand='true')
tool = NavigationToolbar2Tk(c1,frame2,pack_toolbar=False)
tool.update()
tool.pack(side='bottom',fill='both')
out_var = [k,s,PD,FE,IP]
def calculate():
    global path_im
    path_im = tk.StringVar()
    path = os.path.dirname(filevar.get())
    os.chdir(path)
    fil = os.path.basename(filevar.get())
    file = pd.read_csv(fil)
    #Auto Drawdown
    if (var_r1.get() == 2) & (algo.get()==1):
        Pi = file.loc[0,'P']
        file = file.drop(file[file['t']==0].index)
        file = file.reset_index()
        dp = np.log10(file['P']).diff()
        dt = np.log10(file['t']).diff()
        grad = dp/dt
        print(grad)
        diff_grad = grad.diff()
        print(diff_grad)
        val = np.where(abs(diff_grad) < 0.005)
        val = np.array(val)
        mini = val.argmin(); maxi = val.argmax()
        while val[0][mini+2]-val[0][mini] !=2:
            mini +=1
        while abs(val[0][maxi-2]-val[0][maxi]) !=2:
            maxi -=1
        tmin = val[0][mini]-2; tmax=val[0][maxi]
        mtr_x = file.loc[tmin:tmax+1,'t']; mtr_y=file.loc[tmin:tmax+1,'P']
        coeff = np.polyfit(np.log10(mtr_x),mtr_y,1)
        min_p = file['P'].min(); min_t = np.log10(file['t'].min())
        dify = (file['P'].max()-min_p)/3; difx = (np.log10(file['t'].max())-min_t)/2
        ax = fig.add_subplot(1,1,1)
        ax.plot(file['t'],file['P'],'bo')
        ax.plot(mtr_x,coeff[0]*np.log10(mtr_x)+coeff[1],'r--')
        ax.set_xscale('log')
        ax.set_xlabel('Time')
        ax.set_ylabel('Pressure')
        ax.grid(True)
        ax.set_title('Drawdown Test')
        ax.text(min_t+difx,min_p+dify,'%.2fx+%.2f'%(coeff[0],coeff[1]),color='red')
        h2=float(h.get()); rw2=float(rw.get());u2=float(u.get());ct2=float(ct.get())
        q2=float(q.get());phi2=float(phi.get());FVF2=float(FVF.get()) 
        print(u.get(),FVF.get(),h.get(),phi.get(),ct.get(),rw.get(),q.get())
        k2 =-162.6*FVF2*q2*u2/(coeff[0]*h2)
        k.set(str(round(k2,4))+' '+'mD')
        P1 = coeff[0]*np.log10(1)+coeff[1]
        s2 = 1.151*(((P1-Pi)/coeff[0])+np.log10((phi2*u2*ct2*rw2**2)/k2)+3.23)
        s.set(str(round(s2,4)))
        pd2 = 0.869*-coeff[0]*s2
        PD.set(str(round(pd2,2))+' '+'psi')
        fe2 = ((Pi-P1)-pd2)/(Pi-P1)
        FE.set(str(round(fe2,2)))
        IP.set(Pi)
        print(k.get(),s.get())
    #Manual Drawdown
    elif (var_r1.get() == 2) & (algo.get()==2):
        Pi = file.loc[0,'P']
        file = file.drop(file[file['t']==0].index)
        file = file.reset_index()
        ax = fig.add_subplot(1,1,1)
        ax.plot(file['t'],file['P'],'bo')
        ax.set_xscale('log')
        ax.set_xlabel('Time')
        ax.set_ylabel('Pressure')
        ax.grid(True)
        ax.set_title('Drawdown Test')
        def onselect(tmin,tmax):
            indmin,indmax = np.searchsorted(file['t'],(tmin,tmax))
            indmax = min(len(file['t'])-1,indmax)
            global region_x; global region_y
            region_x = file.loc[indmin:indmax+1,'t']
            region_y =  file.loc[indmin:indmax+1,'P']
        span = SpanSelector(ax, onselect,'horizontal',useblit=True,rectprops=dict(alpha=0.5,facecolor='red'))
        coeff = np.polyfit(np.log10(region_x ),region_y,1)
        ax.plot(region_x,coeff[0]*np.log10(region_x)+coeff[1],'r--')
        min_p = file['P'].min(); min_t = np.log10(file['t'].min())
        dify = (file['P'].max()-min_p)/3; difx = (np.log10(file['t'].max())-min_t)/2
        ax.text(min_t+difx,min_p+dify,'%.2fx+%.2f'%(coeff[0],coeff[1]),color='red')
        h2=float(h.get()); rw2=float(rw.get());u2=float(u.get());ct2=float(ct.get())
        q2=float(q.get());phi2=float(phi.get());FVF2=float(FVF.get())
        k2 =-162.6*FVF2*q2*u2/(coeff[0]*h2)
        k.set(str(round(k2,4))+' '+'mD')
        P1 = coeff[0]*np.log10(1)+coeff[1]
        s2 = 1.151*(((P1-Pi)/coeff[0])+np.log10((phi2*u2*ct2*rw2**2)/k2)+3.23)
        s.set(str(round(s2,4)))
        pd2 = 0.869*-coeff[0]*s2
        PD.set(str(round(pd2,2))+' '+'psi')
        fe2 = ((Pi-P1)-pd2)/(Pi-P1)
        FE.set(str(round(fe2,2)))
        IP.set(Pi)
    #Auto buildup
    elif (var_r1.get() == 1) & (algo.get()==1):
        tp1 = tk.StringVar()
        t2 = tk.Toplevel(root)
        t2.title('Enter Production Time')
        t2.geometry('250x50')
        t2.transient(root)
        e2 = tk.Entry(t2,textvariable=tp1)
        e2.pack()
        def close():
            global tp
            tp = float(tp1.get())
            t2.destroy()
        bt = tk.Button(t2, text='Enter',command=close)
        bt.pack()
        root.wait_window(t2)
        Pwf0 = file.loc[0,'P']
        file = file.drop(file[file['t']==0].index)
        file = file.reset_index()
        file['x'] = file['t'].apply(lambda x: (tp+x)/x)
        dp = np.log10(file['P']).diff()
        dx = np.log10(file['x']).diff()
        grad = dp/dx
        diff_grad = grad.diff()
        val = np.where(abs(diff_grad) < 0.005)
        val = np.array(val)
        mini = val.argmin(); maxi = val.argmax()
        while val[0][mini+2]-val[0][mini] !=2:
            mini +=1
        while abs(val[0][maxi-2]-val[0][maxi]) !=2:
            maxi -=1
        xmin = val[0][mini]-2; xmax=val[0][maxi]
        mtr_x = file.loc[xmin:xmax+1,'x']; mtr_y=file.loc[xmin:xmax+1,'P']
        coeff = np.polyfit(np.log10(mtr_x),mtr_y,1)
        Pi = coeff[1]; IP.set(str(round(coeff[1],2))+' '+'psi')
        min_p = file['P'].min(); min_x = np.log10(file['x'].min())
        dify = (file['P'].max()-min_p)/3; difx = (np.log10(file['x'].max())-min_x)/2
        ax = fig.add_subplot(1,1,1)
        ax.plot(file['x'],file['P'],'bo')
        ax.plot(mtr_x,coeff[0]*np.log10(mtr_x)+coeff[1],'r--')
        ax.set_xscale('log')
        ax.set_xlabel('$t_{p}$+t/t')
        ax.set_ylabel('Pressure')
        ax.grid(True)
        ax.set_title('Build Up Test')
        ax.text(min_x+difx,min_p+dify,'%.2fx+%.2f'%(coeff[0],coeff[1]),color='red')
        h2=float(h.get()); rw2=float(rw.get());u2=float(u.get());ct2=float(ct.get())
        q2=float(q.get());phi2=float(phi.get());FVF2=float(FVF.get()) 
        print(u.get(),FVF.get(),h.get(),phi.get(),ct.get(),rw.get(),q.get())
        k2 =-162.6*FVF2*q2*u2/(coeff[0]*h2)
        k.set(str(round(k2,4))+' '+'mD')
        print(k.get())
        P1 = coeff[0]*np.log10((tp+1)/1)+coeff[1]
        print(P1)
        s2 = 1.151*(((P1-Pwf0)/-coeff[0])+np.log10((phi2*u2*ct2*rw2**2)/k2)+3.23)
        s.set(str(round(s2,4)))
        print(s.get())
        pd2 = 0.869*-coeff[0]*s2
        PD.set(str(round(pd2,2))+' '+'psi')
        fe2 = ((Pi-Pwf0)-pd2)/(Pi-Pwf0)
        FE.set(str(round(fe2,2)))
        print(k.get(),s.get())
    #Manual Buildup
    elif (var_r1.get() == 1) & (algo.get()==2):
        tp1 = tk.StringVar()
        t2 = tk.Toplevel(root)
        t2.title('Enter Production Time')
        t2.geometry('250x50')
        t2.transient(root)
        e2 = tk.Entry(t2,textvariable=tp1)
        e2.pack()
        def close():
            global tp
            tp = float(tp1.get())
            t2.destroy()
        bt = tk.Button(t2, text='Enter',command=close)
        bt.pack()
        root.wait_window(t2)
        Pwf0 = file.loc[0,'P']
        file = file.drop(file[file['t']==0].index)
        file = file.reset_index()
        file['x'] = file['t'].apply(lambda x: (tp+x)/x)
        ax = fig.add_subplot(1,1,1)
        ax.plot(file['x'],file['P'],'bo')
        ax.set_xscale('log')
        ax.set_xlabel('$t_{p}$+t/t')
        ax.set_ylabel('Pressure')
        ax.grid(True)
        ax.set_title('Build Up Test')
        region_x = None; region_y = None
        def onselect(xmin,xmax):
            indmin,indmax = np.searchsorted(file['x'],(xmin,xmax))
            indmax = min(len(file['x'])-1,indmax)
            global region_x; global region_y
            region_x = file.loc[indmin:indmax+1,'x']
            region_y =  file.loc[indmin:indmax+1,'x']
        span = SpanSelector(ax,onselect,'horizontal',useblit=True,rectprops=dict(alpha=0.5,facecolor='red'))
        if region_x != None:
            coeff = np.polyfit(np.log10(region_x),region_y,1)
            ax.plot(region_x,coeff[0]*np.log10(region_x)+coeff[1],'r--')
            min_p = file['P'].min(); min_t = np.log10(file['x'].min())
            dify = (file['P'].max()-min_p)/3; difx = (np.log10(file['x'].max())-min_t)/2
            ax.text(min_t+difx,min_p+dify,'%.2fx+%.2f'%(coeff[0],coeff[1]),color='red')
            h2=float(h.get()); rw2=float(rw.get());u2=float(u.get());ct2=float(ct.get())
            q2=float(q.get());phi2=float(phi.get());FVF2=float(FVF.get())
            k2 =-162.6*FVF2*q2*u2/(coeff[0]*h2)
            k.set(str(round(k2,4))+' '+'mD')
            P1 = coeff[0]*np.log10(1)+coeff[1]
            s2 = 1.151*(((P1-Pi)/coeff[0])+np.log10((phi2*u2*ct2*rw2**2)/k2)+3.23)
            s.set(str(round(s2,4)))
            pd2 = 0.869*-coeff[0]*s2
            PD.set(str(round(pd2,2))+' '+'psi')
            fe2 = ((Pi-P1)-pd2)/(Pi-P1)
            FE.set(str(round(fe2,2)))
            IP.set(Pi)    
#Separator3
sep3 = ttk.Separator(root,orient='vertical')
sep3.pack(fill='y',side='left',expand=True)
#DataFrame
resuframe = tk.Frame(root)
resuframe.pack(side='right',expand=True, fill='both')
l5 = tk.Label(resuframe,text='Results',font=headm).grid(row=0,column=0,padx=5,pady=25)
output = ['Permeability:','Skin:','Pressure Drop:','Flow Efficiency:','Initial Pressure:']
for i in range(len(output)):
    tk.Label(resuframe,text=output[i],font=headl2).grid(row=i+1,column=0,padx=10,pady=10)
for i in range(len(out_var)):
    tk.Label(resuframe,textvariable=out_var[i],font=headl2).grid(row=i+1,column=1,padx=5,pady=5)
def clear():
    for i in invar:
        i.set('')
    for j in range(len(out_var)):
        out_var[j].set('')
    filevar.set('')
    fig.clf()
    var_r1.set(0)
l6 = tk.Label(resuframe, text='Created by: Andrew Addo-Yobo').grid(row=9,column=0,padx=5,pady=530)
b2 = tk.Button(resframe,text='Calculate',command=calculate).grid(row=13,column=1,padx=5,pady=5)
b3 = tk.Button(resframe,text='Clear',command=clear).grid(row=13,column=2,padx=5,pady=5)
root.mainloop()
    
    