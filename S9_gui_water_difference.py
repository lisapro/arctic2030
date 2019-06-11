'''
Created on 28. jun. 2017

@author: ELP
'''

import os,sys,datetime 
from PyQt5 import QtWidgets,QtGui, QtCore
from PyQt5.QtWidgets import QTableWidget,QTableWidgetItem
from netCDF4 import Dataset,num2date,date2num,date2index
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.dates as mdates
import tkinter as tk 
from tkinter.filedialog import askopenfilename,askdirectory  
import matplotlib.ticker as ticker
import seaborn as sns
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas)
from matplotlib.backends.backend_qt5agg import (
    NavigationToolbar2QT as NavigationToolbar)
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from dateutil.relativedelta import relativedelta
sns.set() 
root = tk.Tk()
root.withdraw()

class Window(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window) 
        self.figure = plt.figure(figsize=(9.3 ,3), dpi=150,
                        facecolor='None',edgecolor='None')        
        self.canvas = FigureCanvas(self.figure)  
               
        self.directory =  askdirectory()  #self.load_work_directory() 
        self.water_fname = os.path.join(self.directory,'water.nc')
        self.water_base_fname = r'C:\Users\elp\OneDrive - NIVA\BROM_linux_output\Baseline_B\water.nc'
            
        #self.water_base_fname = os.path.join(self.directory,'baseline_water.nc')

        self.fh_water_base = Dataset(self.water_base_fname)
        self.fh_water =  Dataset(self.water_fname)   
        try:
            self.time = self.fh_water.variables['time'][:]
            units = self.fh_water.variables['time'].units     
        except KeyError: 
            self.time = self.fh_water.variables['ocean_time'][:]
            units = self.fh_water.variables['ocean_time'].units                    
        self.format_time = num2date(self.time,units = units,
                                    calendar= 'standard')   
            
        first_year = self.format_time[0].year
        last_year = self.format_time[-1].year
        
        self.fontsize = 12
                
        self.names_vars = [] 
        for names,vars in self.fh_water.variables.items():
            self.names_vars.append(names)
                    
        self.names_vars =  sorted(self.names_vars, key=lambda s: s.lower())  
               
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                        QtWidgets.QSizePolicy.Expanding)
        

        self.toolbar = NavigationToolbar(self.canvas, self)
        self.button = QtWidgets.QPushButton('Plot')
        self.save_button = QtWidgets.QPushButton('Plot and Save_pdf')
        #self.label_choose_var = QtWidgets.QLabel('Choose variable:')        
        self.change_title = QtWidgets.QLineEdit()
        self.symbols = QtWidgets.QLineEdit('Symbols: # \ $ / & % {} ^ *')
        self.checkbox_title = QtWidgets.QCheckBox('Change title ')
        
        self.label_start_year = QtWidgets.QLabel('Start year:') 
        self.combobox_start_year = QtWidgets.QSpinBox()    
        self.combobox_start_year.setRange(first_year, last_year-1)                 
        self.label_stop_year = QtWidgets.QLabel('Stop year:')   
        self.combobox_stop_year = QtWidgets.QSpinBox() 
        self.combobox_stop_year.setRange(first_year+1, last_year+1)                        
        self.qlist_widget = QtWidgets.QListWidget()        
        self.qlist_widget.addItems(self.names_vars)
      
        layout = QtWidgets.QGridLayout()
        
        # fist line 
        layout.addWidget(self.toolbar,0,0,1,4)  
        layout.addWidget(self.checkbox_title,1,0,1,1)         
        layout.addWidget(self.change_title,1,1,1,1)
        layout.addWidget(self.symbols,1,2,1,1)                         
        layout.addWidget(self.button,1,3,1,1)  
               
        #layout.addWidget(self.save_button,1,4,1,1)
        layout.addWidget(self.label_start_year,0,4,1,1)
        layout.addWidget(self.combobox_start_year,1,4,1,1) 
               
        layout.addWidget(self.label_stop_year,0,5,1,1)
        layout.addWidget(self.combobox_stop_year,1,5,1,1)                 

        # third line             
        layout.addWidget(self.qlist_widget,2,0,1,1)      
        layout.addWidget(self.canvas,2,1,1,8)    
                                    
        self.setLayout(layout)  
              
        self.button.released.connect(self.call_show_3fig)  
        #self.save_button.released.connect(self.call_save_3fig) 
    
    def call_show_3fig(self):
        self.action = 'showfigure'
        self.plot_3fig()
     
    def call_save_3fig(self):
        self.action = 'savepdf'
        self.plot_3fig()
        
    def read_var(self):

        self.name =  str(self.qlist_widget.currentItem().text())
        self.long_name = str(self.fh_water.variables[str(self.name)].long_name)
        self.setWindowTitle(str(self.directory[-20:] + self.long_name)) 

        var_water = np.array(
            self.fh_water.variables[self.name][:]).T  
        var_water_base  = np.array(
            self.fh_water_base.variables[self.name][:]).T     
        var_water_dif = var_water-var_water_base              
        data_units = self.fh_water.variables[self.name].units
        if len(self.change_title.text()) < 1: 
            self.change_title.setText(self.name+' '+ data_units)                    
        return var_water_dif,data_units

    def save_to_dir(self,dir_name):
        script_dir = os.path.abspath(os.path.dirname(__file__))
        dir_to_save = os.path.abspath(os.path.join(script_dir,dir_name))
            
        if not os.path.isdir(dir_to_save):
            os.makedirs(dir_to_save)
        #TODO:add directory name
        filename = '{}\Results_brom{}.png'.format(dir_to_save,self.name)       
        #plt.savefig(results_dir+title+'.png')
        plt.savefig(filename, format='png', dpi=300, transparent=True)
           
    def plot_3fig(self):       
        plt.clf() 
             
        self.fh_water =  Dataset(self.water_fname)     
        self.depth_water = np.array(self.fh_water.variables['z_faces'][:])        
        self.min_water = np.amin(self.depth_water)
        self.max_water = np.amax(self.depth_water)
                  
        #ice_time_format  = num2date(ice_time) 
        #print (ice_time_format[0:10])
        try:
            self.time = self.fh_water.variables['time']      
            self.time2 = self.fh_water.variables['time'][:]
            self.time_units = self.fh_water.variables['time'].units
        except KeyError:
            self.time = self.fh_water.variables['ocean_time']   
            self.time2 = self.fh_water.variables['ocean_time'][:]            
            self.time_units = self.fh_water.variables['ocean_time'].units
        self.format_time = num2date(self.time2,units = self.time_units,calendar= 'standard')
         
        #########################
        # Values for time axis  #
        #########################
                
        start_year = self.combobox_start_year.value()
        stop_year = self.combobox_stop_year.value()
        
        to_start = datetime.datetime(start_year,1,1,12,0)
        to_stop= datetime.datetime(stop_year,1,1,12,0)
                
        start = date2index(to_start, self.time,
                            calendar=None, select='nearest')
        stop = date2index(to_stop, self.time,
                            calendar=None, select='nearest')           
        var_water,data_units = self.read_var()

 
        X_water,Y_water = np.meshgrid(self.time2[start:stop],self.depth_water)
        X_water = num2date(X_water,units = self.time_units)   
    
        ice_data = Dataset('Data\Laptev_average_year_3year.nc')    
        ice = ice_data.variables['hice'][:]        
        ice_time = ice_data.variables['time']   
               
        start_ice = date2index(to_start, ice_time,
                            calendar=None, select='nearest')
        stop_ice = date2index(to_stop, ice_time,
                            calendar=None, select='nearest')    
        print (start,start_ice,stop,stop_ice)    
        
        units = ice_data.variables['time'].units
        form_ice_time = num2date(ice_time[start_ice:stop_ice],units = units)
        print ('form_ice_time made ')       
        self.fh_water.close()
        ice_data.close()
        
        gs = gridspec.GridSpec(2, 2,width_ratios=[17, 1], height_ratios=[1, 3])
        gs.update(left=0.15, right= 0.85,top = 0.931,bottom = 0.11,
                           wspace=0.01,hspace=0.01)      
        #add subplots
        ax0 = self.figure.add_subplot(gs[0]) # o2
        #cbax0 = self.figure.add_subplot(gs[1])   
        ax1 = self.figure.add_subplot(gs[2]) 
        cbax1 = self.figure.add_subplot(gs[3])                  
        cmap_water = plt.get_cmap('RdBu_r') 
        
        #TODO: add ice here
        #without interpolation 
        ###CS1 = ax0.pcolor(X,Y,var_ice[:,start:stop],cmap = cmap )#) 3,edgecolor = 'w',
        ####                 #linewidth = 0.000005)
                        
        if self.checkbox_title.isChecked() == True:
            title = self.change_title.text()
            ax0.set_title(title)
        else:                 
            ax0.set_title('MI, '+self.long_name+' ['+ str(data_units)+']')
                   
        def fmt(x, pos):
            a, b = '{:.2e}'.format(x).split('e')
            b = int(b)
            return r'${} \times 10^{{{}}}$'.format(a, b)

     
        CS1 = ax0.stackplot(form_ice_time,ice[start_ice:stop_ice],color = '#7B98A8')
        ax0.set_xlim(to_start,to_stop)   
        ax0.set_ylim(0,max(ice[start_ice:stop_ice]))     
 
        mm = np.max(var_water[:,start:stop])
        mm2 = abs(np.min(var_water[:,start:stop]))
        mm_tot = max(mm,mm2)

        bounds = np.linspace(-mm_tot, mm_tot, 40)
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)        
        
        print (max(mm,mm2))
        CS4 = ax1.pcolor(X_water,Y_water,var_water[:,start:stop],norm = norm, #vmin=-mm_tot, vmax=mm_tot,
                         cmap = cmap_water) 
        from matplotlib.ticker import FormatStrFormatter           
        def add_colorbar(CS,axis,ma1):
            print (ma1)
            if abs(ma1) > 10000 or abs(ma1) < 0.01:
                cb = plt.colorbar(CS,cax = axis,
                format= ticker.FuncFormatter(fmt))
            else: 
                cb = plt.colorbar(CS,cax = axis,
                format = FormatStrFormatter('%.2f'))
            return cb
        

        ma1 = ma.max(var_water[:,start:stop])
        cb1 = add_colorbar(CS4,cbax1,ma1)
              
        dtime = stop-start  
        ### Time ticks ###         
        if dtime >= 367:
            dt =  int((stop - start)/365) #number of years
            time_ticks = [self.format_time[start]+relativedelta(years = n) for n in range(0,dt+1)]   
            

              
        labels = ["Ice thickness \n(m)", "Depth \n(m)"]
        n = 0                        
        for axis in (ax0,ax1): 
            try:
                axis.set_xticks(time_ticks)
            except: NameError
            axis.yaxis.set_label_coords(-0.09, 0.3)
            axis.set_ylabel(labels[n], fontsize = self.fontsize)              
            n=n+1
            

        ax1.set_ylim(self.max_water,self.min_water)
        ax0.set_yticks(np.arange(0.5,2.5,0.5))  
        ax0.set_xticklabels([])  
            
        if dtime >= 367:
            ax1.set_xticklabels(['Jan 1','Jan 2','Jan 3','Jan 4'])  
            for n in time_ticks:  
                ax1.axvline(x=n, ls= '--',  lw = 1 )  
        elif dtime  < 367:
            form = '%b' 
            ax1.xaxis.set_major_formatter(
                mdates.DateFormatter(form)) 
                                           
        if self.action == 'showfigure' : 
            self.canvas.draw()                
        elif self.action == 'savepdf' :
            self.canvas.draw()
            self.save_to_dir('Figures_from_GUI')
                       
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("plastique")
    main = Window()
    main.setStyleSheet("background-color:#dceaed;")
    main.show()  
    sys.exit(app.exec_())    