import numpy as np
import os
import shutil
from scipy import signal
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import obspy
import pyflex


class windows:
    def __init__(self,par):
        # Parameters
        self.par = par
        # Reading WorkDir
        os.chdir(self.par.InvPar["-WorkDir"])

        self.winpar = {}
        # Read ConfigFile
        _winpar = np.genfromtxt('PyFlex.par',dtype='str',
			comments='#',delimiter='=')

        self.event = None
        self.station = None
        self.windows = None

        for i in range(0,_winpar.shape[0]):
	    self.winpar.update({_winpar[i][0].strip():_winpar[i][1].strip()})


        min_period = 1.0 / self.par.freq[-1]
        max_period = 1.0 / self.par.freq[0]
        stalta_waterlevel = float(self.winpar["stalta_waterlevel"])
        tshift_acceptance_level = float(self.winpar["tshift_acceptance_level"])
        tshift_reference = float(self.winpar["tshift_reference"])
        dlna_acceptance_level = float(self.winpar["dlna_acceptance_level"])
        dlna_reference = float(self.winpar["dlna_reference"])
        cc_acceptance_level = float(self.winpar["cc_acceptance_level"])
        s2n_limit = float(self.winpar["s2n_limit"])
        earth_model = str(self.winpar["earth_model"])
        min_surface_wave_velocity = float(self.winpar["min_surface_wave_velocity"])
        max_time_before_first_arrival = float(self.winpar["max_time_before_first_arrival"])
        c_0 = float(self.winpar["c_0"])
        c_1 = float(self.winpar["c_1"])
        c_2 = float(self.winpar["c_2"])
        c_3a = float(self.winpar["c_3a"])
        c_3b = float(self.winpar["c_3b"])
        c_4a = float(self.winpar["c_4a"])
        c_4b = float(self.winpar["c_4b"])
        check_global_data_quality = eval(self.winpar["check_global_data_quality"])
        snr_integrate_base = float(self.winpar["snr_integrate_base"])
        snr_max_base = float(self.winpar["snr_max_base"])
        noise_start_index = eval(self.winpar["noise_start_index"])
        noise_end_index = eval(self.winpar["noise_end_index"])
        signal_start_index = eval(self.winpar["signal_start_index"])
        signal_end_index = -1
        window_weight_fct = None
        window_signal_to_noise_type = str(self.winpar["window_signal_to_noise_type"])
        resolution_strategy = str(self.winpar["resolution_strategy"])

        self.config = pyflex.Config( min_period = min_period\
                               ,max_period = max_period \
                               ,stalta_waterlevel= stalta_waterlevel \
                                ,tshift_acceptance_level = tshift_acceptance_level \
                                , tshift_reference= tshift_reference \
                                , dlna_acceptance_level=dlna_acceptance_level \
                                , dlna_reference=dlna_reference \
                                , cc_acceptance_level= cc_acceptance_level \
                                , s2n_limit= s2n_limit \
                                , earth_model= earth_model\
                                , min_surface_wave_velocity= min_surface_wave_velocity \
                                , max_time_before_first_arrival= max_time_before_first_arrival \
                                , c_0= c_0 \
                                , c_1= c_1 \
                                , c_2= c_2 \
                                , c_3a= c_3a \
                                , c_3b= c_3b \
                                , c_4a= c_4a \
                                , c_4b= c_4b \
                                , check_global_data_quality=check_global_data_quality \
                                , snr_integrate_base= snr_integrate_base\
                                , snr_max_base= snr_max_base \
                                , noise_start_index= noise_start_index \
                                , noise_end_index= noise_end_index\
                                , signal_start_index= signal_start_index \
                                , signal_end_index= signal_end_index \
                                , window_weight_fct= window_weight_fct \
                                , window_signal_to_noise_type= window_signal_to_noise_type \
                                , resolution_strategy= resolution_strategy)



    def CreateInfoStatEvent(self,data,timeoffset):
        self.event = pyflex.Event(latitude=data[0].stats.sac.evla, longitude=data[0].stats.sac.evlo, \
                             depth_in_m=data[0].stats.sac.evdp * 1000,\
                             origin_time=data[0].stats.starttime + timeoffset)
        
        self.station = pyflex.Station(latitude=data[0].stats.sac.stla, longitude=data[0].stats.sac.stlo)
        
    def CreateWindows(self,data,syn,filename="hola",plot=True):
        config = self.config
        self.windows = pyflex.select_windows(data, syn, config, \
                                             event=self.event, station =self.station, \
                                             windows_filename= filename + ".json", plot=plot)

    def LoadWindows(self,data,syn,filename):
        config = self.config
        ws2 = pyflex.WindowSelector(data, syn, config)
        print("Windows before loading: %i" % len(ws2.windows))
        ws2.load(filename + ".json")
        print("Windows after loading: %i" % len(ws2.windows))
        self.windows = ws2.windows
        


        
        
        
