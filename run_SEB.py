#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Processes the input and 
"""
import core # Library of SEB computation
import seb_utils # Library of helper functions (e.g. convert met vars)

import pandas as pd
import numpy as np
import datetime

# Options/Parameters
din="/home/lunet/gytm3/Everest2019/AWS/Logging/"
fs=["south_col.csv","c2.csv","summit_guess.csv"]
end_date=datetime.datetime(year=2019,month=11,day=1)

# Preallocation
out={}
input_data={}
        
for jj in range(len(fs)):
            
            # Update
            print "Computing SEB for: ",fs[jj]

            # Read in data
            fin=din+fs[jj]
            data=pd.read_csv(fin,sep=",",parse_dates=True, index_col=0,\
                             na_values="NAN")

            # Tuncate to ensure no more recent than 'end_date'
            data=data.loc[data.index<=end_date]
            
            # Resample to higher frequency --
            freq="%.0fmin" % (core.ds/60.)
            data=data.resample(freq).interpolate("linear")
                    
            # Assignments (for readability)
            ta=data["T_HMP"].values[:]+273.15
            ta2=data["T_109"].values[:]+273.15
            p=data["PRESS"].values[:]*100.
            rh=data["RH"].values[:]/100.      
            sw=data["SW_IN_AVG"].values[:].astype(np.float)
            sw_o=data["SW_OUT_AVG"].values[:].astype(np.float)
            lw=data["LW_IN_AVG"].values[:]
            
            # Fill missing pressure with LTM
            p[np.isnan(p)]=np.nanmean(p)
            # Wind speed names differ between files - deal with that here. 
            try:
                u=(data["WS_AVG"].values[:]+data["WS_AVG_2"].values[:])/2.
            except:
                u=data["WS_AVG"].values[:]
                
            # Met conversions
            vp=np.squeeze(np.array([seb_utils.SATVP(ta[i])*rh[i] for i in range(len(u))]))
            mr=seb_utils.MIX(p,vp)
            tv=seb_utils.VIRTUAL(ta,mr)
            qa=seb_utils.VP2Q(vp,p)
            rho=seb_utils.RHO(p,tv)
            # Clean input data - at the moment skipping computation when 
            # suspicious values found
            # NaN where big disagreement between temps
            ta_err_idx=np.abs(ta-ta2)>10.
            ta[ta_err_idx]=np.nan
            # Zero Sw when toa is = 0
            # Compute ToA 
            toa,r=seb_utils.sin_toa(data.index.dayofyear.values[:],\
                        data.index.hour.values[:],27.98,86.93)
            sw[toa==0]=0.0
            sw_o[toa==0]=0.0
            
            # Find where SW_IN < SW_OUT. Replace with 24-hour-centred albedo-
            # -based estimate. The albedo-based estimate assumes a maximum 
            # possible albedo of 0.90
            run_albedo=data["SW_OUT_AVG"].rolling(24*30).sum()/\
            data["SW_IN_AVG"].rolling(24*30).sum()
            run_albedo.loc[run_albedo>0.90]=0.90
            sw_idx=sw<sw_o
            sw[sw_idx]=sw_o[sw_idx]/run_albedo.loc[sw_idx]
            
            # Simulate LW -- will be reference for periods when we assume 
            # snow is on the upward pyranometer
            if "summit" not in fs[jj]:
                lw_ref=seb_utils.sim_lw(rh*100.,ta+273.15,lw,toa)
                lw[sw_idx]=lw_ref[sw_idx] # -> slot in the best guess   
            
            # Summarise data being fed to SEB computations
            print "Invalid T = %.0f" % (np.sum(ta_err_idx))  
            print r"Invalid SW = %.0f%%" % (np.sum(sw_idx)/\
                                           np.float(len(sw_idx))*100)    
            print "SW mean: ", np.nanmean(sw)
                                                           
            # Compute the SEB 
            shf,lhf,swn,lwn,seb,ts,melt,sub,nit_log,qg = \
            core.SEB(ta,qa,rho,u,p,sw,lw,ta[0])
               
            # Store in DataFrame for easier computation/plotting 
            # indices 0: SCol and 1: Camp II
            out[jj]=\
            pd.DataFrame({"Melt":melt,"Ablation":melt+sub,"T":ta,\
                          "shf":shf,"lhf":lhf,"sw":swn,"lw":lwn,\
                          "seb":seb,"qg":qg,"qmelt":melt*core.Lf/core.ds,\
                          "wind":u,"rh":rh,"ts":ts,"sub":sub},index=data.index)
    
            # Store input data in case we need to do something with it
            data["run_albedo"]=run_albedo
            input_data[fs[jj].replace(".csv","")]=data
            
            # Update
            print("Computed SEB for file: %s" %fin )