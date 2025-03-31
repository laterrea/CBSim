"""
- Date: March 2025
- Author: Antoine Laterre

@ Project: "CBSim - Carnot Battery Simulator" 

-> Module for post-processing (plots)

"""

#%% IMPORT PACKAGES 

import matplotlib.pyplot as plt
import numpy as np
import CoolProp.CoolProp as CP
from   CoolProp import AbstractState
import CoolProp
import warnings

#%% VCHP

def plot_Th_VCHP(HP):
    
    # Process inputs
    HP.export_states()
    i_ref_hp = HP.i_hp_4
    i_ref_cs = HP.i_hp_cs_ex
    i_ref_hs = HP.i_hp_hs_su
    
    cell_res_long   = 0.05 # K
    size_space_long = int(1/cell_res_long*(HP.T_hp_hs_ex-HP.T_hp_hs_su)+1)
    cell_res_shrt   = 0.10 # K
    size_space_shrt = int(1/cell_res_shrt*(HP.T_hp_cs_su-HP.T_hp_cs_ex)+1)
    
    try:    recup, i_hp_1r, i_hp_3r = True, HP.i_hp_1r, HP.i_hp_3r
    except: recup                   = False
    if HP.p_hp_2 >= CP.PropsSI('P_CRITICAL',HP.parameters['fluid_hp']): trans_crit = True
    else:                                                               trans_crit = False
    
    # Saturation Curve
    T_min = max(min(HP.T_hp)-5,CP.PropsSI('T_MIN',HP.parameters['fluid_hp'])+1)
    T_max = min(max(HP.T_hp)+1,CP.PropsSI('T_MAX',HP.parameters['fluid_hp'])-1,CP.PropsSI('T_CRITICAL',HP.parameters['fluid_hp'])-1)
    temp = np.linspace(T_min,T_max,1000)
    h_lq = np.zeros(len(temp))
    h_vp = np.zeros(len(temp))
    
    my_state_hp = AbstractState("HEOS",HP.parameters['fluid_hp'])
    for i, t in enumerate(temp):
        my_state_hp.update(CoolProp.QT_INPUTS,1,t)
        h_vp[i] = my_state_hp.hmass()-i_ref_hp
        my_state_hp.update(CoolProp.QT_INPUTS,0,t)
        h_lq[i] = my_state_hp.hmass()-i_ref_hp
    
    fig = plt.figure()
    ax1 = fig.subplots(1)
    ax1.plot(h_lq*1e-3,temp-273.15,'-',color='black',linewidth=1,label=HP.parameters['fluid_hp'])
    ax1.plot(h_vp*1e-3,temp-273.15,'-',color='black',linewidth=1)
    ax1.plot(5e-4*(h_lq[len(temp)-1]+h_vp[len(temp)-1]),T_max-273.15,'.',color='black',linewidth=1)
    
    # Heat Pump Diagram
    if recup:
        # Recuperator: cold leg
        hrelt    = np.linspace(HP.i_hp_1,HP.i_hp_1r,size_space_long)
        prelt    = np.linspace(HP.p_hp_1,HP.p_hp_1r,size_space_long)
        Trelt    = CP.PropsSI('T','P',prelt,'H',hrelt,HP.parameters['fluid_hp'])
        ax1.plot((hrelt-i_ref_hp)*1e-3,Trelt-273.15,'--',color='firebrick',linewidth=3)
        # Compression
        pcomp    = np.linspace(HP.p_hp_1r,HP.p_hp_2,size_space_shrt)
        hcomp_is = CP.PropsSI('H','P',pcomp,'S',HP.s_hp_1r,HP.parameters['fluid_hp'])
        hcomp    = HP.i_hp_1r+(hcomp_is-HP.i_hp_1r)/HP.eta_is_cp
        # Recuperator: cold leg
        hreht    = np.linspace(HP.i_hp_3,HP.i_hp_3r,size_space_long)
        preht    = np.linspace(HP.p_hp_3,HP.p_hp_3r,size_space_long)
        Treht    = CP.PropsSI('T','P',preht,'H',hreht,HP.parameters['fluid_hp'])
        ax1.plot((hreht-i_ref_hp)*1e-3,Treht-273.15,'--',color='firebrick',linewidth=3)
        # Expansion
        hexpv    = np.linspace(HP.i_hp_3r,HP.i_hp_4,size_space_shrt)
        pexpv    = np.linspace(HP.p_hp_3r,HP.p_hp_4,size_space_shrt)
    else:
        # Compression
        pcomp    = np.linspace(HP.p_hp_1,HP.p_hp_2,size_space_shrt)
        hcomp_is = CP.PropsSI('H','P',pcomp,'S',HP.s_hp_1,HP.parameters['fluid_hp'])
        hcomp    = HP.i_hp_1+(hcomp_is-HP.i_hp_1)/HP.eta_is_cp
        # Expansion
        hexpv    = np.linspace(HP.i_hp_3,HP.i_hp_4,size_space_shrt)
        pexpv    = np.linspace(HP.p_hp_3,HP.p_hp_4,size_space_shrt)
    # Compression
    Tcomp    = CP.PropsSI('T','P',pcomp,'H',hcomp,    HP.parameters['fluid_hp'])    
    ax1.plot((hcomp-i_ref_hp)*1e-3,Tcomp-273.15,'-',color='firebrick',linewidth=3,label='heat pump cycle')
    
    # Condensation
    if trans_crit:
        hcond    = np.linspace(HP.i_hp_3,HP.i_hp_2,2*size_space_long)
        pcond    = np.linspace(HP.p_hp_3,HP.p_hp_2,2*size_space_long)
    else:
        hcond_lq = np.linspace(HP.i_hp_3, HP.i_hp_3x,size_space_long)
        pcond_lq = np.linspace(HP.p_hp_3, HP.p_hp_3x,size_space_long)
        hcond_st = np.linspace(HP.i_hp_3x,HP.i_hp_2x,size_space_shrt)
        pcond_st = np.linspace(HP.p_hp_3x,HP.p_hp_2x,size_space_shrt)
        hcond_vp = np.linspace(HP.i_hp_2x,HP.i_hp_2, size_space_long)
        pcond_vp = np.linspace(HP.p_hp_2x,HP.p_hp_2, size_space_long)
        hcond    = np.concatenate((hcond_lq,hcond_st,hcond_vp),axis=None)
        pcond    = np.concatenate((pcond_lq,pcond_st,pcond_vp),axis=None)
    Tcond    = CP.PropsSI('T','P',pcond,'H',hcond,    HP.parameters['fluid_hp'])
    ax1.plot((hcond-i_ref_hp)*1e-3,Tcond-273.15,'-',color='firebrick',linewidth=3)
    # Expansion
    Texpv    = CP.PropsSI('T','P',pexpv,'H',hexpv,    HP.parameters['fluid_hp'])
    ax1.plot((hexpv-i_ref_hp)*1e-3,Texpv-273.15,'-',color='firebrick',linewidth=3)
    # Evaporation
    hevap_lq = np.linspace(HP.i_hp_4, HP.i_hp_4x,size_space_long)
    pevap_lq = np.linspace(HP.p_hp_4, HP.p_hp_4x,size_space_long)
    hevap_st = np.linspace(HP.i_hp_4x,HP.i_hp_1x,size_space_shrt)
    pevap_st = np.linspace(HP.p_hp_4x,HP.p_hp_1x,size_space_shrt)
    hevap_vp = np.linspace(HP.i_hp_1x,HP.i_hp_1, size_space_long)
    pevap_vp = np.linspace(HP.p_hp_1x,HP.p_hp_1, size_space_long)
    hevap    = np.concatenate((hevap_lq,hevap_st,hevap_vp),axis=None)
    pevap    = np.concatenate((pevap_lq,pevap_st,pevap_vp),axis=None)
    Tevap    = CP.PropsSI('T','P',pevap,'H',hevap,    HP.parameters['fluid_hp'])
    ax1.plot((hevap-i_ref_hp)*1e-3,Tevap-273.15,'-',color='firebrick',linewidth=3)
    # Mark the points
    ax1.plot((np.asarray(HP.i_hp)-i_ref_hp)*1e-3,
              np.asarray(HP.T_hp          )-273.15,
              'o',color='firebrick',linewidth=4)
    
    # Secondary Fluids Diagram
    pcond_sf   = np.linspace(HP.p_hp_hs_su,HP.p_hp_hs_ex,len(hcond))
    hcond_sf   = np.zeros(len(hcond))
    hcond_sf[0]= HP.i_hp_hs_su
    fra_hs     =(HP.i_hp_2-HP.i_hp_3)/(HP.i_hp_hs_ex-HP.i_hp_hs_su) 
    for i in range(len(hcond_sf)-1):
        hcond_sf[i+1] = hcond_sf[i]+(hcond[i+1]-hcond[i])/fra_hs
    Tcond_sf = CP.PropsSI('T','P',pcond_sf,'H',hcond_sf,HP.fluid_hp_hs)
    if recup:   hcond_sf = (hcond_sf-i_ref_hs)*(max(hcond)-min(hcond))/(max(hcond_sf)-min(hcond_sf))\
                         + (HP.i_hp_3-HP.i_hp_3r)
    else:       hcond_sf = (hcond_sf-i_ref_hs)*(max(hcond)-min(hcond))/(max(hcond_sf)-min(hcond_sf))
    # ax1.plot(hcond_sf*1e-3,Tcond_sf-273.15,'--',color='darkgreen',linewidth=1,label='thermal storage')
    ax1.plot(hcond_sf*1e-3,Tcond_sf-273.15,'--',color='darkgreen',linewidth=1,label='heat sink') # !!!
    ax1.plot(np.asarray([hcond_sf[0],hcond_sf[-1]])*1e-3,
             np.asarray([Tcond_sf[0],Tcond_sf[-1]])-273.15,
             'x',color='darkgreen',linewidth=4)
    
    pevap_sf   = np.linspace(HP.p_hp_cs_ex,HP.p_hp_cs_su,len(hevap))
    hevap_sf   = np.zeros(len(hevap))
    hevap_sf[0]= HP.i_hp_cs_ex
    fra_cs     =(HP.i_hp_1-HP.i_hp_4)/(HP.i_hp_cs_su-HP.i_hp_cs_ex) 
    for i in range(len(hevap_sf)-1):
        hevap_sf[i+1] = hevap_sf[i]+(hevap[i+1]-hevap[i])/fra_cs
    Tevap_sf = CP.PropsSI('T','P',pevap_sf,'H',hevap_sf,HP.fluid_hp_cs)
    hevap_sf = (hevap_sf-i_ref_cs)*(max(hevap)-min(hevap))/(max(hevap_sf)-min(hevap_sf))
    ax1.plot(hevap_sf*1e-3,Tevap_sf-273.15,'--',color='darkslategrey',linewidth=1,label='heat source')
    ax1.plot(np.asarray([hevap_sf[0],hevap_sf[-1]])*1e-3,
             np.asarray([Tevap_sf[0],Tevap_sf[-1]])-273.15,
             'x',color='darkslategrey',linewidth=4)
    
    # Check if pinch is violated
    pp_cs = np.argmin(Tevap_sf-Tevap)
    pp_hs = np.argmin(Tcond-Tcond_sf)
    
    ax1.plot(hevap_sf[pp_cs]*1e-3,Tevap[pp_cs]+HP.parameters['dT_hp_ev_pp']-273.15, 'x', color='darkblue',     linewidth=4)
    ax1.plot(hcond_sf[pp_hs]*1e-3,Tcond[pp_hs]-HP.parameters['dT_hp_cd_pp']-273.15, 'x', color='darkblue',     linewidth=4)
    if round(Tevap_sf[pp_cs]-Tevap[pp_cs],2) < HP.parameters['dT_hp_ev_pp']: 
        warnings.warn(f'In the VCHP evaporator, the pinch is violated: {Tevap_sf[pp_cs]-Tevap[pp_cs]} < {HP.parameters["dT_hp_ev_pp"]}', Warning)
    if round(Tcond[pp_hs]-Tcond_sf[pp_hs],2) < HP.parameters['dT_hp_cd_pp']:
        warnings.warn(f'In the VCHP condenser, the pinch is violated: {Tcond[pp_hs]-Tcond_sf[pp_hs]} < {HP.parameters["dT_hp_cd_pp"]}', Warning)
    
    # Axes
    ax1.legend(loc=(0.05,0.75),frameon=False,labelcolor='linecolor')
    ax1.text(0.055,0.72,'$\mathrm{COP_{HP}}$'+f' = {round(HP.eta_hp_cyclen,2)} [-]',transform=plt.gca().transAxes)
    ax1.set_xlabel('enthalpy [kJ/kg]',fontsize=14, color='dimgray')
    ax1.spines["bottom"].set_bounds(h_lq[0]*1e-3,hcond_sf[-1]*1e-3)
    if recup: ax1.set_xticks(np.round([h_lq[0]*1e-3,hevap_sf[0]*1e-3,hcond_sf[0]*1e-3,hevap_sf[-1]*1e-3,hcond_sf[-1]*1e-3],1))
    else:     ax1.set_xticks(np.round([h_lq[0]*1e-3,hevap_sf[0]*1e-3,                 hevap_sf[-1]*1e-3,hcond_sf[-1]*1e-3],1))
    ax1.spines["bottom"].set_linewidth(1.5)
    ax1.spines["bottom"].set_color('dimgray')
    ax1.set_ylabel('temperature [°C]',fontsize=14, color='dimgray')
    ax1.spines["left"].set_bounds(min(HP.T_hp)-273.15,max(HP.T_hp)-273.15)
    ax1.set_yticks(np.round([min(HP.T_hp)-273.15,HP.T_hp_3-273.15,Tcond_sf[0]-273.15,
                    Tcond_sf[-1]-273.15,Tevap_sf[0]-273.15,Tevap_sf[-1]-273.15,
                    max(HP.T_hp)-273.15],1))
    ax1.spines["left"].set_linewidth(1.5)
    ax1.spines["left"].set_color('dimgray')
    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)
    ax1.tick_params(axis='both', labelsize=12, colors='dimgray', width=1.5, length=8, direction='in')
    ax1.label_outer()
    fig.tight_layout()
    plt.show()
    out = fig
    return out

def plot_Ts_VCHP(HP):
    
    # Process inputs
    HP.export_states()
    i_ref_hp = HP.i_hp_4
    i_ref_cs = HP.i_hp_cs_ex
    i_ref_hs = HP.i_hp_hs_su
    
    cell_res_long   = 0.05 # K
    size_space_long = int(1/cell_res_long*(HP.T_hp_hs_ex-HP.T_hp_hs_su)+1)
    cell_res_shrt   = 0.10 # K
    size_space_shrt = int(1/cell_res_shrt*(HP.T_hp_cs_su-HP.T_hp_cs_ex)+1)
    
    try:    recup, i_hp_1r, i_hp_3r = True, HP.i_hp_1r, HP.i_hp_3r
    except: recup                   = False
    if HP.p_hp_2 >= CP.PropsSI('P_CRITICAL',HP.parameters['fluid_hp']): trans_crit = True
    else:                                                               trans_crit = False
    
    # Saturation Curve
    T_min = max(min(HP.T_hp)-5,CP.PropsSI('T_MIN',HP.parameters['fluid_hp'])+1)
    T_max = min(max(HP.T_hp)+1,CP.PropsSI('T_MAX',HP.parameters['fluid_hp'])-1,CP.PropsSI('T_CRITICAL',HP.parameters['fluid_hp'])-1)
    temp = np.linspace(T_min,T_max,1000)
    s_lq = np.zeros(len(temp))
    s_vp = np.zeros(len(temp))
    
    my_state_hp = AbstractState("HEOS",HP.parameters['fluid_hp'])
    for i, t in enumerate(temp):
        my_state_hp.update(CoolProp.QT_INPUTS,1,t)
        s_vp[i] = my_state_hp.smass()
        my_state_hp.update(CoolProp.QT_INPUTS,0,t)
        s_lq[i] = my_state_hp.smass()
    
    fig = plt.figure()
    ax1 = fig.subplots(1)
    ax1.plot(s_lq*1e-3,temp-273.15,'-',color='black',linewidth=1,label=HP.parameters['fluid_hp'])
    ax1.plot(s_vp*1e-3,temp-273.15,'-',color='black',linewidth=1)
    ax1.plot(5e-4*(s_lq[len(temp)-1]+s_vp[len(temp)-1]),T_max-273.15,'.',color='black',linewidth=1)
    
    # Heat Pump Diagram
    if recup:
        # Recuperator: cold leg
        hrelt    = np.linspace(HP.i_hp_1,HP.i_hp_1r,size_space_long)
        prelt    = np.linspace(HP.p_hp_1,HP.p_hp_1r,size_space_long)
        Trelt    = CP.PropsSI('T','P',prelt,'H',hrelt,HP.parameters['fluid_hp'])
        srelt    = CP.PropsSI('S','P',prelt,'H',hrelt,HP.parameters['fluid_hp'])
        ax1.plot(srelt*1e-3,Trelt-273.15,'--',color='firebrick',linewidth=3)
        # Compression
        pcomp    = np.linspace(HP.p_hp_1r,HP.p_hp_2,size_space_shrt)
        hcomp_is = CP.PropsSI('H','P',pcomp,'S',HP.s_hp_1r,HP.parameters['fluid_hp'])
        hcomp    = HP.i_hp_1r+(hcomp_is-HP.i_hp_1r)/HP.eta_is_cp
        # Recuperator: cold leg
        hreht    = np.linspace(HP.i_hp_3,HP.i_hp_3r,size_space_long)
        preht    = np.linspace(HP.p_hp_3,HP.p_hp_3r,size_space_long)
        Treht    = CP.PropsSI('T','P',preht,'H',hreht,HP.parameters['fluid_hp'])
        sreht    = CP.PropsSI('S','P',preht,'H',hreht,HP.parameters['fluid_hp'])
        ax1.plot(sreht*1e-3,Treht-273.15,'--',color='firebrick',linewidth=3)
        # Expansion
        hexpv    = np.linspace(HP.i_hp_3r,HP.i_hp_4,size_space_shrt)
        pexpv    = np.linspace(HP.p_hp_3r,HP.p_hp_4,size_space_shrt)
    else:
        # Compression
        pcomp    = np.linspace(HP.p_hp_1,HP.p_hp_2,size_space_shrt)
        hcomp_is = CP.PropsSI('H','P',pcomp,'S',HP.s_hp_1,HP.parameters['fluid_hp'])
        hcomp    = HP.i_hp_1+(hcomp_is-HP.i_hp_1)/HP.eta_is_cp
        # Expansion
        hexpv    = np.linspace(HP.i_hp_3,HP.i_hp_4,size_space_shrt)
        pexpv    = np.linspace(HP.p_hp_3,HP.p_hp_4,size_space_shrt)
    # Compression
    Tcomp    = CP.PropsSI('T','P',pcomp,'H',hcomp,    HP.parameters['fluid_hp'])
    scomp    = CP.PropsSI('S','P',pcomp,'H',hcomp,    HP.parameters['fluid_hp'])
    ax1.plot(scomp*1e-3,Tcomp-273.15,'-',color='firebrick',linewidth=3,label='heat pump cycle')
    
    # Condensation
    if trans_crit:
        hcond    = np.linspace(HP.i_hp_3,HP.i_hp_2,2*size_space_long)
        pcond    = np.linspace(HP.p_hp_3,HP.p_hp_2,2*size_space_long)
    else:
        hcond_lq = np.linspace(HP.i_hp_3, HP.i_hp_3x,size_space_long)
        pcond_lq = np.linspace(HP.p_hp_3, HP.p_hp_3x,size_space_long)
        hcond_st = np.linspace(HP.i_hp_3x,HP.i_hp_2x,size_space_shrt)
        pcond_st = np.linspace(HP.p_hp_3x,HP.p_hp_2x,size_space_shrt)
        hcond_vp = np.linspace(HP.i_hp_2x,HP.i_hp_2, size_space_long)
        pcond_vp = np.linspace(HP.p_hp_2x,HP.p_hp_2, size_space_long)
        hcond    = np.concatenate((hcond_lq,hcond_st,hcond_vp),axis=None)
        pcond    = np.concatenate((pcond_lq,pcond_st,pcond_vp),axis=None)
    Tcond    = CP.PropsSI('T','P',pcond,'H',hcond,    HP.parameters['fluid_hp'])
    scond    = CP.PropsSI('S','P',pcond,'H',hcond,    HP.parameters['fluid_hp'])
    ax1.plot(scond*1e-3,Tcond-273.15,'-',color='firebrick',linewidth=3)
    # Expansion
    Texpv    = CP.PropsSI('T','P',pexpv,'H',hexpv,    HP.parameters['fluid_hp'])
    sexpv    = CP.PropsSI('S','P',pexpv,'H',hexpv,    HP.parameters['fluid_hp'])
    ax1.plot(sexpv*1e-3,Texpv-273.15,'-',color='firebrick',linewidth=3)
    # Evaporation
    hevap_lq = np.linspace(HP.i_hp_4, HP.i_hp_4x,size_space_long)
    pevap_lq = np.linspace(HP.p_hp_4, HP.p_hp_4x,size_space_long)
    hevap_st = np.linspace(HP.i_hp_4x,HP.i_hp_1x,size_space_shrt)
    pevap_st = np.linspace(HP.p_hp_4x,HP.p_hp_1x,size_space_shrt)
    hevap_vp = np.linspace(HP.i_hp_1x,HP.i_hp_1, size_space_long)
    pevap_vp = np.linspace(HP.p_hp_1x,HP.p_hp_1, size_space_long)
    hevap    = np.concatenate((hevap_lq,hevap_st,hevap_vp),axis=None)
    pevap    = np.concatenate((pevap_lq,pevap_st,pevap_vp),axis=None)
    Tevap    = CP.PropsSI('T','P',pevap,'H',hevap,    HP.parameters['fluid_hp'])
    sevap    = CP.PropsSI('S','P',pevap,'H',hevap,    HP.parameters['fluid_hp'])
    ax1.plot(sevap*1e-3,Tevap-273.15,'-',color='firebrick',linewidth=3)
    # Mark the points
    ax1.plot(np.asarray(HP.s_hp)*1e-3,
             np.asarray(HP.T_hp)-273.15,
             'o',color='firebrick',linewidth=4)
    
    # Axes
    ax1.legend(loc=(0.05,0.855),frameon=False,labelcolor='linecolor')
    ax1.text(0.055,0.82,'$\mathrm{COP_{HP}}$'+f' = {round(HP.eta_hp_cyclen,2)} [-]',transform=plt.gca().transAxes)
    ax1.set_xlabel('entropy [kJ/kg/K]',fontsize=14, color='dimgray')
    ax1.spines["bottom"].set_bounds(s_lq[0]*1e-3,scomp[-1]*1e-3)
    if recup: ax1.set_xticks(np.round([s_lq[0]*1e-3,HP.s_hp_3*1e-3,HP.s_hp_3r*1e-3,HP.s_hp_4*1e-3,HP.s_hp_2*1e-3],3))
    else:     ax1.set_xticks(np.round([s_lq[0]*1e-3,HP.s_hp_3*1e-3,                HP.s_hp_4*1e-3,HP.s_hp_2*1e-3],3))
    ax1.spines["bottom"].set_linewidth(1.5)
    ax1.spines["bottom"].set_color('dimgray')
    ax1.set_ylabel('temperature [°C]',fontsize=14, color='dimgray')
    ax1.spines["left"].set_bounds(min(HP.T_hp)-273.15,max(HP.T_hp)-273.15)
    if trans_crit:
        if recup: 
            ax1.set_yticks(np.round([min(HP.T_hp)-273.15,HP.T_hp_1r-273.15,
                                     HP.T_hp_3-273.15,HP.T_hp_1-273.15,
                                     HP.T_hp_3r-273.15,max(HP.T_hp)-273.15],1))
        else:
            ax1.set_yticks(np.round([min(HP.T_hp)-273.15,
                                     HP.T_hp_3-273.15,HP.T_hp_1-273.15,
                                     max(HP.T_hp)-273.15],1))
    else:
        if recup: 
            ax1.set_yticks(np.round([min(HP.T_hp)-273.15,HP.T_hp_1r-273.15,
                                     HP.T_hp_3-273.15,HP.T_hp_1-273.15,HP.T_hp_2x-273.15,
                                     HP.T_hp_3r-273.15,max(HP.T_hp)-273.15],1))
        else:
            ax1.set_yticks(np.round([min(HP.T_hp)-273.15,
                                     HP.T_hp_3-273.15,HP.T_hp_1-273.15,HP.T_hp_2x-273.15,
                                     max(HP.T_hp)-273.15],1))
    ax1.spines["left"].set_linewidth(1.5)
    ax1.spines["left"].set_color('dimgray')
    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)
    ax1.tick_params(axis='both', labelsize=12, colors='dimgray', width=1.5, length=8, direction='in')
    ax1.label_outer()
    fig.tight_layout()
    plt.show()
    out = fig
    return out

#%% ORC

def plot_Th_ORC(HE):
    
    # Process inputs
    HE.export_states()
    i_ref_he = HE.i_he_1
    i_ref_cs = HE.i_he_cs_su
    i_ref_hs = HE.i_he_hs_ex
    
    cell_res_long   = 0.05 # K
    size_space_long = int(1/cell_res_long*(HE.T_he_hs_su-HE.T_he_hs_ex)+1)
    cell_res_shrt   = 0.10 # K
    size_space_shrt = int(1/cell_res_shrt*(HE.T_he_cs_ex-HE.T_he_cs_su)+1)
    
    try:    recup, i_he_2r, i_he_4r = True, HE.i_he_2r, HE.i_he_4r
    except: recup                   = False
    if HE.p_he_2 >= CP.PropsSI('P_CRITICAL',HE.parameters['fluid_he']): trans_crit = True
    else:                                                               trans_crit = False
    
    # Saturation Curve
    T_min = max(min(HE.T_he)-2,CP.PropsSI('T_MIN',HE.parameters['fluid_he'])+1)
    # T_max = min(HE.T_he_hs_su, CP.PropsSI('T_MAX',HE.parameters['fluid_he'])-1,CP.PropsSI('T_CRITICAL',HE.parameters['fluid_he'])-1)
    T_max = min(max(HE.T_he)+5,CP.PropsSI('T_MAX',HE.parameters['fluid_he'])-1,CP.PropsSI('T_CRITICAL',HE.parameters['fluid_he'])-1)
    temp = np.linspace(T_min,T_max,1000)
    h_lq = np.zeros(len(temp))
    h_vp = np.zeros(len(temp))
    
    my_state_he = AbstractState("HEOS",HE.parameters['fluid_he'])
    for i, t in enumerate(temp):
        my_state_he.update(CoolProp.QT_INPUTS,1,t)
        h_vp[i] = my_state_he.hmass()-i_ref_he
        my_state_he.update(CoolProp.QT_INPUTS,0,t)
        h_lq[i] = my_state_he.hmass()-i_ref_he
        
    fig = plt.figure()
    ax1 = fig.subplots(1)
    ax1.plot(h_lq*1e-3,temp-273.15,'-',color='black',linewidth=1,label=HE.parameters['fluid_he'])
    ax1.plot(h_vp*1e-3,temp-273.15,'-',color='black',linewidth=1)
    ax1.plot(5e-4*(h_lq[len(temp)-1]+h_vp[len(temp)-1]),T_max-273.15,'.',color='black',linewidth=1)
    
    # Heat Engine Diagram
    if recup:
        # Recuperator: cold leg
        hrelt    = np.linspace(HE.i_he_2,HE.i_he_2r,size_space_long)
        prelt    = np.linspace(HE.p_he_2,HE.p_he_2r,size_space_long)
        Trelt    = CP.PropsSI('T','P',prelt,'H',hrelt,HE.parameters['fluid_he'])
        ax1.plot((hrelt-i_ref_he)*1e-3,Trelt-273.15,'--',color='darkblue',linewidth=3)
        # Evporator        
        if trans_crit:
            hevap    = np.linspace(HE.i_he_2r,HE.i_he_3,2*size_space_long)
            pevap    = np.linspace(HE.p_he_2r,HE.p_he_3,2*size_space_long)
        else:
            hevap_lq = np.linspace(HE.i_he_2r,HE.i_he_2x,size_space_long)
            pevap_lq = np.linspace(HE.p_he_2r,HE.p_he_2x,size_space_long)
            hevap_st = np.linspace(HE.i_he_2x,HE.i_he_3x,size_space_shrt)
            pevap_st = np.linspace(HE.p_he_2x,HE.p_he_3x,size_space_shrt)
            hevap_vp = np.linspace(HE.i_he_3x,HE.i_he_3, size_space_long)
            pevap_vp = np.linspace(HE.p_he_3x,HE.p_he_3, size_space_long)
            hevap    = np.concatenate((hevap_lq,hevap_st,hevap_vp),axis=None)
            pevap    = np.concatenate((pevap_lq,pevap_st,pevap_vp),axis=None)
        # Recuperator: hot leg
        hreht    = np.linspace(HE.i_he_4r,HE.i_he_4,size_space_long)
        preht    = np.linspace(HE.p_he_4r,HE.p_he_4,size_space_long)
        Treht    = CP.PropsSI('T','P',preht,'H',hreht,HE.parameters['fluid_he'])
        ax1.plot((hreht-i_ref_he)*1e-3,Treht-273.15,'--',color='darkblue',linewidth=3)
        # Condenser
        hcond_lq = np.linspace(HE.i_he_1, HE.i_he_1x,size_space_long)
        pcond_lq = np.linspace(HE.p_he_1, HE.p_he_1x,size_space_long)
        hcond_st = np.linspace(HE.i_he_1x,HE.i_he_4x,size_space_shrt)
        pcond_st = np.linspace(HE.p_he_1x,HE.p_he_4x,size_space_shrt)
        hcond_vp = np.linspace(HE.i_he_4x,HE.i_he_4r,size_space_long)
        pcond_vp = np.linspace(HE.p_he_4x,HE.p_he_4r,size_space_long)
        hcond    = np.concatenate((hcond_lq,hcond_st,hcond_vp),axis=None)
        pcond    = np.concatenate((pcond_lq,pcond_st,pcond_vp),axis=None)
        
    else:
        # Evporator
        if trans_crit:
            hevap    = np.linspace(HE.i_he_2,HE.i_he_3,2*size_space_long)
            pevap    = np.linspace(HE.p_he_2,HE.p_he_3,2*size_space_long)
        else:
            hevap_lq = np.linspace(HE.i_he_2, HE.i_he_2x,size_space_long)
            pevap_lq = np.linspace(HE.p_he_2, HE.p_he_2x,size_space_long)
            hevap_st = np.linspace(HE.i_he_2x,HE.i_he_3x,size_space_shrt)
            pevap_st = np.linspace(HE.p_he_2x,HE.p_he_3x,size_space_shrt)
            hevap_vp = np.linspace(HE.i_he_3x,HE.i_he_3, size_space_long)
            pevap_vp = np.linspace(HE.p_he_3x,HE.p_he_3, size_space_long)
            hevap    = np.concatenate((hevap_lq,hevap_st,hevap_vp),axis=None)
            pevap    = np.concatenate((pevap_lq,pevap_st,pevap_vp),axis=None)
        # Condenser
        hcond_lq = np.linspace(HE.i_he_1, HE.i_he_1x,size_space_long)
        pcond_lq = np.linspace(HE.p_he_1, HE.p_he_1x,size_space_long)
        hcond_st = np.linspace(HE.i_he_1x,HE.i_he_4x,size_space_shrt)
        pcond_st = np.linspace(HE.p_he_1x,HE.p_he_4x,size_space_shrt)
        hcond_vp = np.linspace(HE.i_he_4x,HE.i_he_4, size_space_long)
        pcond_vp = np.linspace(HE.p_he_4x,HE.p_he_4, size_space_long)
        hcond    = np.concatenate((hcond_lq,hcond_st,hcond_vp),axis=None)
        pcond    = np.concatenate((pcond_lq,pcond_st,pcond_vp),axis=None)
    # Condenser
    Tcond    = CP.PropsSI('T','P',pcond,'H',hcond,HE.parameters['fluid_he'])
    ax1.plot((hcond-i_ref_he)*1e-3,Tcond-273.15,'-',color='darkblue',linewidth=3)
    # Pump
    ppump    = np.linspace(HE.p_he_1,HE.p_he_2,size_space_shrt)
    hpump_is = CP.PropsSI('H','P',ppump,'S',HE.s_he_1,HE.parameters['fluid_he'])
    hpump    = HE.i_he_1+(hpump_is-HE.i_he_1)/HE.eta_is_pm
    Tpump    = CP.PropsSI('T','P',ppump,'H',hpump,    HE.parameters['fluid_he'])    
    ax1.plot((hpump-i_ref_he)*1e-3,Tpump-273.15,'-',color='darkblue',linewidth=3)
    # Evporator
    Tevap    = CP.PropsSI('T','P',pevap,'H',hevap,HE.parameters['fluid_he'])
    ax1.plot((hevap-i_ref_he)*1e-3,Tevap-273.15,'-',color='darkblue',linewidth=3)
    # Expander
    pexpa    = np.linspace(HE.p_he_4,HE.p_he_3,size_space_shrt)
    hexpa_is = CP.PropsSI('H','P',pexpa,'S',HE.s_he_3,HE.parameters['fluid_he'])
    hexpa    = HE.i_he_3-HE.eta_is_ex*(HE.i_he_3-hexpa_is)
    Texpa    = CP.PropsSI('T','P',pexpa,'H',hexpa,    HE.parameters['fluid_he'])    
    ax1.plot((hexpa-i_ref_he)*1e-3,Texpa-273.15,'-',color='darkblue',linewidth=3,label='heat engine cycle')
    # Mark the points
    ax1.plot((np.asarray(HE.i_he)-i_ref_he)*1e-3,
              np.asarray(HE.T_he          )-273.15,
              'o',color='darkblue',linewidth=4)
    
    # Secondary Fluids Diagram
    pevap_sf   = np.linspace(HE.p_he_hs_ex,HE.p_he_hs_su,len(hevap))
    hevap_sf   = np.zeros(len(hevap))
    hevap_sf[0]= HE.i_he_hs_ex
    if recup:   fra_hs = (HE.i_he_3-HE.i_he_2r)/(HE.i_he_hs_su-HE.i_he_hs_ex) 
    else:       fra_hs = (HE.i_he_3-HE.i_he_2 )/(HE.i_he_hs_su-HE.i_he_hs_ex) 
    for i in range(len(hevap_sf)-1):
        hevap_sf[i+1] = hevap_sf[i]+(hevap[i+1]-hevap[i])/fra_hs
    Tevap_sf = CP.PropsSI('T','P',pevap_sf,'H',hevap_sf,HE.fluid_he_hs)
    if recup:   hevap_sf = (hevap_sf-i_ref_hs)*(max(hevap)-min(hevap))/(max(hevap_sf)-min(hevap_sf))\
                         + (HE.i_he_2r-HE.i_he_1)
    else:       hevap_sf = (hevap_sf-i_ref_hs)*(max(hevap)-min(hevap))/(max(hevap_sf)-min(hevap_sf))\
                         + (HE.i_he_2 -HE.i_he_1)
    # ax1.plot(hevap_sf*1e-3,Tevap_sf-273.15,'--',color='darkgreen',linewidth=1,label='thermal storage')
    ax1.plot(hevap_sf*1e-3,Tevap_sf-273.15,'--',color='darkgreen',linewidth=1,label='heat source') # !!!
    ax1.plot(np.asarray([hevap_sf[0],hevap_sf[-1]])*1e-3,
              np.asarray([Tevap_sf[0],Tevap_sf[-1]])-273.15,
              'x',color='darkgreen',linewidth=4)
    
    pcond_sf   = np.linspace(HE.p_he_cs_su,HE.p_he_cs_ex,len(hcond))
    hcond_sf   = np.zeros(len(hcond))
    hcond_sf[0]= HE.i_he_cs_su
    if recup:   fra_cs = (HE.i_he_4r-HE.i_he_1)/(HE.i_he_cs_ex-HE.i_he_cs_su) 
    else:       fra_cs = (HE.i_he_4 -HE.i_he_1)/(HE.i_he_cs_ex-HE.i_he_cs_su) 
    for i in range(len(hcond_sf)-1):
        hcond_sf[i+1] = hcond_sf[i]+(hcond[i+1]-hcond[i])/fra_cs
    Tcond_sf = CP.PropsSI('T','P',pcond_sf,'H',hcond_sf,HE.fluid_he_cs)
    hcond_sf = (hcond_sf-i_ref_cs)*(max(hcond)-min(hcond))/(max(hcond_sf)-min(hcond_sf))
    ax1.plot(hcond_sf*1e-3,Tcond_sf-273.15,'--',color='darkslategrey',linewidth=1,label='heat sink')
    ax1.plot(np.asarray([hcond_sf[0],hcond_sf[-1]])*1e-3,
              np.asarray([Tcond_sf[0],Tcond_sf[-1]])-273.15,
              'x',color='darkslategrey',linewidth=4)
    
    # Check if pinch is violated
    pp_hs = np.argmin(Tevap_sf-Tevap)
    pp_cs = np.argmin(Tcond-Tcond_sf)
    ax1.plot(hevap_sf[pp_hs]*1e-3,Tevap[pp_hs]+HE.parameters['dT_he_ev_pp']-273.15, 'x', color='firebrick',     linewidth=4)
    ax1.plot(hcond_sf[pp_cs]*1e-3,Tcond[pp_cs]-HE.parameters['dT_he_cd_pp']-273.15, 'x', color='firebrick',     linewidth=4)
    if round(Tevap_sf[pp_hs]-Tevap[pp_hs],2) < HE.parameters['dT_he_ev_pp']: 
        warnings.warn(f'In the ORC evaporator, the pinch is violated: {Tevap_sf[pp_hs]-Tevap[pp_hs]} < {HE.parameters["dT_he_ev_pp"]}', Warning)
    if round(Tcond[pp_cs]-Tcond_sf[pp_cs],2) < HE.parameters['dT_he_cd_pp']:
        warnings.warn(f'In the ORC condenser, the pinch is violated: {Tcond[pp_cs]-Tcond_sf[pp_cs]} < {HE.parameters["dT_he_cd_pp"]}', Warning)

    # Axes
    ax1.legend(loc=(0.05,0.75),frameon=False,labelcolor='linecolor')
    ax1.text(0.055,0.72,'$\mathrm{\eta_{ORC}}$'+f' = {round(HE.eta_he_cyclen*100,2)} [%]',transform=plt.gca().transAxes)
    ax1.set_xlabel('enthalpy [kJ/kg]',fontsize=14, color='dimgray')
    ax1.spines["bottom"].set_bounds(hcond_sf[0]*1e-3,hevap_sf[-1]*1e-3)
    ax1.set_xticks(np.round([hcond_sf[0]*1e-3,hevap_sf[0]*1e-3,hcond_sf[-1]*1e-3,hevap_sf[-1]*1e-3],1))
    ax1.spines["bottom"].set_linewidth(1.5)
    ax1.spines["bottom"].set_color('dimgray')
    ax1.set_ylabel('temperature [°C]',fontsize=14, color='dimgray')
    ax1.spines["left"].set_bounds(Tcond_sf[0]-273.15,max(Tevap_sf[-1],max(HE.T_he))-273.15)
    ax1.set_yticks(np.round([min(HE.T_he)-273.15,HE.T_he_1x-273.15,
                    max(HE.T_he)-273.15,Tevap_sf[0]-273.15,Tevap_sf[-1]-273.15,
                    Tcond_sf[-1]-273.15,Tcond_sf[0]-273.15],1))
    ax1.spines["left"].set_linewidth(1.5)
    ax1.spines["left"].set_color('dimgray')
    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)
    ax1.tick_params(axis='both', labelsize=12, colors='dimgray', width=1.5, length=8, direction='in')
    ax1.label_outer()
    fig.tight_layout()
    plt.show()
    out = fig
    return out

def plot_Ts_ORC(HE):
    
    # Process inputs
    HE.export_states()
    i_ref_he = HE.i_he_1
    i_ref_cs = HE.i_he_cs_su
    i_ref_hs = HE.i_he_hs_ex
    
    cell_res_long   = 0.05 # K
    size_space_long = int(1/cell_res_long*(HE.T_he_hs_su-HE.T_he_hs_ex)+1)
    cell_res_shrt   = 0.10 # K
    size_space_shrt = int(1/cell_res_shrt*(HE.T_he_cs_ex-HE.T_he_cs_su)+1)
    
    try:    recup, i_he_2r, i_he_4r = True, HE.i_he_2r, HE.i_he_4r
    except: recup                   = False
    if HE.p_he_2 >= CP.PropsSI('P_CRITICAL',HE.parameters['fluid_he']): trans_crit = True
    else:                                                               trans_crit = False
    
    # Saturation Curve
    T_min = max(min(HE.T_he)-2,CP.PropsSI('T_MIN',HE.parameters['fluid_he'])+1)
    # T_max = min(HE.T_he_hs_su, CP.PropsSI('T_MAX',HE.parameters['fluid_he'])-1,CP.PropsSI('T_CRITICAL',HE.parameters['fluid_he'])-1)
    T_max = min(max(HE.T_he)+0,CP.PropsSI('T_MAX',HE.parameters['fluid_he'])-1,CP.PropsSI('T_CRITICAL',HE.parameters['fluid_he'])-1)
    temp = np.linspace(T_min,T_max,1000)
    s_lq = np.zeros(len(temp))
    s_vp = np.zeros(len(temp))
    
    my_state_he = AbstractState("HEOS",HE.parameters['fluid_he'])
    for i, t in enumerate(temp):
        my_state_he.update(CoolProp.QT_INPUTS,1,t)
        s_vp[i] = my_state_he.smass()
        my_state_he.update(CoolProp.QT_INPUTS,0,t)
        s_lq[i] = my_state_he.smass()
        
    fig = plt.figure()
    ax1 = fig.subplots(1)
    ax1.plot(s_lq*1e-3,temp-273.15,'-',color='black',linewidth=1,label=HE.parameters['fluid_he'])
    ax1.plot(s_vp*1e-3,temp-273.15,'-',color='black',linewidth=1)
    ax1.plot(5e-4*(s_lq[len(temp)-1]+s_vp[len(temp)-1]),T_max-273.15,'.',color='black',linewidth=1)
    
    # Heat Engine Diagram
    if recup:
        # Recuperator: cold leg
        hrelt    = np.linspace(HE.i_he_2,HE.i_he_2r,size_space_long)
        prelt    = np.linspace(HE.p_he_2,HE.p_he_2r,size_space_long)
        Trelt    = CP.PropsSI('T','P',prelt,'H',hrelt,HE.parameters['fluid_he'])
        srelt    = CP.PropsSI('S','P',prelt,'H',hrelt,HE.parameters['fluid_he'])
        ax1.plot(srelt*1e-3,Trelt-273.15,'--',color='darkblue',linewidth=3)
        # Evporator        
        if trans_crit:
            hevap    = np.linspace(HE.i_he_2r,HE.i_he_3,2*size_space_long)
            pevap    = np.linspace(HE.p_he_2r,HE.p_he_3,2*size_space_long)
        else:
            hevap_lq = np.linspace(HE.i_he_2r,HE.i_he_2x,size_space_long)
            pevap_lq = np.linspace(HE.p_he_2r,HE.p_he_2x,size_space_long)
            hevap_st = np.linspace(HE.i_he_2x,HE.i_he_3x,size_space_shrt)
            pevap_st = np.linspace(HE.p_he_2x,HE.p_he_3x,size_space_shrt)
            hevap_vp = np.linspace(HE.i_he_3x,HE.i_he_3, size_space_long)
            pevap_vp = np.linspace(HE.p_he_3x,HE.p_he_3, size_space_long)
            hevap    = np.concatenate((hevap_lq,hevap_st,hevap_vp),axis=None)
            pevap    = np.concatenate((pevap_lq,pevap_st,pevap_vp),axis=None)
        # Recuperator: hot leg
        hreht    = np.linspace(HE.i_he_4r,HE.i_he_4,size_space_long)
        preht    = np.linspace(HE.p_he_4r,HE.p_he_4,size_space_long)
        Treht    = CP.PropsSI('T','P',preht,'H',hreht,HE.parameters['fluid_he'])
        sreht    = CP.PropsSI('S','P',preht,'H',hreht,HE.parameters['fluid_he'])
        ax1.plot(sreht*1e-3,Treht-273.15,'--',color='darkblue',linewidth=3)
        # Condenser
        hcond_lq = np.linspace(HE.i_he_1, HE.i_he_1x,size_space_long)
        pcond_lq = np.linspace(HE.p_he_1, HE.p_he_1x,size_space_long)
        hcond_st = np.linspace(HE.i_he_1x,HE.i_he_4x,size_space_shrt)
        pcond_st = np.linspace(HE.p_he_1x,HE.p_he_4x,size_space_shrt)
        hcond_vp = np.linspace(HE.i_he_4x,HE.i_he_4r,size_space_long)
        pcond_vp = np.linspace(HE.p_he_4x,HE.p_he_4r,size_space_long)
        hcond    = np.concatenate((hcond_lq,hcond_st,hcond_vp),axis=None)
        pcond    = np.concatenate((pcond_lq,pcond_st,pcond_vp),axis=None)
    else:
        # Evporator
        if trans_crit:
            hevap    = np.linspace(HE.i_he_2,HE.i_he_3,2*size_space_long)
            pevap    = np.linspace(HE.p_he_2,HE.p_he_3,2*size_space_long)
        else:
            hevap_lq = np.linspace(HE.i_he_2, HE.i_he_2x,size_space_long)
            pevap_lq = np.linspace(HE.p_he_2, HE.p_he_2x,size_space_long)
            hevap_st = np.linspace(HE.i_he_2x,HE.i_he_3x,size_space_shrt)
            pevap_st = np.linspace(HE.p_he_2x,HE.p_he_3x,size_space_shrt)
            hevap_vp = np.linspace(HE.i_he_3x,HE.i_he_3, size_space_long)
            pevap_vp = np.linspace(HE.p_he_3x,HE.p_he_3, size_space_long)
            hevap    = np.concatenate((hevap_lq,hevap_st,hevap_vp),axis=None)
            pevap    = np.concatenate((pevap_lq,pevap_st,pevap_vp),axis=None)
        # Condenser
        hcond_lq = np.linspace(HE.i_he_1, HE.i_he_1x,size_space_long)
        pcond_lq = np.linspace(HE.p_he_1, HE.p_he_1x,size_space_long)
        hcond_st = np.linspace(HE.i_he_1x,HE.i_he_4x,size_space_shrt)
        pcond_st = np.linspace(HE.p_he_1x,HE.p_he_4x,size_space_shrt)
        hcond_vp = np.linspace(HE.i_he_4x,HE.i_he_4, size_space_long)
        pcond_vp = np.linspace(HE.p_he_4x,HE.p_he_4, size_space_long)
        hcond    = np.concatenate((hcond_lq,hcond_st,hcond_vp),axis=None)
        pcond    = np.concatenate((pcond_lq,pcond_st,pcond_vp),axis=None)
    # Condenser
    Tcond    = CP.PropsSI('T','P',pcond,'H',hcond,HE.parameters['fluid_he'])
    scond    = CP.PropsSI('S','P',pcond,'H',hcond,HE.parameters['fluid_he'])
    ax1.plot(scond*1e-3,Tcond-273.15,'-',color='darkblue',linewidth=3)
    # Pump
    ppump    = np.linspace(HE.p_he_1,HE.p_he_2,size_space_shrt)
    hpump_is = CP.PropsSI('H','P',ppump,'S',HE.s_he_1,HE.parameters['fluid_he'])
    hpump    = HE.i_he_1+(hpump_is-HE.i_he_1)/HE.eta_is_pm
    Tpump    = CP.PropsSI('T','P',ppump,'H',hpump,    HE.parameters['fluid_he'])    
    spump    = CP.PropsSI('S','P',ppump,'H',hpump,    HE.parameters['fluid_he'])    
    ax1.plot(spump*1e-3,Tpump-273.15,'-',color='darkblue',linewidth=3)
    # Evporator
    Tevap    = CP.PropsSI('T','P',pevap,'H',hevap,HE.parameters['fluid_he'])
    sevap    = CP.PropsSI('S','P',pevap,'H',hevap,HE.parameters['fluid_he'])
    ax1.plot(sevap*1e-3,Tevap-273.15,'-',color='darkblue',linewidth=3)
    # Expander
    pexpa    = np.linspace(HE.p_he_4,HE.p_he_3,size_space_shrt)
    hexpa_is = CP.PropsSI('H','P',pexpa,'S',HE.s_he_3,HE.parameters['fluid_he'])
    hexpa    = HE.i_he_3-HE.eta_is_ex*(HE.i_he_3-hexpa_is)
    Texpa    = CP.PropsSI('T','P',pexpa,'H',hexpa,    HE.parameters['fluid_he'])    
    sexpa    = CP.PropsSI('S','P',pexpa,'H',hexpa,    HE.parameters['fluid_he'])    
    ax1.plot(sexpa*1e-3,Texpa-273.15,'-',color='darkblue',linewidth=3,label='heat engine cycle')
    # Mark the points
    ax1.plot(np.asarray(HE.s_he)*1e-3,
             np.asarray(HE.T_he)-273.15,
             'o',color='darkblue',linewidth=4)
    
    # Axes
    ax1.legend(loc=(0.05,0.855),frameon=False,labelcolor='linecolor')
    ax1.text(0.055,0.82,'$\mathrm{\eta_{ORC}}$'+f' = {round(HE.eta_he_cyclen*100,2)} [%]',transform=plt.gca().transAxes)
    ax1.set_xlabel('entropy [kJ/kg/K]',fontsize=14, color='dimgray')
    ax1.spines["bottom"].set_bounds(min(HE.s_he)*1e-3,max(HE.s_he)*1e-3)
    
    ax1.set_xticks(np.round([min(HE.s_he)*1e-3,
                             max(HE.s_he)*1e-3],3))
    
    ax1.spines["bottom"].set_linewidth(1.5)
    ax1.spines["bottom"].set_color('dimgray')
    ax1.set_ylabel('temperature [°C]',fontsize=14, color='dimgray')
    ax1.spines["left"].set_bounds(min(HE.T_he)-273.15,max(HE.T_he)-273.15)
    if trans_crit:
        if recup:
            ax1.set_yticks(np.round([min(HE.T_he)-273.15,HE.T_he_1x-273.15,
                                     HE.T_he_2r-273.15,HE.T_he_4-273.15,
                                     HE.T_he_4r-273.15,
                                     max(HE.T_he)-273.15],1))
        else:
            ax1.set_yticks(np.round([min(HE.T_he)-273.15,HE.T_he_1x-273.15,
                                     HE.T_he_4-273.15,
                                     max(HE.T_he)-273.15],1))
    else:
        if recup:
            ax1.set_yticks(np.round([min(HE.T_he)-273.15,HE.T_he_1x-273.15,
                                     HE.T_he_2x-273.15,HE.T_he_4-273.15,
                                     HE.T_he_2r-273.15,HE.T_he_4r-273.15,
                                     max(HE.T_he)-273.15],1))
        else:
            ax1.set_yticks(np.round([min(HE.T_he)-273.15,HE.T_he_1x-273.15,
                                     HE.T_he_2x-273.15,HE.T_he_4-273.15,
                                     max(HE.T_he)-273.15],1))
    
    ax1.spines["left"].set_linewidth(1.5)
    ax1.spines["left"].set_color('dimgray')
    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)
    ax1.tick_params(axis='both', labelsize=12, colors='dimgray', width=1.5, length=8, direction='in')
    ax1.label_outer()
    fig.tight_layout()
    plt.show()
    out = fig
    return out

#%% CB

def plot_Ts_VCHP_ORC_SHTES_norm(CB):
    
    fig = plt.figure()
    ax1 = fig.subplots(1)
    
    # --- Process inputs ------------------------------------------------------
    CB.export_states()  
    cell_res_long   = 0.05 # K
    cell_res_shrt   = 0.10 # K
    size_space_long_hp = int(1/cell_res_long*(CB.my_HP.T_hp_hs_ex-CB.my_HP.T_hp_hs_su)+1)
    size_space_shrt_hp = int(1/cell_res_shrt*(CB.my_HP.T_hp_cs_su-CB.my_HP.T_hp_cs_ex)+1)
    size_space_long_he = int(1/cell_res_long*(CB.my_HE.T_he_hs_su-CB.my_HE.T_he_hs_ex)+1)
    size_space_shrt_he = int(1/cell_res_shrt*(CB.my_HE.T_he_cs_ex-CB.my_HE.T_he_cs_su)+1)
    # -> Recuperated?
    try:    recup_hp, i_hp_1r, i_hp_3r = True, CB.my_HP.i_hp_1r, CB.my_HP.i_hp_3r
    except: recup_hp                   = False
    try:    recup_he, i_he_2r, i_he_4r = True, CB.my_HE.i_he_2r, CB.my_HE.i_he_4r
    except: recup_he                   = False
    # -> Transcritical?
    if CB.my_HE.p_he_2 >= CP.PropsSI('P_CRITICAL',CB.my_HE.parameters['fluid_he']): trans_crit_he = True
    else:                                                                           trans_crit_he = False
    if CB.my_HP.p_hp_2 >= CP.PropsSI('P_CRITICAL',CB.my_HP.parameters['fluid_hp']): trans_crit_hp = True
    else:                                                                           trans_crit_hp = False
    # -> Saturation curve for heat pump
    my_state_hp = AbstractState("HEOS",CB.parameters['fluid_hp'])
    T_min_hp = max(min(CB.my_HP.T_hp)-5,CP.PropsSI('T_MIN',CB.my_HP.parameters['fluid_hp'])+1)
    T_max_hp = min(max(CB.my_HP.T_hp)+5,CP.PropsSI('T_MAX',CB.my_HP.parameters['fluid_hp'])-1,CP.PropsSI('T_CRITICAL',CB.my_HP.parameters['fluid_hp'])-0.5)
    temp_hpf = np.linspace(T_max_hp-1,T_max_hp,1000)
    if recup_hp: temp_hpl = np.linspace(CB.my_HP.T_hp_3r-10,T_max_hp,1000)
    else:        temp_hpl = np.linspace(CB.my_HP.T_hp_3 -10,T_max_hp,1000)
    iil      = np.searchsorted(temp_hpl, temp_hpf)
    temp_hpl = np.insert(temp_hpl, iil, temp_hpf)
    s_lq_hp  = np.zeros(len(temp_hpl))
    for i, t in enumerate(temp_hpl):
        my_state_hp.update(CoolProp.QT_INPUTS,0,t)
        s_lq_hp[i] = my_state_hp.smass()
    temp_hpv = np.linspace(T_min_hp,  T_max_hp,1000)
    iiv      = np.searchsorted(temp_hpv, temp_hpf)
    temp_hpv = np.insert(temp_hpv, iiv, temp_hpf)
    s_vp_hp  = np.zeros(len(temp_hpv))
    for i, t in enumerate(temp_hpv):
        my_state_hp.update(CoolProp.QT_INPUTS,1,t)
        s_vp_hp[i] = my_state_hp.smass()
    # -> Saturation curve for heat engine
    T_min_he = max(min(CB.my_HE.T_he)-5,CP.PropsSI('T_MIN',CB.my_HE.parameters['fluid_he'])+1)
    T_max_he = min(max(CB.my_HE.T_he)+5,CP.PropsSI('T_MAX',CB.my_HE.parameters['fluid_he'])-1,CP.PropsSI('T_CRITICAL',CB.my_HE.parameters['fluid_he'])-0.5)
    temp_he  = np.linspace(T_min_he,  T_max_he,1000)
    temp_hef = np.linspace(T_max_he-1,T_max_he,1000)
    ii       = np.searchsorted(temp_he, temp_hef)
    temp_he  = np.insert(temp_he, ii, temp_hef)
    s_lq_he  = np.zeros(len(temp_he))
    s_vp_he  = np.zeros(len(temp_he))
    my_state_he = AbstractState("HEOS",CB.parameters['fluid_he'])
    for i, t in enumerate(temp_he):
        my_state_he.update(CoolProp.QT_INPUTS,1,t)
        s_vp_he[i] = my_state_he.smass()
        my_state_he.update(CoolProp.QT_INPUTS,0,t)
        s_lq_he[i] = my_state_he.smass()
    # -> Plot surves
    s_sat_hp_min = CB.my_HP.s_hp_3#min(s_lq_hp)
    s_sat_hp_max = CB.my_HP.s_hp_2#max(s_vp_hp)
    ds_sat_hp    = s_sat_hp_max-s_sat_hp_min
    if recup_he:    s_sat_he_min = CB.my_HE.s_he_2r#min(s_lq_he)
    else:           s_sat_he_min = CB.my_HE.s_he_2#min(s_lq_he)
    s_sat_he_max = CB.my_HE.s_he_3#max(s_vp_he)
    ds_sat_he    = s_sat_he_max-s_sat_he_min
    ax1.plot((s_lq_hp-s_sat_hp_min)/ds_sat_hp,temp_hpl-273.15,'--',color='black',linewidth=0.8, label=CB.parameters['fluid_hp'])
    ax1.plot((s_vp_hp-s_sat_hp_min)/ds_sat_hp,temp_hpv-273.15,'--',color='black',linewidth=0.8)
    if CB.parameters['fluid_he'] != CB.parameters['fluid_hp']:
        ax1.plot((s_lq_he-s_sat_he_min)/ds_sat_he,temp_he-273.15,'-.',color='black',linewidth=0.8, label=CB.parameters['fluid_he'])
        ax1.plot((s_vp_he-s_sat_he_min)/ds_sat_he,temp_he-273.15,'-.',color='black',linewidth=0.8)
    
    # --- Heat Pump Diagram ---------------------------------------------------
    if recup_hp:
        # -> Recuperator: cold leg
        hrelt_hp    = np.linspace(CB.my_HP.i_hp_1,CB.my_HP.i_hp_1r,size_space_long_hp)
        prelt_hp    = np.linspace(CB.my_HP.p_hp_1,CB.my_HP.p_hp_1r,size_space_long_hp)
        Trelt_hp    = CP.PropsSI('T','P',prelt_hp,'H',hrelt_hp,CB.my_HP.parameters['fluid_hp'])
        srelt_hp    = CP.PropsSI('S','P',prelt_hp,'H',hrelt_hp,CB.my_HP.parameters['fluid_hp'])
        ax1.plot((srelt_hp-s_sat_hp_min)/ds_sat_hp,Trelt_hp-273.15,'--',color='firebrick',linewidth=3)
        # -> Compression
        pcomp_hp    = np.linspace(CB.my_HP.p_hp_1r,CB.my_HP.p_hp_2,size_space_shrt_hp)
        hcomp_hp_is = CP.PropsSI('H','P',pcomp_hp,'S',CB.my_HP.s_hp_1r,CB.my_HP.parameters['fluid_hp'])
        hcomp_hp    = CB.my_HP.i_hp_1r+(hcomp_hp_is-CB.my_HP.i_hp_1r)/CB.my_HP.eta_is_cp
        # -> Recuperator: cold leg
        hreht_hp    = np.linspace(CB.my_HP.i_hp_3,CB.my_HP.i_hp_3r,size_space_long_hp)
        preht_hp    = np.linspace(CB.my_HP.p_hp_3,CB.my_HP.p_hp_3r,size_space_long_hp)
        Treht_hp    = CP.PropsSI('T','P',preht_hp,'H',hreht_hp,CB.my_HP.parameters['fluid_hp'])
        sreht_hp    = CP.PropsSI('S','P',preht_hp,'H',hreht_hp,CB.my_HP.parameters['fluid_hp'])
        ax1.plot((sreht_hp-s_sat_hp_min)/ds_sat_hp,Treht_hp-273.15,'--',color='firebrick',linewidth=3)
        # -> Expansion
        hexpv_hp    = np.linspace(CB.my_HP.i_hp_3r,CB.my_HP.i_hp_4,size_space_shrt_hp)
        pexpv_hp    = np.linspace(CB.my_HP.p_hp_3r,CB.my_HP.p_hp_4,size_space_shrt_hp)
    else:
        # -> Compression
        pcomp_hp    = np.linspace(CB.my_HP.p_hp_1,CB.my_HP.p_hp_2,size_space_shrt_hp)
        hcomp_hp_is = CP.PropsSI('H','P',pcomp_hp,'S',CB.my_HP.s_hp_1,CB.my_HP.parameters['fluid_hp'])
        hcomp_hp    = CB.my_HP.i_hp_1+(hcomp_hp_is-CB.my_HP.i_hp_1)/CB.my_HP.eta_is_cp
        # -> Expansion
        hexpv_hp    = np.linspace(CB.my_HP.i_hp_3,CB.my_HP.i_hp_4,size_space_shrt_hp)
        pexpv_hp    = np.linspace(CB.my_HP.p_hp_3,CB.my_HP.p_hp_4,size_space_shrt_hp)
    # -> Compression
    Tcomp_hp        = CP.PropsSI('T','P',pcomp_hp,'H',hcomp_hp,    CB.my_HP.parameters['fluid_hp'])
    scomp_hp        = CP.PropsSI('S','P',pcomp_hp,'H',hcomp_hp,    CB.my_HP.parameters['fluid_hp'])
    ax1.plot((scomp_hp-s_sat_hp_min)/ds_sat_hp,Tcomp_hp-273.15,'-',color='firebrick',linewidth=3)
    # -> Condensation
    if trans_crit_hp:
        hcond_hp    = np.linspace(CB.my_HP.i_hp_3,CB.my_HP.i_hp_2,2*size_space_long_hp)
        pcond_hp    = np.linspace(CB.my_HP.p_hp_3,CB.my_HP.p_hp_2,2*size_space_long_hp)
    else:
        hcond_hp_lq = np.linspace(CB.my_HP.i_hp_3, CB.my_HP.i_hp_3x,size_space_long_hp)
        pcond_hp_lq = np.linspace(CB.my_HP.p_hp_3, CB.my_HP.p_hp_3x,size_space_long_hp)
        hcond_hp_st = np.linspace(CB.my_HP.i_hp_3x,CB.my_HP.i_hp_2x,size_space_shrt_hp)
        pcond_hp_st = np.linspace(CB.my_HP.p_hp_3x,CB.my_HP.p_hp_2x,size_space_shrt_hp)
        hcond_hp_vp = np.linspace(CB.my_HP.i_hp_2x,CB.my_HP.i_hp_2, size_space_long_hp)
        pcond_hp_vp = np.linspace(CB.my_HP.p_hp_2x,CB.my_HP.p_hp_2, size_space_long_hp)
        hcond_hp    = np.concatenate((hcond_hp_lq,hcond_hp_st,hcond_hp_vp),axis=None)
        pcond_hp    = np.concatenate((pcond_hp_lq,pcond_hp_st,pcond_hp_vp),axis=None)
    Tcond_hp    = CP.PropsSI('T','P',pcond_hp,'H',hcond_hp,    CB.my_HP.parameters['fluid_hp'])
    scond_hp    = CP.PropsSI('S','P',pcond_hp,'H',hcond_hp,    CB.my_HP.parameters['fluid_hp'])
    ax1.plot((scond_hp-s_sat_hp_min)/ds_sat_hp,Tcond_hp-273.15,'-',color='firebrick',linewidth=3)
    # -> Expansion
    Texpv_hp    = CP.PropsSI('T','P',pexpv_hp,'H',hexpv_hp,    CB.my_HP.parameters['fluid_hp'])
    sexpv_hp    = CP.PropsSI('S','P',pexpv_hp,'H',hexpv_hp,    CB.my_HP.parameters['fluid_hp'])
    ax1.plot((sexpv_hp-s_sat_hp_min)/ds_sat_hp,Texpv_hp-273.15,'-',color='firebrick',linewidth=3)
    # -> Evaporation
    hevap_hp_lq = np.linspace(CB.my_HP.i_hp_4, CB.my_HP.i_hp_4x,size_space_long_hp)
    pevap_hp_lq = np.linspace(CB.my_HP.p_hp_4, CB.my_HP.p_hp_4x,size_space_long_hp)
    hevap_hp_st = np.linspace(CB.my_HP.i_hp_4x,CB.my_HP.i_hp_1x,size_space_shrt_hp)
    pevap_hp_st = np.linspace(CB.my_HP.p_hp_4x,CB.my_HP.p_hp_1x,size_space_shrt_hp)
    hevap_hp_vp = np.linspace(CB.my_HP.i_hp_1x,CB.my_HP.i_hp_1, size_space_long_hp)
    pevap_hp_vp = np.linspace(CB.my_HP.p_hp_1x,CB.my_HP.p_hp_1, size_space_long_hp)
    hevap_hp    = np.concatenate((hevap_hp_lq,hevap_hp_st,hevap_hp_vp),axis=None)
    pevap_hp    = np.concatenate((pevap_hp_lq,pevap_hp_st,pevap_hp_vp),axis=None)
    Tevap_hp    = CP.PropsSI('T','P',pevap_hp,'H',hevap_hp,    CB.my_HP.parameters['fluid_hp'])
    sevap_hp    = CP.PropsSI('S','P',pevap_hp,'H',hevap_hp,    CB.my_HP.parameters['fluid_hp'])
    ax1.plot((sevap_hp-s_sat_hp_min)/ds_sat_hp,Tevap_hp-273.15,'-',color='firebrick',linewidth=3)
    # -> Mark the points
    ax1.plot((np.asarray(CB.my_HP.s_hp)-s_sat_hp_min)/ds_sat_hp,
              np.asarray(CB.my_HP.T_hp)-273.15,
              'o',color='firebrick',linewidth=4)
    
    # --- Heat Engine Diagram -------------------------------------------------
    if recup_he:
        # -> Recuperator: cold leg
        hrelt_he    = np.linspace(CB.my_HE.i_he_2,CB.my_HE.i_he_2r,size_space_long_he)
        prelt_he    = np.linspace(CB.my_HE.p_he_2,CB.my_HE.p_he_2r,size_space_long_he)
        Trelt_he    = CP.PropsSI('T','P',prelt_he,'H',hrelt_he,CB.my_HE.parameters['fluid_he'])
        srelt_he    = CP.PropsSI('S','P',prelt_he,'H',hrelt_he,CB.my_HE.parameters['fluid_he'])
        ax1.plot((srelt_he-s_sat_he_min)/ds_sat_he,Trelt_he-273.15,'--',color='darkblue',linewidth=3)
        # -> Evporator        
        if trans_crit_he:
            hevap_he    = np.linspace(CB.my_HE.i_he_2r,CB.my_HE.i_he_3,2*size_space_long_he)
            pevap_he    = np.linspace(CB.my_HE.p_he_2r,CB.my_HE.p_he_3,2*size_space_long_he)
        else:
            hevap_he_lq = np.linspace(CB.my_HE.i_he_2r,CB.my_HE.i_he_2x,size_space_long_he)
            pevap_he_lq = np.linspace(CB.my_HE.p_he_2r,CB.my_HE.p_he_2x,size_space_long_he)
            hevap_he_st = np.linspace(CB.my_HE.i_he_2x,CB.my_HE.i_he_3x,size_space_shrt_he)
            pevap_he_st = np.linspace(CB.my_HE.p_he_2x,CB.my_HE.p_he_3x,size_space_shrt_he)
            hevap_he_vp = np.linspace(CB.my_HE.i_he_3x,CB.my_HE.i_he_3, size_space_long_he)
            pevap_he_vp = np.linspace(CB.my_HE.p_he_3x,CB.my_HE.p_he_3, size_space_long_he)
            hevap_he    = np.concatenate((hevap_he_lq,hevap_he_st,hevap_he_vp),axis=None)
            pevap_he    = np.concatenate((pevap_he_lq,pevap_he_st,pevap_he_vp),axis=None)
        # -> Recuperator: hot leg
        hreht_he    = np.linspace(CB.my_HE.i_he_4r,CB.my_HE.i_he_4,size_space_long_he)
        preht_he    = np.linspace(CB.my_HE.p_he_4r,CB.my_HE.p_he_4,size_space_long_he)
        Treht_he    = CP.PropsSI('T','P',preht_he,'H',hreht_he,CB.my_HE.parameters['fluid_he'])
        sreht_he    = CP.PropsSI('S','P',preht_he,'H',hreht_he,CB.my_HE.parameters['fluid_he'])
        ax1.plot((sreht_he-s_sat_he_min)/ds_sat_he,Treht_he-273.15,'--',color='darkblue',linewidth=3)
        # -> Condenser
        hcond_he_lq = np.linspace(CB.my_HE.i_he_1, CB.my_HE.i_he_1x,size_space_long_he)
        pcond_he_lq = np.linspace(CB.my_HE.p_he_1, CB.my_HE.p_he_1x,size_space_long_he)
        hcond_he_st = np.linspace(CB.my_HE.i_he_1x,CB.my_HE.i_he_4x,size_space_shrt_he)
        pcond_he_st = np.linspace(CB.my_HE.p_he_1x,CB.my_HE.p_he_4x,size_space_shrt_he)
        hcond_he_vp = np.linspace(CB.my_HE.i_he_4x,CB.my_HE.i_he_4r,size_space_long_he)
        pcond_he_vp = np.linspace(CB.my_HE.p_he_4x,CB.my_HE.p_he_4r,size_space_long_he)
        hcond_he    = np.concatenate((hcond_he_lq,hcond_he_st,hcond_he_vp),axis=None)
        pcond_he    = np.concatenate((pcond_he_lq,pcond_he_st,pcond_he_vp),axis=None)
    else:
        # -> Evporator
        if trans_crit_he:
            hevap_he    = np.linspace(CB.my_HE.i_he_2,CB.my_HE.i_he_3,2*size_space_long_he)
            pevap_he    = np.linspace(CB.my_HE.p_he_2,CB.my_HE.p_he_3,2*size_space_long_he)
        else:
            hevap_he_lq = np.linspace(CB.my_HE.i_he_2, CB.my_HE.i_he_2x,size_space_long_he)
            pevap_he_lq = np.linspace(CB.my_HE.p_he_2, CB.my_HE.p_he_2x,size_space_long_he)
            hevap_he_st = np.linspace(CB.my_HE.i_he_2x,CB.my_HE.i_he_3x,size_space_shrt_he)
            pevap_he_st = np.linspace(CB.my_HE.p_he_2x,CB.my_HE.p_he_3x,size_space_shrt_he)
            hevap_he_vp = np.linspace(CB.my_HE.i_he_3x,CB.my_HE.i_he_3, size_space_long_he)
            pevap_he_vp = np.linspace(CB.my_HE.p_he_3x,CB.my_HE.p_he_3, size_space_long_he)
            hevap_he    = np.concatenate((hevap_he_lq,hevap_he_st,hevap_he_vp),axis=None)
            pevap_he    = np.concatenate((pevap_he_lq,pevap_he_st,pevap_he_vp),axis=None)
        # Condenser
        hcond_he_lq = np.linspace(CB.my_HE.i_he_1, CB.my_HE.i_he_1x,size_space_long_he)
        pcond_he_lq = np.linspace(CB.my_HE.p_he_1, CB.my_HE.p_he_1x,size_space_long_he)
        hcond_he_st = np.linspace(CB.my_HE.i_he_1x,CB.my_HE.i_he_4x,size_space_shrt_he)
        pcond_he_st = np.linspace(CB.my_HE.p_he_1x,CB.my_HE.p_he_4x,size_space_shrt_he)
        hcond_he_vp = np.linspace(CB.my_HE.i_he_4x,CB.my_HE.i_he_4, size_space_long_he)
        pcond_he_vp = np.linspace(CB.my_HE.p_he_4x,CB.my_HE.p_he_4, size_space_long_he)
        hcond_he    = np.concatenate((hcond_he_lq,hcond_he_st,hcond_he_vp),axis=None)
        pcond_he    = np.concatenate((pcond_he_lq,pcond_he_st,pcond_he_vp),axis=None)
    # Condenser
    Tcond_he    = CP.PropsSI('T','P',pcond_he,'H',hcond_he,CB.my_HE.parameters['fluid_he'])
    scond_he    = CP.PropsSI('S','P',pcond_he,'H',hcond_he,CB.my_HE.parameters['fluid_he'])
    ax1.plot((scond_he-s_sat_he_min)/ds_sat_he,Tcond_he-273.15,'-',color='darkblue',linewidth=3)
    # Pump
    ppump_he    = np.linspace(CB.my_HE.p_he_1,CB.my_HE.p_he_2,size_space_shrt_he)
    hpump_he_is = CP.PropsSI('H','P',ppump_he,'S',CB.my_HE.s_he_1,CB.my_HE.parameters['fluid_he'])
    hpump_he    = CB.my_HE.i_he_1+(hpump_he_is-CB.my_HE.i_he_1)/CB.my_HE.eta_is_pm
    Tpump_he    = CP.PropsSI('T','P',ppump_he,'H',hpump_he,    CB.my_HE.parameters['fluid_he'])    
    spump_he    = CP.PropsSI('S','P',ppump_he,'H',hpump_he,    CB.my_HE.parameters['fluid_he'])    
    ax1.plot((spump_he-s_sat_he_min)/ds_sat_he,Tpump_he-273.15,'-',color='darkblue',linewidth=3)
    # Evporator
    Tevap_he    = CP.PropsSI('T','P',pevap_he,'H',hevap_he,CB.my_HE.parameters['fluid_he'])
    sevap_he    = CP.PropsSI('S','P',pevap_he,'H',hevap_he,CB.my_HE.parameters['fluid_he'])
    ax1.plot((sevap_he-s_sat_he_min)/ds_sat_he,Tevap_he-273.15,'-',color='darkblue',linewidth=3)
    # Expander
    pexpa_he    = np.linspace(CB.my_HE.p_he_4,CB.my_HE.p_he_3,size_space_shrt_he)
    hexpa_he_is = CP.PropsSI('H','P',pexpa_he,'S',CB.my_HE.s_he_3,CB.my_HE.parameters['fluid_he'])
    hexpa_he    = CB.my_HE.i_he_3-CB.my_HE.eta_is_ex*(CB.my_HE.i_he_3-hexpa_he_is)
    Texpa_he    = CP.PropsSI('T','P',pexpa_he,'H',hexpa_he,    CB.my_HE.parameters['fluid_he'])    
    sexpa_he    = CP.PropsSI('S','P',pexpa_he,'H',hexpa_he,    CB.my_HE.parameters['fluid_he'])    
    ax1.plot((sexpa_he-s_sat_he_min)/ds_sat_he,Texpa_he-273.15,'-',color='darkblue',linewidth=3)
    # Mark the points
    ax1.plot((np.asarray(CB.my_HE.s_he)-s_sat_he_min)/ds_sat_he,
              np.asarray(CB.my_HE.T_he)-273.15,
             'o',color='darkblue',linewidth=4)
    
    # --- Seconfary Fluids Diagram --------------------------------------------
    ax1.plot((np.asarray([CB.my_HP.s_hp_2,CB.my_HP.s_hp_3])-s_sat_hp_min)/ds_sat_hp,np.asarray([CB.T_st_ht,   CB.T_st_lt   ])-273.15,'--',color='darkgreen',    zorder=-5,linewidth=1.2, label='storage')
    ax1.plot((np.asarray([CB.my_HP.s_hp_2,CB.my_HP.s_hp_3])-s_sat_hp_min)/ds_sat_hp,np.asarray([CB.T_st_ht,   CB.T_st_lt   ])-273.15,'x', color='darkgreen',    zorder=-5,linewidth=1.2)
    ax1.plot((np.asarray([CB.my_HP.s_hp_4,CB.my_HP.s_hp_1])-s_sat_hp_min)/ds_sat_hp,np.asarray([CB.T_hp_cs_ex,CB.T_hp_cs_su])-273.15,'--',color='darkslategrey',zorder=-5,linewidth=1.2)
    ax1.plot((np.asarray([CB.my_HP.s_hp_4,CB.my_HP.s_hp_1])-s_sat_hp_min)/ds_sat_hp,np.asarray([CB.T_hp_cs_ex,CB.T_hp_cs_su])-273.15,'x', color='darkslategrey',zorder=-5,linewidth=1.2)
    ax1.plot((np.asarray([CB.my_HE.s_he_4,CB.my_HE.s_he_1])-s_sat_he_min)/ds_sat_he,np.asarray([CB.T_he_cs_ex,CB.T_he_cs_su])-273.15,'--',color='darkslategrey',zorder=-5,linewidth=1.2, label='source/sink')
    ax1.plot((np.asarray([CB.my_HE.s_he_4,CB.my_HE.s_he_1])-s_sat_he_min)/ds_sat_he,np.asarray([CB.T_he_cs_ex,CB.T_he_cs_su])-273.15,'x', color='darkslategrey',zorder=-5,linewidth=1.2)
    
    # --- Axes ----------------------------------------------------------------
    T_min = min([  15,min(CB.my_HP.T_hp)-273.15,
                      min(CB.my_HE.T_he)-273.15])
    T_max = max([ 180,max(CB.my_HP.T_hp)-273.15,
                      max(CB.my_HE.T_he)-273.15])
    
    s_min = min([-0.2,(min(CB.my_HP.s_hp)-s_sat_hp_min)/ds_sat_hp,
                      (min(CB.my_HE.s_he)-s_sat_he_min)/ds_sat_he])
    s_max = max([ 1.2,(max(CB.my_HP.s_hp)-s_sat_hp_min)/ds_sat_hp,
                      (max(CB.my_HE.s_he)-s_sat_he_min)/ds_sat_he])
    
    ax1.plot(s_min-0.02,T_min-1,'o',color='w')
    ax1.plot(s_min-0.02,T_max-1,'o',color='w')
    ax1.plot(s_max+0.02,T_min+1,'o',color='w')
    ax1.plot(s_max+0.02,T_max+1,'o',color='w')
    
    # ax1.legend(loc='upper left', frameon=False)
    ax1.legend(loc=(0.05, 0.75), frameon=False)
    ax1.set_xlabel('$\mathrm{(s-s_{TES,min})/(s_{TES,max}-s_{TES,min})}$ [-]',fontsize=13, color='dimgray', loc='right')
    ax1.spines["bottom"].set_bounds(s_min,s_max)
    ax1.set_xticks(np.round([s_min,0.0,0.2,0.4,0.6,0.8,1.0,s_max],2))
    ax1.spines["bottom"].set_linewidth(1.5)
    ax1.spines["bottom"].set_color('dimgray')
    ax1.set_ylabel('T [°C]',fontsize=13, color='dimgray', rotation=0, loc='top')
    ax1.yaxis.set_label_coords(0.0,1.0)
    ax1.spines["left"].set_bounds(T_min,T_max)
    if trans_crit_hp and trans_crit_he:
        ax1.set_yticks([T_min,CB.my_HE.T_he_1x-273.15,CB.T_st_ht-273.15,CB.T_st_lt-273.15,CB.my_HP.T_hp_1x-273.15,T_max])
    elif trans_crit_hp:
        ax1.set_yticks([T_min,CB.my_HE.T_he_1x-273.15,CB.T_st_ht-273.15,CB.T_st_lt-273.15,CB.my_HP.T_hp_1x-273.15,CB.my_HE.T_he_2x-273.15,T_max])
        
    elif trans_crit_he:
        ax1.set_yticks([T_min,CB.my_HE.T_he_1x-273.15,CB.T_st_ht-273.15,CB.T_st_lt-273.15,CB.my_HP.T_hp_1x-273.15,CB.my_HP.T_hp_2x-273.15,T_max])
    else:    
        ax1.set_yticks([T_min,CB.my_HE.T_he_1x-273.15,CB.T_st_ht-273.15,CB.T_st_lt-273.15,CB.my_HP.T_hp_1x-273.15,CB.my_HE.T_he_2x-273.15,CB.my_HP.T_hp_2x-273.15,T_max])
    ax1.spines["left"].set_linewidth(1.5)
    ax1.spines["left"].set_color('dimgray')
    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)
    ax1.tick_params(axis='both', labelsize=12, colors='dimgray', width=1.5, length=8, direction='in')
    
    ax1.label_outer()
    fig.tight_layout()
    plt.show()
    
    out = fig
    return out
