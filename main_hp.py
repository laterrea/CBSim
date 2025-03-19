"""
- Date: March 2025
- Author: Antoine Laterre

@ Project: "CBSim - Carnot Battery Simulator" 

-> Main script for the simulation of vapour compression heat pumps

"""

#%% IMPORT PACKAGES

import CoolProp.CoolProp as CP
import time
import os
import sys
SIMULATOR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SIMULATOR,'src'))
import _module_heat_pump   as HP

#%% TEST FUNCTION

def eval_HP(inputs,hp):
    # 0 - Inputs --------------------------------------------------------------
    if hp == 'SBVCHP':
        T_hp_cs, dT_hp_cs_gl, fluid_hp_cs,\
        T_hp_hs, dT_hp_hs_gl, fluid_hp_hs,\
        fluid_hp,version,     post_proc, p_0, T_0 = inputs
    if hp == 'SRVCHP':
        T_hp_cs, dT_hp_cs_gl, fluid_hp_cs,\
        T_hp_hs, dT_hp_hs_gl, fluid_hp_hs,\
        fluid_hp,version,     post_proc, p_0, T_0, epsilon = inputs
    if hp == 'TBVCHP':
        T_hp_cs, dT_hp_cs_gl, fluid_hp_cs,\
        T_hp_hs, dT_hp_hs_gl, fluid_hp_hs,\
        fluid_hp,version,     post_proc, p_0, T_0, m_rat = inputs
    if hp == 'TRVCHP':
        T_hp_cs, dT_hp_cs_gl, fluid_hp_cs,\
        T_hp_hs, dT_hp_hs_gl, fluid_hp_hs,\
        fluid_hp,version,     post_proc, p_0, T_0, m_rat, epsilon = inputs
    T_0 = T_0 + 273.15
    # 1.1 - Heat source -------------------------------------------------------
    p_hp_cs_su      = 1.0e+5
    T_hp_cs_su      = T_hp_cs + 273.15
    p_hp_cs_ex      = 1.0e+5
    T_hp_cs_ex      = T_hp_cs + 273.15 - dT_hp_cs_gl
    i_hp_cs_su      = CP.PropsSI('H','T',T_hp_cs_su,
                                     'P',p_hp_cs_su,fluid_hp_cs)
    i_hp_cs_ex      = CP.PropsSI('H','T',T_hp_cs_ex,
                                     'P',p_hp_cs_ex,fluid_hp_cs)
    m_hp_cs         = 0
    # 1.2 - Heat sink ---------------------------------------------------------
    p_hp_hs_su      = 2.5e+5
    T_hp_hs_su      = T_hp_hs + 273.15
    p_hp_hs_ex      = 2.5e+5
    T_hp_hs_ex      = T_hp_hs_su + dT_hp_hs_gl
    i_hp_hs_su      = CP.PropsSI('H','T',T_hp_hs_su,
                                     'P',p_hp_hs_su,fluid_hp_hs)
    i_hp_hs_ex      = CP.PropsSI('H','T',T_hp_hs_ex,
                                     'P',p_hp_hs_ex,fluid_hp_hs)
    m_hp_hs         = 0.3333333333333333
    # 1.3 - Heat pump ---------------------------------------------------------
    P_hp            = 5e+3
    # 1.4 - Heat pump parameters ----------------------------------------------    
    eta_max_cp      = 0.75
    
    dT_hp_ev_pp     = 5.0
    dT_hp_cd_pp     = 3.0
    dT_hp_ev_sh     = 5.0
    dT_hp_cd_sc     = 40.0
    dp_hp_ev        = 0.0
    dp_hp_cd        = 0.0 
    dp_hp_rg_vp     = 0.0
    dp_hp_rg_lq     = 0.0
    
    # 1.5 - Prepare the dictionnary -------------------------------------------
    my_HP = {'eta_max_cp':   eta_max_cp,
             
             'dT_hp_ev_pp':  dT_hp_ev_pp,
             'dT_hp_cd_pp':  dT_hp_cd_pp,
             'dT_hp_ev_sh':  dT_hp_ev_sh, 
             'dT_hp_cd_sc':  dT_hp_cd_sc,
             'dp_hp_ev':     dp_hp_ev,
             'dp_hp_cd':     dp_hp_cd,
             'version':      version,
             'mode':         'sink',
             'p_ref':        p_hp_cs_ex,
             'T_ref':        T_hp_cs_ex,
             'p_0':          p_0,
             'T_0':          T_0,
             'fluid_hp':     fluid_hp,
             }  
    if hp == 'SRVCHP' or hp == 'TRVCHP':
        my_HP['epsilon']     = epsilon 
        my_HP['dp_hp_rg_vp'] = dp_hp_rg_vp
        my_HP['dp_hp_rg_lq'] = dp_hp_rg_lq
    if hp == 'TBVCHP' or hp == 'TRVCHP':
       my_HP['m_rat']        = m_rat 
    # 1.6 - Prepare the options -----------------------------------------------
    print_flag, plot_flag = False, False
    debug                 = False
    exergy                = True
    options               = {'plot_flag':     plot_flag,
                             'print_flag':    print_flag,
                             'debug':         debug,
                             'exergy':        exergy,}   
    # 2.1 - Process the inputs ------------------------------------------------
    inputs      = p_hp_cs_su, i_hp_cs_su, p_hp_cs_ex, i_hp_cs_ex, m_hp_cs, fluid_hp_cs, \
                  p_hp_hs_su, i_hp_hs_su, p_hp_hs_ex, i_hp_hs_ex, m_hp_hs, fluid_hp_hs, \
                  P_hp
    parameters  = my_HP
    # 2.2 - Run the evaluation of the heat pump -------------------------------
    if hp == 'SBVCHP':   my_HP = HP.SBVCHP (inputs, parameters, options)
    if hp == 'SRVCHP':   my_HP = HP.SRVCHP (inputs, parameters, options)
    if hp == 'TBVCHP':   my_HP = HP.TBVCHP (inputs, parameters, options)
    if hp == 'TRVCHP':   my_HP = HP.TRVCHP (inputs, parameters, options)
    my_HP.evaluate()
    # 3.0 - Make the plots ----------------------------------------------------
    if post_proc:
        print('--- Performance for '+hp+' -----------------------------------')
        my_HP.export_states()
        print('eta_cyclen    [%]: '+str(round(my_HP.eta_hp_cyclen*100,2)))
        print('comp. ratio   [-]: '+str(round(my_HP.comp_ratio,       2)))
        print('p_min       [bar]: '+str(round(my_HP.p_hp_1*1e-5,3)))
        print('p_max       [bar]: '+str(round(my_HP.p_hp_2*1e-5,3)))
        if hp == 'TBVCHP' or hp == 'TRVCHP':
            print('m_rat [(kg/s)_hs/(kg/s)_cs]: '+str(round((my_HP.m_rat),6)))
        if hp == 'SBVCHP' or hp == 'SRVCHP':
            print('subcooling    [K]: '+str(round(my_HP.T_hp_3x-my_HP.T_hp_3,2)))
        if version != 'thermodynamic_full':
            print('Q_in_hp_en   [kW]: '+str(round(my_HP.m_hp*my_HP.q_in_hp_en *1e-3,2)))
            print('source     [kg/s]: '+str(round(my_HP.m_hp_cs,                    2)))
            print('W_in_hp_en   [kW]: '+str(round(my_HP.W_in_hp_en *1e-3,2)))
            print('Q_out_hp_en  [kW]: '+str(round(my_HP.Q_out_hp_en*1e-3,2)))
            print('sink       [kg/s]: '+str(round(my_HP.m_hp_hs,         2)))
        from _module_plots import plot_Th_VCHP
        plot_Th_VCHP(my_HP)
        from _module_plots import plot_Ts_VCHP
        plot_Ts_VCHP(my_HP)

#%% SINGLE POINT ANALYSIS 

p_0             = 1e+5
T_0             = 15

T_hp_cs         = 15
dT_hp_cs_gl     = 10
T_hp_hs         = 60
dT_hp_hs_gl     = 60

fluid_hp_cs     = 'H2O'
fluid_hp_hs     = 'H2O'

version='thermodynamic_full'
# version= 'operational_light'
post_proc       = True

fluid_sbhp      = 'R1234ze(Z)'
fluid_srhp      = 'R1234ze(Z)'
fluid_tbhp      = 'R143a'
fluid_trhp      = 'R143a'

epsilon         = 0.80
m_rat           = 0.0

# --- SBVCHP ------------------------------------------------------------------

inputs_sbvchp = T_hp_cs,    dT_hp_cs_gl,    fluid_hp_cs,\
                T_hp_hs,    dT_hp_hs_gl,    fluid_hp_hs,\
                fluid_sbhp, version, post_proc, p_0, T_0

sta = time.time()
out = eval_HP(inputs_sbvchp,'SBVCHP')
end = time.time()
print(f'took {round(end-sta,3)} s')

# --- SRVCHP ------------------------------------------------------------------

inputs_srvchp = T_hp_cs,    dT_hp_cs_gl,    fluid_hp_cs,\
                T_hp_hs,    dT_hp_hs_gl,    fluid_hp_hs,\
                fluid_srhp, version, post_proc, p_0, T_0, epsilon

sta = time.time()
out = eval_HP(inputs_srvchp,'SRVCHP')
end = time.time()
print(f'took {round(end-sta,3)} s')

# --- TBVCHP ------------------------------------------------------------------

inputs_tbvchp = T_hp_cs,    dT_hp_cs_gl,    fluid_hp_cs,\
                T_hp_hs,    dT_hp_hs_gl,    fluid_hp_hs,\
                fluid_tbhp, version, post_proc, p_0, T_0, m_rat

sta = time.time()
out = eval_HP(inputs_tbvchp,'TBVCHP')
end = time.time()
print(f'took {round(end-sta,3)} s')

# --- TRVCHP ------------------------------------------------------------------

inputs_trvchp = T_hp_cs,    dT_hp_cs_gl,    fluid_hp_cs,\
                T_hp_hs,    dT_hp_hs_gl,    fluid_hp_hs,\
                fluid_trhp, version, post_proc, p_0, T_0, m_rat, epsilon

sta = time.time()
out = eval_HP(inputs_trvchp,'TRVCHP')
end = time.time()
print(f'took {round(end-sta,3)} s')