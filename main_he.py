"""
- Date: March 2025
- Author: Antoine Laterre

@ Project: "CBSim - Carnot Battery Simulator" 

-> Main script for the simulation of heat engines

"""

#%% IMPORT PACKAGES

import CoolProp.CoolProp as CP
import time
import os
import sys
SIMULATOR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SIMULATOR,'src'))
import _module_heat_engine as HE

#%% TEST FUNCTION

def eval_HE(inputs,he):
    # 0 - Inputs --------------------------------------------------------------
    if he == 'SBORC':
        T_he_cs, dT_he_cs_gl, fluid_he_cs,\
        T_he_hs, dT_he_hs_gl, fluid_he_hs,\
        fluid_he,version,     post_proc, p_0, T_0                 = inputs
    if he == 'SRORC':
        T_he_cs, dT_he_cs_gl, fluid_he_cs,\
        T_he_hs, dT_he_hs_gl, fluid_he_hs,\
        fluid_he,version,     post_proc, p_0, T_0,        epsilon = inputs
    if he == 'TBORC':
        T_he_cs, dT_he_cs_gl, fluid_he_cs,\
        T_he_hs, dT_he_hs_gl, fluid_he_hs,\
        fluid_he,version,     post_proc, p_0, T_0, m_rat          = inputs
    if he == 'TRORC':
        T_he_cs, dT_he_cs_gl, fluid_he_cs,\
        T_he_hs, dT_he_hs_gl, fluid_he_hs,\
        fluid_he,version,     post_proc, p_0, T_0, m_rat, epsilon = inputs
    T_0 = T_0 + 273.15
    # 1.1 - Heat sink ---------------------------------------------------------
    p_he_cs_su      = 1.0e+5
    T_he_cs_su      = T_he_cs + 273.15
    p_he_cs_ex      = 1.0e+5
    T_he_cs_ex      = T_he_cs + 273.15 + dT_he_cs_gl
    i_he_cs_su      = CP.PropsSI('H','T',T_he_cs_su,
                                     'P',p_he_cs_su,fluid_he_cs)
    i_he_cs_ex      = CP.PropsSI('H','T',T_he_cs_ex,
                                     'P',p_he_cs_ex,fluid_he_cs)
    m_he_cs         = 0
    # 1.2 - Heat source ----------------------------------------------------
    p_he_hs_su      = 2.5e+5
    T_he_hs_su      = T_he_hs + 273.15
    p_he_hs_ex      = 2.5e+5
    T_he_hs_ex      = T_he_hs + 273.15 - dT_he_hs_gl
    i_he_hs_su      = CP.PropsSI('H','T',T_he_hs_su,
                                     'P',p_he_hs_su,fluid_he_hs)
    i_he_hs_ex      = CP.PropsSI('H','T',T_he_hs_ex,
                                     'P',p_he_hs_ex,fluid_he_hs)
    m_he_hs         = 0.3333333333333333
    # 1.3 - Heat engine--------------------------------------------------------
    P_he            = 5e+3
    # 1.4 - Heat engine parameters --------------------------------------------    
    eta_max_ex      = 0.75
    eta_pm          = 0.60
    dT_he_cd_pp     = 5.0
    dT_he_ev_pp     = 3.0
    dT_he_ev_sh     = 1.0
    dT_he_cd_sc     = 3.0
    dp_he_ev        = 0.0
    dp_he_cd        = 0.0
    dp_he_rg_vp     = 0.0
    dp_he_rg_lq     = 0.0
    wet_ex          = 0
    # 1.5 - Prepare the dictionnary -------------------------------------------
    my_HE = {'eta_max_ex':   eta_max_ex,
             'eta_pm':       eta_pm,
             'dT_he_cd_pp':  dT_he_cd_pp,
             'dT_he_ev_pp':  dT_he_ev_pp,
             'dT_he_ev_sh':  dT_he_ev_sh, 
             'dT_he_cd_sc':  dT_he_cd_sc,
             'dp_he_ev':     dp_he_ev,
             'dp_he_cd':     dp_he_cd,
             'version':      version,
             'mode':         'source',
             'p_ref':        p_he_cs_su,
             'T_ref':        T_he_cs_su,
             'p_0':          p_he_cs_su,
             'T_0':          T_he_cs_su,
             'fluid_he':     fluid_he,
             'wet_ex':       wet_ex,}
    if he == 'SRORC' or he == 'TRORC':
       my_HE['epsilon']  = epsilon 
       my_HE['dp_he_rg_lq'] = dp_he_rg_lq
       my_HE['dp_he_rg_vp'] = dp_he_rg_vp
    if he == 'TBORC' or he == 'TRORC':
       my_HE['m_rat']    = m_rat 
    # 1.6 - Prepare the options -----------------------------------------------
    print_flag, plot_flag = False, False
    debug                 = False
    exergy                = True
    options               = {'plot_flag':     plot_flag,
                             'print_flag':    print_flag,
                             'debug':         debug,
                             'exergy':        exergy,}    
    # 2.1 - Process the inputs ------------------------------------------------
    inputs      = p_he_cs_su, i_he_cs_su, p_he_cs_ex, i_he_cs_ex, m_he_cs, fluid_he_cs, \
                  p_he_hs_su, i_he_hs_su, p_he_hs_ex, i_he_hs_ex, m_he_hs, fluid_he_hs, \
                  P_he  
    parameters  = my_HE
    # 2.2 - Run the evaluation of the heat engine -----------------------------
    if he == 'SBORC':  my_HE = HE.SBORC(inputs, parameters, options)
    if he == 'SRORC':  my_HE = HE.SRORC(inputs, parameters, options)
    if he == 'TBORC':  my_HE = HE.TBORC(inputs, parameters, options)
    if he == 'TRORC':  my_HE = HE.TRORC(inputs, parameters, options)
    my_HE.evaluate()
    # 3.0 - Make the plots ----------------------------------------------------    
    if post_proc:
        print('--- Performance for '+he+' -----------------------------------')
        print('eta_cyclen   [%]: '+str(round(my_HE.eta_he_cyclen*100,3)))
        print('exp. ratio   [-]: '+str(round(my_HE.exp_ratio,        2)))
        print('p_min_he   [bar]: '+str(round(my_HE.p_he_1*1e-5,3)))
        print('p_max_he   [bar]: '+str(round(my_HE.p_he_2*1e-5,3)))
        if he == 'TBORC' or he == 'TRORC':
            print('m_rat [(kg/s)_hs/(kg/s)_cs]: '+str(round((my_HE.m_rat),6)))
        if version != 'thermodynamic_full':
            print('Q_in_he_en  [kW]: '+str(round(my_HE.m_he*my_HE.q_in_he_en *1e-3,2)))
            print('source    [kg/s]: '+str(round(my_HE.m_he_hs,                    2)))
            print('W_net_he_en [kW]: '+str(round(my_HE.m_he*my_HE.w_net_he_en*1e-3,2)))
            print('W_pmp_he_en [kW]: '+str(round(my_HE.m_he*my_HE.w_in_he_en *1e-3,2)))
            print('W_exp_he_en [kW]: '+str(round(my_HE.m_he*my_HE.w_out_he_en*1e-3,2)))
            print('Q_out_he_en [kW]: '+str(round(my_HE.m_he*my_HE.q_out_he_en*1e-3,2)))
            print('sink      [kg/s]: '+str(round(my_HE.m_he_cs,                    2)))
            
    from _module_plots import plot_Th_ORC
    plot_Th_ORC(my_HE)
    from _module_plots import plot_Ts_ORC
    plot_Ts_ORC(my_HE)

#%% SINGLE POINT ANALYSIS 

p_0         = 1e+5
T_0         = 15

T_he_cs     = 15
dT_he_cs_gl = 10
T_he_hs     = 120
dT_he_hs_gl = 60 

fluid_sbhe  = 'R1234ze(E)'
fluid_srhe  = 'R227ea'
fluid_tbhe  = 'R744'
fluid_trhe  = 'R744'

fluid_he_cs ='H2O'
fluid_he_hs ='H2O'

# version='operational_light'
version='thermodynamic_full'
post_proc=True

epsilon = 0.80
m_rat   = 0

# --- SBORC -------------------------------------------------------------------

inputs_sborc = T_he_cs, dT_he_cs_gl, fluid_he_cs,\
                T_he_hs, dT_he_hs_gl, fluid_he_hs,\
                fluid_sbhe,version,     post_proc, p_0, T_0

sta = time.time()
out = eval_HE(inputs_sborc,'SBORC')
end = time.time()
print(f'took {round(end-sta,3)} s')

# --- SRORC -------------------------------------------------------------------

inputs_srorc = T_he_cs, dT_he_cs_gl, fluid_he_cs,\
                T_he_hs, dT_he_hs_gl, fluid_he_hs,\
                fluid_srhe,version,     post_proc, p_0, T_0, epsilon

sta = time.time()
out = eval_HE(inputs_srorc,'SRORC')
end = time.time()
print(f'took {round(end-sta,3)} s')

# --- TBORC -------------------------------------------------------------------

inputs_tborc = T_he_cs, dT_he_cs_gl, fluid_he_cs,\
                T_he_hs, dT_he_hs_gl, fluid_he_hs,\
                fluid_tbhe,version,     post_proc, p_0, T_0, m_rat

sta = time.time()
out = eval_HE(inputs_tborc,'TBORC')
end = time.time()
print(f'took {round(end-sta,3)} s')

# --- TRORC -------------------------------------------------------------------

inputs_trorc = T_he_cs, dT_he_cs_gl, fluid_he_cs,\
                T_he_hs, dT_he_hs_gl, fluid_he_hs,\
                fluid_trhe,version,     post_proc, p_0, T_0, m_rat, epsilon

sta = time.time()
out = eval_HE(inputs_trorc,'TRORC')
end = time.time()
print(f'took {round(end-sta,3)} s')