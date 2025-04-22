"""
- Date: March 2025
- Author: Antoine Laterre

@ Project: "CBSim - Carnot Battery Simulator" 

-> Main script for the simulation of Carnot batteries

"""

#%% IMPORT PACKAGES

import CoolProp.CoolProp as CP
import time
import os
import sys
SIMULATOR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SIMULATOR,'src'))
import _module_carnot_battery as CB

#%% TEST THE CARNOT BATTERY

def eval_CB(T_hp_cs, dT_hp_cs_gl, T_he_cs, dT_he_cs_gl, T_st_ht, dT_st_sp, 
            dT_hp_cd_sc, dT_hp_ev_sh, dT_he_cd_sc, dT_he_ev_sh,
            fluid_hp, fluid_he, fluid_hp_cs, fluid_he_cs, fluid_st, 
            cb_type, version, post_proc):
    # 1.1 - Heat source -------------------------------------------------------
    p_hp_cs_su      = 1e+5
    T_hp_cs_su      = T_hp_cs + 273.15
    T_hp_cs_ex      = T_hp_cs_su - dT_hp_cs_gl
    m_hp_cs         = 1.0
    i_hp_cs_su      = CP.PropsSI('H','T',T_hp_cs_su,'P',p_hp_cs_su,fluid_hp_cs)
    i_hp_cs_ex      = CP.PropsSI('H','T',T_hp_cs_ex,'P',p_hp_cs_su,fluid_hp_cs)
    P_hp            = 1e+3
    # 1.2 - Heat sink ---------------------------------------------------------
    p_he_cs_su      = 1e+5 
    T_he_cs_su      = T_he_cs + 273.15
    T_he_cs_ex      = T_he_cs_su + dT_he_cs_gl
    m_he_cs         = 1.0
    i_he_cs_su      = CP.PropsSI('H','T',T_he_cs_su,'P',p_he_cs_su,fluid_he_cs)
    i_he_cs_ex      = CP.PropsSI('H','T',T_he_cs_ex,'P',p_he_cs_su,fluid_he_cs)
    P_he            = 1e+3
    # 1.3 - Carnot battery design ---------------------------------------------
    T_st_ht         = T_st_ht + 273.15    
    # 1.4 - Carnot battery parameters -----------------------------------------    
    eta_max_cp      = 0.75
    eta_max_ex      = 0.75
    eta_pm          = 0.60
    dT_hp_ev_pp     = 5.0
    dT_hp_cd_pp     = 3.0
    dT_he_cd_pp     = 5.0
    dT_he_ev_pp     = 3.0
    dp_hp_ev        = 0.00e+5
    dp_hp_cd        = 0.00e+5
    dp_he_ev        = 0.00e+5
    dp_he_cd        = 0.00e+5
    p_st_ht         = 2.5e+5
    p_st_lt         = 2.5e+5
    epsilon         = 0.80
    # 1.5 - Prepare the dictionnary -------------------------------------------
    my_CB_params = {
        'p_hp_cs_su':   p_hp_cs_su, 
        'i_hp_cs_su':   i_hp_cs_su,
        'i_hp_cs_ex':   i_hp_cs_ex,
        'p_he_cs_su':   p_he_cs_su, 
        'i_he_cs_su':   i_he_cs_su,
        'i_he_cs_ex':   i_he_cs_ex,
        'm_hp_cs':      m_hp_cs,
        'm_he_cs':      m_he_cs, 
        'p_st_ht':      p_st_ht, 
        'p_st_lt':      p_st_lt, 
        'T_st_ht':      T_st_ht, 
        'dT_st_sp':     dT_st_sp,
        'eta_max_cp':   eta_max_cp,
        'eta_max_ex':   eta_max_ex,
        'eta_pm':       eta_pm, 
        'dT_hp_ev_pp':  dT_hp_ev_pp,
        'dT_hp_cd_pp':  dT_hp_cd_pp,
        'dT_he_ev_pp':  dT_he_ev_pp,
        'dT_he_cd_pp':  dT_he_cd_pp,
        'dT_he_ev_sh':  dT_he_ev_sh,
        'dT_hp_ev_sh':  dT_hp_ev_sh, 
        'dT_he_cd_sc':  dT_he_cd_sc,
        'dT_hp_cd_sc':  dT_hp_cd_sc,
        'dp_hp_ev':     dp_hp_ev,
        'dp_hp_cd':     dp_hp_cd,
        'dp_hp_rg_lq':  (dp_hp_ev+dp_hp_cd)/2,
        'dp_hp_rg_vp':  (dp_hp_ev+dp_hp_cd)/2,
        'epsilon_hp':   epsilon,
        'dp_he_ev':     dp_he_ev,
        'dp_he_cd':     dp_he_cd,
        'dp_he_rg_lq':  (dp_he_ev+dp_he_cd)/2,
        'dp_he_rg_vp':  (dp_he_ev+dp_he_cd)/2,
        'epsilon_he':   epsilon,
        'm_hp_st_max':  0e+0,
        'm_he_st_max':  0e+0,
        'version':      version,
        'mode_hp':      True,
        'mode_he':      True,
        'mode':         'source',
        'p_ref':        p_he_cs_su,
        'T_ref':        T_he_cs_su,
        'p_0':          p_he_cs_su,
        'T_0':          T_he_cs_su,
        'fluid_hp':     fluid_hp,
        'fluid_he':     fluid_he,
        'fluid_st':     fluid_st,
        'wet_ex':       0,
        'm_rat_hp':     0,
        'm_rat_he':     0,
        }
    # 1.6 - Prepare the options -----------------------------------------------
    print_flag, plot_flag = False, False
    debug                 = False
    exergy                = True
    options               = {'plot_flag':     plot_flag,
                             'print_flag':    print_flag,
                             'debug':         debug,
                             'exergy':        exergy,} 
    # 2 - Evaluate the CB _____________________________________________________
    # 2.1 - Process the inputs ------------------------------------------------
    inputs =    p_hp_cs_su, i_hp_cs_su, p_hp_cs_su, i_hp_cs_ex, m_hp_cs, fluid_hp_cs, \
                p_he_cs_su, i_he_cs_su, p_he_cs_su, i_he_cs_ex, m_he_cs, fluid_he_cs, \
                P_hp, P_he
    parameters            = my_CB_params
    # 2.2 - Run the evaluation of the Carnot Battery --------------------------
    if cb_type == 0: my_CB = CB.SBVCHP_SBORC_STES2T(inputs, parameters, options)
    if cb_type == 1: my_CB = CB.SBVCHP_SRORC_STES2T(inputs, parameters, options)
    if cb_type == 2: my_CB = CB.SRVCHP_SBORC_STES2T(inputs, parameters, options)
    if cb_type == 3: my_CB = CB.SRVCHP_SRORC_STES2T(inputs, parameters, options)
    if cb_type == 4: my_CB = CB.SBVCHP_TBORC_STES2T(inputs, parameters, options)
    if cb_type == 5: my_CB = CB.TBVCHP_SBORC_STES2T(inputs, parameters, options)
    if cb_type == 6: my_CB = CB.SRVCHP_TBORC_STES2T(inputs, parameters, options)
    if cb_type == 7: my_CB = CB.TBVCHP_SRORC_STES2T(inputs, parameters, options)
    if cb_type == 8: my_CB = CB.TBVCHP_TBORC_STES2T(inputs, parameters, options)
    if cb_type == 9: my_CB = CB.SBVCHP_TRORC_STES2T(inputs, parameters, options)
    if cb_type == 10:my_CB = CB.TRVCHP_SBORC_STES2T(inputs, parameters, options)
    if cb_type == 11:my_CB = CB.SRVCHP_TRORC_STES2T(inputs, parameters, options)
    if cb_type == 12:my_CB = CB.TRVCHP_SRORC_STES2T(inputs, parameters, options)
    if cb_type == 13:my_CB = CB.TBVCHP_TRORC_STES2T(inputs, parameters, options)
    if cb_type == 14:my_CB = CB.TRVCHP_TBORC_STES2T(inputs, parameters, options)
    if cb_type == 15:my_CB = CB.TRVCHP_TRORC_STES2T(inputs, parameters, options)
    my_CB.evaluate()
    
    # 3 - Make the plots ______________________________________________________
    if post_proc and not my_CB.error:
        print('\n === HEAT PUMP === ')
        print('eta_cyclen [%]: '+str(round(my_CB.my_HP.eta_hp_cyclen*100,2)))
        print('eta_totex  [%]: '+str(round(my_CB.my_HP.eta_hp_totex *100,2)))
        print('comp.ratio [-]: '+str(round(my_CB.my_HP.comp_ratio,       2)))
        print('p_min_hp [bar]: '+str(round(my_CB.my_HP.p_hp_1*1e-5,3)))
        print('p_max_hp [bar]: '+str(round(my_CB.my_HP.p_hp_2*1e-5,3)))
        print(' === HEAT ENGINE === ')
        print('eta_cyclen [%]: '+str(round(my_CB.my_HE.eta_he_cyclen*100,2)))
        print('eta_totex  [%]: '+str(round(my_CB.my_HE.eta_he_totex *100,2)))
        print('exp. ratio [-]: '+str(round(my_CB.my_HE.exp_ratio,        2)))
        print('p_min_he [bar]: '+str(round(my_CB.my_HE.p_he_1*1e-5,3)))
        print('p_max_he [bar]: '+str(round(my_CB.my_HE.p_he_2*1e-5,3)))
        print(' === CARNOT BATTERY === ')
        print('eta_P2P [%]:  '+str(my_CB.eta_cb_elec*100))
        print('eta_II [%]:   '+str(my_CB.eta_cb_exer*100))
        print('[kWh_th/m³]:  '+str(my_CB.E_dens_th/3.6e+6))
        print('[kWh_el/m³]:  '+str(my_CB.E_dens_el/3.6e+6))
        print('psi_hp_carnot [%]: '+str(my_CB.my_HP.psi_hp_carnot*100))
        print('psi_he_carnot [%]: '+str(my_CB.my_HE.psi_he_carnot*100))
        
        if version != 'thermodynamic_full':
            print(' === HP OPERATIONS === ')
            print('Pel  [kW]:   '+str(my_CB.my_HP.W_in_hp_en*1e-3))
            print('Qin  [kW]:   '+str(my_CB.my_HP.Q_in_hp_en*1e-3))
            print('Qout [kW]:   '+str(my_CB.my_HP.Q_out_hp_en*1e-3))
            print('src  [kg/s]: '+str(my_CB.my_HP.m_hp_cs))
            print('wf   [kg/s]: '+str(my_CB.my_HP.m_hp))
            print('sink [kg/s]: '+str(my_CB.my_HP.m_hp_hs))
            print(' === HE OPERATIONS === ')
            print('Pel  [kW]:   '+str(my_CB.my_HE.W_net_he_en*1e-3))
            print('Qin  [kW]:   '+str(my_CB.my_HE.Q_in_he_en*1e-3))
            print('Qout [kW]:   '+str(my_CB.my_HE.Q_out_he_en*1e-3))
            print('src  [kg/s]: '+str(my_CB.my_HE.m_he_hs))
            print('wf   [kg/s]: '+str(my_CB.my_HE.m_he))
            print('sink [kg/s]: '+str(my_CB.my_HE.m_he_cs))
            
        from _module_plots import plot_Th_VCHP
        plot_Th_VCHP(my_CB.my_HP)
        from _module_plots import plot_Th_ORC
        plot_Th_ORC(my_CB.my_HE)
        from _module_plots import plot_Ts_VCHP_ORC_SHTES_norm
        plot_Ts_VCHP_ORC_SHTES_norm(my_CB)

# -----------------------------------------------------------------------------

T_st_ht     = 120
dT_st_sp    = 60

T_hp_cs     = 15
T_he_cs     = 15

dT_hp_cs_gl = 10
dT_he_cs_gl = 10

dT_hp_ev_sh = 5.0
dT_hp_cd_sc = 40.0
dT_he_ev_sh = 1.0
dT_he_cd_sc = 3.0

version='thermodynamic_full'
post_proc=True

fluid_hp_cs = 'H2O'
fluid_he_cs = 'H2O'
fluid_st    = 'H2O'

#%% SINGLE POINT ANALYSIS 

# --- SBVCHP + SBORC ----------------------------------------------------------
cb_type     = 0
fluid_hp    = 'R1234ze(Z)'
fluid_he    = 'R1234ze(E)'

sta = time.time()
out = eval_CB(T_hp_cs, 
              dT_hp_cs_gl, 
              T_he_cs, 
              dT_he_cs_gl, 
              T_st_ht, 
              dT_st_sp, 
              dT_hp_cd_sc, 
              dT_hp_ev_sh, 
              dT_he_cd_sc, 
              dT_he_ev_sh,
              fluid_hp, 
              fluid_he, 
              fluid_hp_cs, 
              fluid_he_cs, 
              fluid_st,
              cb_type,
              version, 
              post_proc)
end = time.time()
print(f'took {round(end-sta,3)} s')

# --- SRVCHP + SRORC ----------------------------------------------------------
cb_type     = 3
fluid_hp    = 'R1234ze(Z)'
fluid_he    = 'R227ea'

sta = time.time()
out = eval_CB(T_hp_cs, 
              dT_hp_cs_gl, 
              T_he_cs, 
              dT_he_cs_gl, 
              T_st_ht, 
              dT_st_sp, 
              dT_hp_cd_sc, 
              dT_hp_ev_sh, 
              dT_he_cd_sc, 
              dT_he_ev_sh,
              fluid_hp, 
              fluid_he, 
              fluid_hp_cs, 
              fluid_he_cs, 
              fluid_st,
              cb_type,
              version, 
              post_proc)
end = time.time()
print(f'took {round(end-sta,3)} s')

# --- TBVCHP + TBORC ----------------------------------------------------------
cb_type     = 8
fluid_hp    = 'R143a'
fluid_he    = 'R744'

sta = time.time()
out = eval_CB(T_hp_cs, 
              dT_hp_cs_gl, 
              T_he_cs, 
              dT_he_cs_gl, 
              T_st_ht, 
              dT_st_sp, 
              dT_hp_cd_sc, 
              dT_hp_ev_sh, 
              dT_he_cd_sc, 
              dT_he_ev_sh,
              fluid_hp, 
              fluid_he, 
              fluid_hp_cs, 
              fluid_he_cs, 
              fluid_st,
              cb_type,
              version, 
              post_proc)
end = time.time()
print(f'took {round(end-sta,3)} s')

# --- TRVCHP + TRORC ----------------------------------------------------------
cb_type     = 15
fluid_hp    = 'R143a'
fluid_he    = 'R744'

sta = time.time()
out = eval_CB(T_hp_cs, 
              dT_hp_cs_gl, 
              T_he_cs, 
              dT_he_cs_gl, 
              T_st_ht, 
              dT_st_sp, 
              dT_hp_cd_sc, 
              dT_hp_ev_sh, 
              dT_he_cd_sc, 
              dT_he_ev_sh,
              fluid_hp, 
              fluid_he, 
              fluid_hp_cs, 
              fluid_he_cs, 
              fluid_st,
              cb_type,
              version, 
              post_proc)
end = time.time()
print(f'took {round(end-sta,3)} s')
