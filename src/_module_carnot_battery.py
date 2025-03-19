"""
- Date: March 2025
- Author: Antoine Laterre

@ Project: "CBSim - Carnot Battery Simulator" 

-> Module containing the Carnot batteries models

"""

#%% IMPORT PACKAGES 

import CoolProp.CoolProp as CP

import os
import sys
SIMULATOR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(SIMULATOR,'COMPONENTS'))
import _module_heat_pump   as HP
import _module_heat_engine as HE

#%% 0: SBVCHP_SBORC_STES2T

class SBVCHP_SBORC_STES2T:
    """
    Class for the simulation of Carnot Batteries based on:
        - Subcritical Basic Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Subcritical Basic Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    def __init__(self, inputs, parameters, options):
        """
        Create a SBVCHP_SBORC_SHS2T object.
        """
        
        p_hp_cs_su, i_hp_cs_su, p_hp_cs_ex, i_hp_cs_ex, m_hp_cs, fluid_hp_cs, \
        p_he_cs_su, i_he_cs_su, p_he_cs_ex, i_he_cs_ex, m_he_cs, fluid_he_cs, \
        P_hp, P_he = inputs        
        
        # --- Cold heat source ------------------------------------------------
        self.p_hp_cs_su  = p_hp_cs_su # [Pa]   cold heat source supply pressure
        self.i_hp_cs_su  = i_hp_cs_su # [J/kg] cold heat source supply enthalpy
        self.p_hp_cs_ex  = p_hp_cs_ex # [Pa]   cold heat source exit   pressure
        self.i_hp_cs_ex  = i_hp_cs_ex # [J/kg] cold heat source exit   enthalpy
        self.m_hp_cs     = m_hp_cs    # [kg/s] cold heat source mass flow rate
        self.fluid_hp_cs = fluid_hp_cs# [-]    cold heat source fluid
        self.T_hp_cs_su  = CP.PropsSI('T','P',self.p_hp_cs_su,
                                          'H',self.i_hp_cs_su,self.fluid_hp_cs)            
        self.T_hp_cs_ex  = CP.PropsSI('T','P',self.p_hp_cs_ex,
                                          'H',self.i_hp_cs_ex,self.fluid_hp_cs)
        
        # --- Cold heat sink --------------------------------------------------
        self.p_he_cs_su  = p_he_cs_su # [Pa]   cold heat sink   supply pressure 
        self.i_he_cs_su  = i_he_cs_su # [J/kg] cold heat sink   supply enthalpy
        self.p_he_cs_ex  = p_he_cs_ex # [Pa]   cold heat sink   exit   pressure 
        self.i_he_cs_ex  = i_he_cs_ex # [J/kg] cold heat sink   exit enthalpy
        self.m_he_cs     = m_he_cs    # [kg/s] cold heat sink   mass flow rate
        self.fluid_he_cs = fluid_he_cs# [-]    cold heat sink   fluid
        self.T_he_cs_su  = CP.PropsSI('T','P',self.p_he_cs_su,
                                          'H',self.i_he_cs_su,self.fluid_he_cs)
        self.T_he_cs_ex  = CP.PropsSI('T','P',self.p_he_cs_ex,
                                          'H',self.i_he_cs_ex,self.fluid_he_cs)
        
        # --- Carnot battery --------------------------------------------------
        self.P_hp        = P_hp
        self.P_he        = P_he
        
        # --- Parameters and options ------------------------------------------
        self.parameters  = parameters
        self.options     = options
    
    def evaluate(self):
        """
        Evaluate the SBVCHP_SBORC_SHS2T cycle for the given boundary 
        conditions.
        """
        
        self.error = True
        
        if self.parameters['version'] not in ['thermodynamic_full',
                                              'operational_light',
                                              'operational_full']:
            raise ValueError('An inconsistency was detected!\
                              In: '+os.path.join(os.path.abspath(__file__),
          'SBVCHP_SBORC_SHS2T','evaluate: parameters["version"] is not valid'))
        
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light'\
        or self.parameters['version'] == 'operational_full':
            if self.options['debug']:
                self.evaluate_cycle()
                self.error = False
            else:
                try:    
                    self.error = False
                    self.evaluate_cycle()
                    self.check_consistency()
                except: self.error = True

        if self.error == True: 
            raise ValueError('An inconsistency was detected in the cycle!\
                              In: '+os.path.join(os.path.abspath(__file__),
                              'SBVCHP_SBORC_SHS2T','check_consistency'))
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.SBVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.SBORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()
        
    # ======================================================================= #
    # ==== SUB-FUNCTIONS FOR CARNOT BATTERY MODEL =========================== #
    # ======================================================================= #
    
    def check_consistency(self):
        """
        Check that the results are consistent.
        """
        if self.my_HP.error:                self.error = True
        if self.my_HE.error:                self.error = True
        # if self.T_st_lt < self.T_hp_cs_ex:  self.error = True
        if self.T_st_ht < self.T_hp_cs_su:  self.error = True
    
    def retrieve_kpi(self):
        """
        Retrieve the energy KPI's of the Carnot battery.
        """
        
        # Replace with the retrieve KPI's from the subfunctions.
    
        self.eta_cb_elec    = self.my_HP.eta_hp_cyclen*self.my_HE.eta_he_cyclen
        # net HP exergy input: # W/kg/s|hp_sh
        pow_hp_supplex      = self.my_HP.pow_hp_supplex/self.my_HP.f_hp_cs
        # net HE exergy input: # W/kg/s|hp_sh
        w_net_he_en         = self.my_HE.w_net_he_en   /self.my_HE.f_he_hs\
                            * self.my_HP.f_hp_hs/self.my_HP.f_hp_cs
        self.eta_cb_exer    = w_net_he_en/pow_hp_supplex
        
        self.E_dens         = (self.i_st_ht-self.i_st_lt)             # J_th/kg
        self.E_dens_th      = (self.i_st_ht-self.i_st_lt)\
                            / (self.v_st_lt+self.v_st_ht)             # J_th/m³
        self.E_dens_el      = self.E_dens_th*self.my_HE.eta_he_cyclen # J_el/m³
        
        self.P_su_spe       = self.my_HP.w_in_hp_en/self.my_HP.f_hp_cs # W/kg/s
        
    def export_states(self):
        """
        Export the thermodynamic states of the Carnot battery necessary to draw
        the T-s and p-h diagrams.
        """
        
        self.my_HP.export_states()
        self.my_HE.export_states()
        
        self.p_cb = self.my_HP.p_hp + self.my_HE.p_he
        self.T_cb = self.my_HP.T_hp + self.my_HE.T_he
        self.i_cb = self.my_HP.i_hp + self.my_HE.i_he
        self.s_cb = self.my_HP.s_hp + self.my_HE.s_he
        self.x_cb = self.my_HP.x_hp + self.my_HE.x_he
        self.e_cb = self.my_HP.e_hp + self.my_HE.e_he
        
        self.T_st = self.T_st_lt,self.T_st_ht,self.T_hp_cs_su,self.T_hp_cs_ex,\
                    self.T_he_cs_su,self.T_he_cs_ex
        self.i_st = self.i_st_lt,self.i_st_ht,self.i_hp_cs_su,self.i_hp_cs_ex,\
                    self.i_he_cs_su,self.i_he_cs_ex
        
        out = (self.p_cb,self.T_cb,self.i_cb,self.s_cb,self.x_cb,self.e_cb), \
              (self.T_st,self.i_st)
        return out
    
    def get_operational_params(self):
        """
        Evaluates the operational parameters of the Carnot battery.
        """
        if self.parameters['mode_hp']:     
            self.V_hp_st = self.my_HP.m_hp_hs*self.v_st_lt # [m³/s]
            self.m_hp_cs = self.my_HP.m_hp_cs # [kg/s]
            self.m_hp_hs = self.my_HP.m_hp_hs # [kg/s]
        if not self.parameters['mode_hp']: 
            self.V_hp_st = 0
            self.m_hp_cs = 0
            self.m_hp_hs = 0
        if self.parameters['mode_he']:     
            self.V_he_st = self.my_HE.m_he_hs*self.v_st_lt # [m³/s]
            self.m_he_cs = self.my_HE.m_he_cs # [kg/s]
            self.m_he_hs = self.my_HE.m_he_hs # [kg/s]
        if not self.parameters['mode_he']: 
            self.V_he_st = 0
            self.m_he_cs = 0
            self.m_he_hs = 0

#%% 1: SBVCHP_SRORC_STES2T

class SBVCHP_SRORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Subcritical Basic Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Subcritical Recuperated Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'], 
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.SBVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'dp_he_rg_lq':  self.parameters['dp_he_rg_lq'],
                              'dp_he_rg_vp':  self.parameters['dp_he_rg_vp'],
                              'epsilon':      self.parameters['epsilon_he'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.SRORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()

#%% 2: SRVCHP_SBORC_STES2T

class SRVCHP_SBORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Subcritical Recuperated Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Subcritical Basic Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'dp_hp_rg_lq':  self.parameters['dp_hp_rg_lq'],
                              'dp_hp_rg_vp':  self.parameters['dp_hp_rg_vp'],
                              'epsilon':      self.parameters['epsilon_hp'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.SRVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.SBORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()
    
#%% 3: SRVCHP_SRORC_STES2T

class SRVCHP_SRORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Subcritical Recuperated Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Subcritical Recuperated Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'], 
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'dp_hp_rg_lq':  self.parameters['dp_hp_rg_lq'],
                              'dp_hp_rg_vp':  self.parameters['dp_hp_rg_vp'],
                              'epsilon':      self.parameters['epsilon_hp'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.SRVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'dp_he_rg_lq':  self.parameters['dp_he_rg_lq'],
                              'dp_he_rg_vp':  self.parameters['dp_he_rg_vp'],
                              'epsilon':      self.parameters['epsilon_he'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.SRORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()
            
#%% 4: SBVCHP_TBORC_STES2T

class SBVCHP_TBORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Subcritical Basic Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Transcritical Basic Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------

        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.SBVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()

        # --- HEAT ENGINE CYCLE -----------------------------------------------

        if self.parameters['mode_he']:

            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],
                              'm_rat':        self.parameters['m_rat_he'],}

            options_he     = self.options
            self.my_HE     = HE.TBORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()

        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()

#%% 5: TBVCHP_SBORC_STES2T

class TBVCHP_SBORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Transcritical Basic Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Subcritical Basic Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],
                              'm_rat':        self.parameters['m_rat_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.TBVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'dp_he_rg_lq':  self.parameters['dp_he_rg_lq'],
                              'dp_he_rg_vp':  self.parameters['dp_he_rg_vp'],
                              'epsilon':      self.parameters['epsilon_he'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.SBORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()

#%% 6: SRVCHP_TBORC_STES2T

class SRVCHP_TBORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Subcritical Recuperated Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Transcritical Basic Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'dp_hp_rg_lq':  self.parameters['dp_hp_rg_lq'],
                              'dp_hp_rg_vp':  self.parameters['dp_hp_rg_vp'],
                              'epsilon':      self.parameters['epsilon_hp'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.SRVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],
                              'm_rat':        self.parameters['m_rat_he'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.TBORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()

#%% 7: TBVCHP_SRORC_STES2T

class TBVCHP_SRORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Transcritical Basic Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Transcritical Basic Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],
                              'm_rat':        self.parameters['m_rat_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.TBVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'], 
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'dp_he_rg_lq':  self.parameters['dp_he_rg_lq'],
                              'dp_he_rg_vp':  self.parameters['dp_he_rg_vp'],
                              'epsilon':      self.parameters['epsilon_he'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.SRORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()

#%% 8: TBVCHP_TBORC_STES2T

class TBVCHP_TBORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Transcritical Basic Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Transcritical Basic Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],
                              'm_rat':        self.parameters['m_rat_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.TBVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],
                              'm_rat':        self.parameters['m_rat_he'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.TBORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()

#%% 9: SBVCHP_TRORC_STES2T

class SBVCHP_TRORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Subcritical Basic Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Transcritical Recuperated Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']

        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.SBVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'dp_he_rg_lq':  self.parameters['dp_he_rg_lq'],
                              'dp_he_rg_vp':  self.parameters['dp_he_rg_vp'],
                              'epsilon':      self.parameters['epsilon_he'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],
                              'm_rat':        self.parameters['m_rat_he'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.TRORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()

#%% 10: TRVCHP_SBORC_STES2T

class TRVCHP_SBORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Transcritical Recuperated Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Subcritical Basic Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'dp_hp_rg_lq':  self.parameters['dp_hp_rg_lq'],
                              'dp_hp_rg_vp':  self.parameters['dp_hp_rg_vp'],
                              'epsilon':      self.parameters['epsilon_hp'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],
                              'm_rat':        self.parameters['m_rat_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.TRVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'dp_he_rg_lq':  self.parameters['dp_he_rg_lq'],
                              'dp_he_rg_vp':  self.parameters['dp_he_rg_vp'],
                              'epsilon':      self.parameters['epsilon_he'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.SBORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()

#%% 11: SRVCHP_TRORC_STES2T

class SRVCHP_TRORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Subcritical Recuperated Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Transcritical Recuperated Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'dp_hp_rg_lq':  self.parameters['dp_hp_rg_lq'],
                              'dp_hp_rg_vp':  self.parameters['dp_hp_rg_vp'],
                              'epsilon':      self.parameters['epsilon_hp'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.SRVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'], 
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'dp_he_rg_lq':  self.parameters['dp_he_rg_lq'],
                              'dp_he_rg_vp':  self.parameters['dp_he_rg_vp'],
                              'epsilon':      self.parameters['epsilon_he'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],
                              'm_rat':        self.parameters['m_rat_he'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.TRORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()

#%% 12: TRVCHP_SRORC_STES2T

class TRVCHP_SRORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Transcritical Recuperated Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Subcritical Recuperated Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'dp_hp_rg_lq':  self.parameters['dp_hp_rg_lq'],
                              'dp_hp_rg_vp':  self.parameters['dp_hp_rg_vp'],
                              'epsilon':      self.parameters['epsilon_hp'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],
                              'm_rat':        self.parameters['m_rat_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.TRVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'dp_he_rg_lq':  self.parameters['dp_he_rg_lq'],
                              'dp_he_rg_vp':  self.parameters['dp_he_rg_vp'],
                              'epsilon':      self.parameters['epsilon_he'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.SRORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()


#%% 13: TBVCHP_TRORC_STES2T

class TBVCHP_TRORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Transcritical Basic Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Transcritical Recuperated Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],
                              'm_rat':        self.parameters['m_rat_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.TBVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'dp_he_rg_lq':  self.parameters['dp_he_rg_lq'],
                              'dp_he_rg_vp':  self.parameters['dp_he_rg_vp'],
                              'epsilon':      self.parameters['epsilon_he'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],
                              'm_rat':        self.parameters['m_rat_he'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.TRORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()

#%% 14: TRVCHP_TBORC_STES2T

class TRVCHP_TBORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Transcritical Recuperated Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Transcritical Basic Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'dp_hp_rg_lq':  self.parameters['dp_hp_rg_lq'],
                              'dp_hp_rg_vp':  self.parameters['dp_hp_rg_vp'],
                              'epsilon':      self.parameters['epsilon_hp'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],
                              'm_rat':        self.parameters['m_rat_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.TRVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'dp_he_rg_lq':  self.parameters['dp_he_rg_lq'],
                              'dp_he_rg_vp':  self.parameters['dp_he_rg_vp'],
                              'epsilon':      self.parameters['epsilon_he'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],
                              'm_rat':        self.parameters['m_rat_he'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.TBORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()

#%% 15: TRVCHP_TRORC_STES2T

class TRVCHP_TRORC_STES2T(SBVCHP_SBORC_STES2T):
    """
    Class for the simulation of Carnot Batteries based on:
        - Transcritical Recuperated Vapor Compression Heat Pumps
        - Sensible Thermal Energy Storage in 2 Tanks
        - Transcritical Recuperated Organic Rankine Cycles

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (Carnot battery boundaries and 
                                        operational parameters).
    'parameters':
        -> Dictionary with the model parameters.
    'options':
        -> Dictionary with the simulation options.
        
    Methods
    -------
    'evaluate':
        -> Compute and returns the model outputs.
    """
    
    # ======================================================================= #
    # ==== CARNOT BATTERY MODEL ============================================= #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the source
            * 'sink':   available heat flow rate at the sink
            * 'power':  electrical power input/output
            * 'speed':  compressor/expander speed (only for operational_full)
        
        Evaluation time
        ---------------
            ~ x.x to x.x ms per run (thermodynamic)
            ~ x.x to x.x ms per run (operational_light)
            ~ x.x to x.x ms per run (operational_full)
        """        
        
        # --- STORAGE SYSTEM --------------------------------------------------
        
        self.T_st_ht = self.parameters['T_st_ht']
        self.dT_st_sp= self.parameters['dT_st_sp']
        self.T_st_lt = self.T_st_ht - self.dT_st_sp
        self.p_st_ht = self.parameters['p_st_ht']
        self.p_st_lt = self.parameters['p_st_lt']
        self.i_st_lt = CP.PropsSI('H','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.v_st_lt = CP.PropsSI('D','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_lt = CP.PropsSI('S','T',self.T_st_lt,'P',self.p_st_lt,
                                          self.parameters['fluid_st'])
        self.i_st_ht = CP.PropsSI('H','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        self.v_st_ht = CP.PropsSI('D','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])**(-1)
        self.s_st_ht = CP.PropsSI('S','T',self.T_st_ht,'P',self.p_st_ht,
                                          self.parameters['fluid_st'])
        
        self.m_hp_st = self.parameters['m_hp_st_max']
        self.m_he_st = self.parameters['m_he_st_max']
        
        # --- HEAT PUMP CYCLE -------------------------------------------------
        
        if self.parameters['mode_hp']:
        
            inputs_hp      = self.p_hp_cs_su,self.i_hp_cs_su,self.p_hp_cs_ex, \
                             self.i_hp_cs_ex,self.m_hp_cs,   self.fluid_hp_cs,\
                             self.p_st_lt,   self.i_st_lt,   self.p_st_ht,    \
                             self.i_st_ht,   self.m_hp_st,                    \
                             self.parameters['fluid_st'], self.P_hp
                          
            parameters_hp  = {'eta_max_cp':   self.parameters['eta_max_cp'],
                              'dT_hp_cd_pp':  self.parameters['dT_hp_cd_pp'],
                              'dT_hp_ev_pp':  self.parameters['dT_hp_ev_pp'],
                              'dT_hp_ev_sh':  self.parameters['dT_hp_ev_sh'], 
                              'dT_hp_cd_sc':  self.parameters['dT_hp_cd_sc'],
                              'dp_hp_ev':     self.parameters['dp_hp_ev'],
                              'dp_hp_cd':     self.parameters['dp_hp_cd'],
                              'dp_hp_rg_lq':  self.parameters['dp_hp_rg_lq'],
                              'dp_hp_rg_vp':  self.parameters['dp_hp_rg_vp'],
                              'epsilon':      self.parameters['epsilon_hp'],
                              'p_ref':        self.parameters['p_ref'],
                              'T_ref':        self.parameters['T_ref'],
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_hp':     self.parameters['fluid_hp'],
                              'm_rat':        self.parameters['m_rat_hp'],}
            
            options_hp     = self.options
            
            self.my_HP     = HP.TRVCHP(inputs_hp, parameters_hp, options_hp)
            self.my_HP.evaluate()
        
        # --- HEAT ENGINE CYCLE -----------------------------------------------
        
        if self.parameters['mode_he']:
        
            inputs_he      = self.p_he_cs_su,self.i_he_cs_su,self.p_he_cs_ex, \
                             self.i_he_cs_ex,self.m_he_cs,   self.fluid_he_cs,\
                             self.p_st_ht,   self.i_st_ht,   self.p_st_lt,    \
                             self.i_st_lt,    self.m_he_st,                   \
                             self.parameters['fluid_st'], self.P_he 
            
            parameters_he  = {'eta_max_ex':   self.parameters['eta_max_ex'],
                              'eta_pm':       self.parameters['eta_pm'],
                              'dT_he_cd_pp':  self.parameters['dT_he_cd_pp'],
                              'dT_he_ev_pp':  self.parameters['dT_he_ev_pp'],
                              'dT_he_ev_sh':  self.parameters['dT_he_ev_sh'], 
                              'dT_he_cd_sc':  self.parameters['dT_he_cd_sc'],
                              'dp_he_ev':     self.parameters['dp_he_ev'],
                              'dp_he_cd':     self.parameters['dp_he_cd'],
                              'dp_he_rg_lq':  self.parameters['dp_he_rg_lq'],
                              'dp_he_rg_vp':  self.parameters['dp_he_rg_vp'],
                              'epsilon':      self.parameters['epsilon_he'],
                              'p_ref':        self.parameters['p_st_lt'],
                              'T_ref':        self.T_st_lt,
                              'p_0':          self.parameters['p_0'],
                              'T_0':          self.parameters['T_0'],
                              'version':      self.parameters['version'],
                              'mode':         self.parameters['mode'],
                              'fluid_he':     self.parameters['fluid_he'],
                              'wet_ex':       self.parameters['wet_ex'],
                              'm_rat':        self.parameters['m_rat_he'],}
            
            options_he     = self.options
            
            self.my_HE     = HE.TRORC(inputs_he, parameters_he, options_he)
            self.my_HE.evaluate()
        
        # --- CARNOT BATTERY ANALYSIS -----------------------------------------
        
        self.retrieve_kpi()
        
        if self.parameters['version'] == 'operational_full'\
        or self.parameters['version'] == 'operational_light':
            self.get_operational_params()