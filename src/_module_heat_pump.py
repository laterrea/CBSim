"""
- Date: March 2025
- Author: Antoine Laterre

@ Project: "CBSim - Carnot Battery Simulator" 

-> Module containing the heat pumps models

"""

#%% IMPORT PACKAGES

import scipy.optimize
import numpy as np

import CoolProp
import CoolProp.CoolProp as CP
from   CoolProp import AbstractState

import os
import sys
SIMULATOR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(SIMULATOR,'COMPONENTS'))

#%% MODEL: SBVCHP

class SBVCHP:
    """
    Class for the simulation of Subcritical Basic Vapor Compression Heat Pumps.

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (heat pump boundaries and operational 
                                        parameters).
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
        Create a SBVCHP object.
        """
        
        p_hp_cs_su, i_hp_cs_su, p_hp_cs_ex, i_hp_cs_ex, m_hp_cs, fluid_hp_cs, \
        p_hp_hs_su, i_hp_hs_su, p_hp_hs_ex, i_hp_hs_ex, m_hp_hs, fluid_hp_hs, \
        P_hp = inputs  
        
        # --- Heat pump -------------------------------------------------------
        self.P_hp        = P_hp       # [W]   heat pump power input
        backend          = 'TTSE'     # HEOS for Helm. Eq. Of State (slower)
        self.state_hp    = AbstractState(backend,parameters['fluid_hp'])
        self.state_hp.unspecify_phase()
        
        # --- Heat source -----------------------------------------------------
        self.p_hp_cs_su  = p_hp_cs_su # [Pa]   cold heat source supply pressure
        self.i_hp_cs_su  = i_hp_cs_su # [J/kg] cold heat source supply enthalpy
        self.p_hp_cs_ex  = p_hp_cs_ex # [Pa]   cold heat source exit   pressure
        self.i_hp_cs_ex  = i_hp_cs_ex # [J/kg] cold heat source exit   enthalpy
        self.m_hp_cs     = m_hp_cs    # [kg/s] cold heat source mass flow rate
        self.fluid_hp_cs = fluid_hp_cs# [-]    cold heat source fluid
        self.state_cs    = AbstractState(backend,self.fluid_hp_cs)
        self.state_cs.unspecify_phase()
        
        self.state_cs.update(CoolProp.HmassP_INPUTS,self.i_hp_cs_su,
                                                    self.p_hp_cs_su)
        self.T_hp_cs_su  = self.state_cs.T()
        self.s_hp_cs_su  = self.state_cs.smass()
        self.state_cs.update(CoolProp.HmassP_INPUTS,self.i_hp_cs_ex,
                                                    self.p_hp_cs_ex)
        self.T_hp_cs_ex  = self.state_cs.T()
        self.s_hp_cs_ex  = self.state_cs.smass()
        
        if self.T_hp_cs_su < self.T_hp_cs_ex:
            raise ValueError('The glide on the evaporator must be positive!\
                             In: '+os.path.join(os.path.abspath(__file__),
                             'SBVCHP','__init__'))
        
        # --- Heat sink -------------------------------------------------------
        self.p_hp_hs_su  = p_hp_hs_su # [Pa]   hot  heat sink   supply pressure
        self.i_hp_hs_su  = i_hp_hs_su # [J/kg] hot  heat sink   supply enthalpy
        self.p_hp_hs_ex  = p_hp_hs_ex # [Pa]   hot  heat sink   exit pressure
        self.i_hp_hs_ex  = i_hp_hs_ex # [J/kg] hot  heat sink   exit   enthalpy
        self.m_hp_hs     = m_hp_hs    # [kg/s] hot  heat sink   mass flow rate
        self.fluid_hp_hs = fluid_hp_hs# [-]    hot  heat sink   fluid
        self.state_hs    = AbstractState(backend,self.fluid_hp_hs)
        self.state_hs.unspecify_phase()
        
        self.state_hs.update(CoolProp.HmassP_INPUTS,self.i_hp_hs_su,
                                                    self.p_hp_hs_su)
        self.T_hp_hs_su  = self.state_hs.T()
        self.s_hp_hs_su  = self.state_hs.smass()
        self.state_hs.update(CoolProp.HmassP_INPUTS,self.i_hp_hs_ex,
                                                    self.p_hp_hs_ex)
        self.T_hp_hs_ex  = self.state_hs.T()
        self.s_hp_hs_ex  = self.state_hs.smass()

        if self.T_hp_hs_ex < self.T_hp_hs_su:
            raise ValueError('The glide on the condenser must be positive!\
                             In: '+os.path.join(os.path.abspath(__file__),
                             'SBVCHP','__init__'))
        
        # --- Parameters and options ------------------------------------------
        self.parameters  = parameters
        self.options     = options
    
    def evaluate(self):
        """
        Evaluate the SBVCHP cycle for the given boundary conditions.
        """
        
        self.error = True
        
        if self.parameters['version'] not in ['thermodynamic_full',
                                              'operational_light']:
            raise ValueError('An inconsistency was detected!\
                              In: '+os.path.join(os.path.abspath(__file__),
                'SBVCHP','evaluate: parameters["version"] is not valid'))
        
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
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
                              'SBVCHP','check_consistency'))
    
    # ======================================================================= #
    # ==== HEAT PUMP MODEL ================================================== #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the cold source
            * 'sink':   available heat flow rate at the hot sink
            * 'power':  electrical power supplied to the compressor
        
        """        
        
        # --- Cycle computation -----------------------------------------------
        
        # ..... two-phase .....................................................
        self.T_hp_4x = min(  self.T_hp_cs_su
                           - self.parameters['dT_hp_ev_pp']
                           - self.parameters['dT_hp_ev_sh'],
                             self.T_hp_cs_ex
                           - self.parameters['dT_hp_ev_pp'])
        if self.T_hp_4x > self.state_hp.T_critical():
            raise ValueError('An inconsistency was detected in the cycle!\
                              In: '+os.path.join(os.path.abspath(__file__),
                              'SBVCHP','evaluate_cycle:\
                    the evaporation temperature is above the critical point!'))
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.QT_INPUTS,0.5,self.T_hp_4x)
        self.p_hp_4x = self.state_hp.p()
        self.state_hp.unspecify_phase()
        # ..... saturated (gas) ...............................................
        self.p_hp_1x = self.p_hp_4x*(1-self.parameters['dp_hp_ev']) 
        self.x_hp_1x = 1.0
        self.state_hp.update(CoolProp.PQ_INPUTS,self.p_hp_1x,self.x_hp_1x)
        self.T_hp_1x = self.state_hp.T()
        self.i_hp_1x = self.state_hp.hmass()
        self.s_hp_1x = self.state_hp.smass()
        # ..... gas ...........................................................
        self.T_hp_1  = self.T_hp_1x+self.parameters['dT_hp_ev_sh']
        self.p_hp_1  = self.p_hp_1x
        if  self.p_hp_1 < self.state_hp.p_critical()\
        and self.T_hp_1 < self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_gas)
        if  self.p_hp_1 < self.state_hp.p_critical()\
        and self.T_hp_1 > self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_hp.update(CoolProp.PT_INPUTS,self.p_hp_1,self.T_hp_1)
        self.i_hp_1  = self.state_hp.hmass()
        self.x_hp_1  = self.state_hp.Q()
        self.s_hp_1  = self.state_hp.smass()
        self.v_hp_1  = self.state_hp.rhomass()**(-1)
        self.state_hp.unspecify_phase()
        # ..... saturated (gas) ...............................................
        self.p_hp_2x = self.find_p()
        self.x_hp_2x = 1.0
        self.state_hp.update(CoolProp.PQ_INPUTS,self.p_hp_2x,self.x_hp_2x)
        self.T_hp_2x = self.state_hp.T()
        self.i_hp_2x = self.state_hp.hmass()       
        self.s_hp_2x = self.state_hp.smass()
        # ..... gas ...........................................................
        self.p_hp_2  = self.p_hp_2x
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_hp.update(CoolProp.PSmass_INPUTS,
                                    self.p_hp_2,self.s_hp_1)
            self.i_hp_2s = self.state_hp.hmass()
            self.i_hp_2  = self.i_hp_1 \
                         +(self.i_hp_2s-self.i_hp_1)\
                         / self.parameters['eta_max_cp']
            self.m_hp    = 0
            self.P_cp    = 0
            self.eta_is_cp = self.parameters['eta_max_cp']
            self.cp_flag = 'na'
            self.state_hp.update(CoolProp.HmassP_INPUTS,
                                    self.i_hp_2,self.p_hp_2)
            self.T_hp_2  = self.state_hp.T()
        self.s_hp_2  = self.state_hp.smass()
        self.x_hp_2  = self.state_hp.Q()
        self.v_hp_2  = self.state_hp.rhomass()**(-1)
        # ..... saturated (liquid) ............................................
        self.p_hp_3x = self.p_hp_2x*(1-self.parameters['dp_hp_cd'])
        self.x_hp_3x = 0.0
        self.state_hp.update(CoolProp.PQ_INPUTS,self.p_hp_3x,self.x_hp_3x)
        self.T_hp_3x = self.state_hp.T()
        self.i_hp_3x = self.state_hp.hmass()
        self.s_hp_3x = self.state_hp.smass()
        # ..... liquid ........................................................
        self.T_hp_3  = self.T_hp_3x-self.parameters['dT_hp_cd_sc']
        self.p_hp_3  = self.p_hp_3x
        self.state_hp.specify_phase(CoolProp.iphase_liquid)
        self.state_hp.update(CoolProp.PT_INPUTS,self.p_hp_3,self.T_hp_3)
        self.i_hp_3  = self.state_hp.hmass()
        self.s_hp_3  = self.state_hp.smass()
        self.x_hp_3  = self.state_hp.Q()
        self.state_hp.unspecify_phase()
        # ..... two-phase .....................................................
        self.p_hp_4  = self.p_hp_4x 
        self.i_hp_4  = self.i_hp_3
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.HmassP_INPUTS,self.i_hp_4,self.p_hp_4)
        self.T_hp_4  = self.state_hp.T()
        self.s_hp_4  = self.state_hp.smass()
        self.x_hp_4  = self.state_hp.Q()
        self.state_hp.unspecify_phase()
        self.i_hp_4x = self.i_hp_4
        self.T_hp_4x = self.T_hp_4
        self.s_hp_4x = self.s_hp_4
        self.x_hp_4x = self.x_hp_4
        
        # --- Evaluate the operational conditions -----------------------------
        
        if self.parameters['version'] == 'operational_light':
            self.m_hp    = self.find_mass_flow()
            self.P_cp    = self.m_hp    * (self.i_hp_2 - self.i_hp_1)
            self.m_hp_cs = self.find_flow_secondary_cs()
            self.m_hp_hs = self.find_flow_secondary_hs()
        
        # --- Evaluate the Key Performance Indicators -------------------------
        
        self.retrieve_kpi()
    
    # ======================================================================= #
    # ==== SUB-FUNCTIONS FOR HEAT PUMP MODEL ================================ #
    # ======================================================================= #
    
    def find_p(self):
        """
        Find the condensing pressure in the SBVCHP based on the pinch-point and 
        for the selected subcooling.
        Note that phase change is not permitted for the secondary fluid.
        """      
        # --- Minimum condenser pressure (firm) ------------------------------- 
        T_hp_3x_min = self.T_hp_hs_su \
                    + self.parameters['dT_hp_cd_pp'] \
                    + self.parameters['dT_hp_cd_sc']
        if T_hp_3x_min > self.state_hp.T_critical():
            raise ValueError('An inconsistency was detected in the cycle!\
                              In: '+os.path.join(os.path.abspath(__file__),
                              'SBVCHP','find_p:\
            the min. condensation temperature is above the critical point!'))
        self.state_hp.update(CoolProp.QT_INPUTS,0.0,T_hp_3x_min)  
        self.p_hp_3x_min = self.state_hp.p()
        self.p_hp_2x_min = self.p_hp_3x_min/(1-self.parameters['dp_hp_cd'])
        # --- Maximum condenser pressure (indicative) -------------------------
        T_hp_2x_max = self.T_hp_hs_ex \
                    + self.parameters['dT_hp_cd_pp']
        if T_hp_2x_max >= self.state_hp.T_critical() - 1:
            T_hp_2x_max = self.state_hp.T_critical() - 1
        self.state_hp.update(CoolProp.QT_INPUTS,1.0,T_hp_2x_max)
        p_hp_2x_max = self.state_hp.p()
        # --- Find the pressure level -----------------------------------------
        guess       = p_hp_2x_max
        try:
            p_hp_2x =    scipy.optimize.fsolve(self.resi_p, guess)[0]
        except:
            self.solutions_p = []
            self.residuals_p = []
            res = scipy.optimize.least_squares(self.resi_p, guess, 
                                        bounds = (self.p_hp_2x_min,
                                    max(min(self.p_hp_2x_min*1.5,
                                            self.state_hp.p_critical()*0.99),
                                            self.state_hp.p_critical()*0.99)))
            p_hp_2x = res.x
            imin = np.argmin(self.residuals_p)
            if p_hp_2x != self.solutions_p[imin]:
                p_hp_2x = self.solutions_p[imin]
        
        if p_hp_2x < self.p_hp_2x_min\
        or p_hp_2x > self.state_hp.p_critical(): 
            self.solutions_p = []
            self.residuals_p = []
            res = scipy.optimize.least_squares(self.resi_p, guess, 
                                        bounds = (self.p_hp_2x_min,
                                    max(min(self.p_hp_2x_min*1.5,
                                            self.state_hp.p_critical()*0.99),
                                            self.state_hp.p_critical()*0.99)))
            p_hp_2x = res.x
            imin = np.argmin(self.residuals_p)
            if p_hp_2x != self.solutions_p[imin]:
                p_hp_2x = self.solutions_p[imin]
                
        try:    p_hp_2x = p_hp_2x[0]
        except: p_hp_2x = p_hp_2x
        self.resi = self.resi_p(p_hp_2x)
        return p_hp_2x
    
    def resi_p(self,iter_var):
        """
        Residuals of pinch-point method in condenser.  
        Used to find the condensing pressure in the SBVCHP.
        
        Cells discretization: 0.25K per cell (for primary fluid)
        """
        # --- Guess values ----------------------------------------------------
        p_hp_2x = iter_var
        # --- Computation -----------------------------------------------------
        p_hp_2  = p_hp_2x
        if self.parameters['version'] == 'thermodynamic_full' \
        or self.parameters['version'] == 'operational_light':
            self.state_hp.update(CoolProp.PSmass_INPUTS,p_hp_2,self.s_hp_1)
            h_hp_2s  = self.state_hp.hmass()
            h_hp_2   = self.i_hp_1 \
                     + (h_hp_2s-self.i_hp_1)/self.parameters['eta_max_cp']
            self.state_hp.update(CoolProp.HmassP_INPUTS,h_hp_2,p_hp_2)
            T_hp_2   = self.state_hp.T()
        self.state_hp.update(CoolProp.PQ_INPUTS,p_hp_2x,1.0)
        h_hp_2x  = self.state_hp.hmass()
        T_hp_2x  = self.state_hp.T()
        p_hp_3x  = p_hp_2x*(1-self.parameters['dp_hp_cd'])
        self.state_hp.update(CoolProp.PQ_INPUTS,p_hp_3x,0.0)
        T_hp_3x  = self.state_hp.T()
        h_hp_3x  = self.state_hp.hmass()
        T_hp_3   = T_hp_3x-self.parameters['dT_hp_cd_sc']
        p_hp_3   = p_hp_3x
        if  p_hp_3 < self.state_hp.p_critical()\
        and T_hp_3 < self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_liquid)
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_3,T_hp_3)
        h_hp_3   = self.state_hp.hmass()
        self.state_hp.unspecify_phase()
        
        # --> Energy balance on condenser -------------------------------------
        frac_hs = (h_hp_2-h_hp_3)/(self.i_hp_hs_ex-self.i_hp_hs_su)
        T_hi    = T_hp_2x - self.parameters['dT_hp_cd_pp']
        self.state_hs.update(CoolProp.PT_INPUTS,self.p_hp_hs_ex,T_hi)
        h_hi    = self.state_hs.hmass()
        d_hi    = frac_hs*(h_hi-self.i_hp_hs_su)-(h_hp_2x-h_hp_3)
        T_wi    = T_hp_2x
        # .......... check if pinch violated ..................................
        cell_resolution = 1.0 # K
        n_elem = max(int(1/cell_resolution*(self.T_hp_hs_ex-T_hi)+1),3)
        p_hp_wf_vec  = np.linspace(p_hp_2x,        p_hp_2,         n_elem)
        h_hp_wf_vec  = np.linspace(h_hp_2x,        h_hp_2,         n_elem)
        p_hp_hs_vec  = np.linspace(self.p_hp_hs_ex,self.p_hp_hs_ex,n_elem)
        h_hp_hs_vec  = np.linspace(h_hi,           self.i_hp_hs_ex,n_elem)
        T_hp_wf_vec  = np.zeros(n_elem)
        T_hp_hs_vec  = np.zeros(n_elem)
        for i, T in enumerate(T_hp_wf_vec):
            self.state_hp.update(CoolProp.HmassP_INPUTS,h_hp_wf_vec[i],
                                                        p_hp_wf_vec[i])
            self.state_hs.update(CoolProp.HmassP_INPUTS,h_hp_hs_vec[i],
                                                        p_hp_hs_vec[i])
            T_hp_wf_vec[i] = self.state_hp.T()
            T_hp_hs_vec[i] = self.state_hs.T()
        real_pp = round(min(T_hp_wf_vec-T_hp_hs_vec),3)
        trgt_pp = self.parameters['dT_hp_cd_pp']
        chck_pp = real_pp < trgt_pp
        # .......... discretized model ........................................
        if chck_pp:
            self.m_rat      = frac_hs
            cell_resolution = 0.25 # K
            # ======> NEW <====================================================
            n_el_lq         = int(1/cell_resolution*abs((T_hp_3 -T_hp_3x))+2)
            n_el_st         = int(1/cell_resolution*abs((T_hp_3x-T_hp_2x))+2)
            n_el_vp         = int(1/cell_resolution*abs((T_hp_2x-T_hp_2 ))+2)  
            p_hp_wf_lq_vec  = np.linspace(p_hp_3, p_hp_3x,n_el_lq)
            p_hp_wf_st_vec  = np.linspace(p_hp_3x,p_hp_2x,n_el_st)
            p_hp_wf_vp_vec  = np.linspace(p_hp_2x,p_hp_2, n_el_vp)
            p_hp_wf_vec     = np.concatenate((p_hp_wf_lq_vec,
                                              p_hp_wf_st_vec,
                                              p_hp_wf_vp_vec),axis=None)
            h_hp_wf_lq_vec  = np.linspace(h_hp_3, h_hp_3x,n_el_lq)
            h_hp_wf_st_vec  = np.linspace(h_hp_3x,h_hp_2x,n_el_st)
            h_hp_wf_vp_vec  = np.linspace(h_hp_2x,h_hp_2, n_el_vp)
            h_hp_wf_vec     = np.concatenate((h_hp_wf_lq_vec,
                                              h_hp_wf_st_vec,
                                              h_hp_wf_vp_vec),axis=None)
            n_elem          = len(p_hp_wf_vec)
            # --- working fluid -----------------------------------------------
            T_hp_wf_vec     = np.zeros(n_elem)
            for i, h in enumerate(h_hp_wf_vec):
                self.state_hp.update(CoolProp.HmassP_INPUTS,h,p_hp_wf_vec[i])
                T_hp_wf_vec[i] = self.state_hp.T()
            # --- secondary fluid ---------------------------------------------
            p_hp_hs_vec   = np.linspace(self.p_hp_hs_su,self.p_hp_hs_ex,n_elem)
            T_hp_hs_vec   = np.zeros(n_elem)
            T_hp_hs_vec[0]= self.T_hp_hs_su
            h_hp_hs_vec   = np.zeros(n_elem)
            h_hp_hs_vec[0]= self.i_hp_hs_su
            for i in range(len(h_hp_hs_vec)-1):
                h_hp_hs_vec[i+1] = h_hp_hs_vec[i]\
                                  +(h_hp_wf_vec[i+1]-h_hp_wf_vec[i])/self.m_rat
                self.state_hs.update(CoolProp.HmassP_INPUTS,h_hp_hs_vec[i+1],
                                                            p_hp_hs_vec[i+1])
                T_hp_hs_vec[i+1] = self.state_hs.T()
            # --- find pinch --------------------------------------------------
            imin = np.argmin(T_hp_wf_vec-T_hp_hs_vec)
            T_hi = T_hp_hs_vec[imin]
            T_wi = T_hp_wf_vec[imin]
            d_hi = T_hp_wf_vec[imin]-T_hp_hs_vec[imin]\
                                    -self.parameters['dT_hp_cd_pp']

        # --- Final balance ---------------------------------------------------
        out     = d_hi
        self.solutions_p.append(p_hp_2x.copy())
        self.residuals_p.append(d_hi)
        self.T_wi = T_wi
        self.T_hi = T_hi
        return out

    def find_mass_flow(self):
        """
        Find the mass flow rate for the selected operational conditions.
        """
        gue     = 1
        m_hp    = scipy.optimize.fsolve(self.resi_mass_flow, gue)[0]
        if m_hp < 0:    m_hp = 0
        return m_hp
    
    def resi_mass_flow(self,iter_var):
        """
        Residuals of different parameters to control the heat pump.
        """
        # --- Evaluate the massflow -------------------------------------------
        cheat_m_hp       = iter_var
        # --- Check if match the demanded power input -------------------------
        if  self.parameters['mode'] == 'power':
            resi         = self.P_hp - cheat_m_hp*(self.i_hp_2-self.i_hp_1)
        # --- Check if match the maximum heat input ---------------------------
        if  self.parameters['mode'] == 'source':
            q_available  = self.i_hp_cs_su  - self.i_hp_cs_ex
            Q_available  = self.m_hp_cs     * q_available
            q_required   = self.i_hp_1      - self.i_hp_4
            m_required   = Q_available/q_required
            resi         = cheat_m_hp - m_required
        # --- Check if match the maximum heat output --------------------------
        if  self.parameters['mode'] == 'sink':
            q_available  = self.i_hp_hs_ex  - self.i_hp_hs_su
            Q_available  = self.m_hp_hs     * q_available
            q_required   = self.i_hp_2      - self.i_hp_3
            m_required   = Q_available/q_required
            resi         = cheat_m_hp - m_required
        return resi
    
    def find_flow_secondary_hs(self):
        """
        Find the mass flow rate of the hot secondary fluid to match the 
        specified temperature profile based on the refrigerant flow rate.
        """
        q_available  = self.i_hp_2      - self.i_hp_3
        Q_available  = self.m_hp        * q_available
        q_required   = self.i_hp_hs_ex  - self.i_hp_hs_su
        if q_required > 0:  m_required = Q_available/q_required
        else:               m_required = 0
        return m_required

    def find_flow_secondary_cs(self):
        """
        Find the mass flow rate of the cold secondary fluid to match the 
        specified temperature profile based on the refrigerant flow rate.
        """
        q_available  = self.i_hp_1      - self.i_hp_4
        Q_available  = self.m_hp        * q_available
        q_required   = self.i_hp_cs_su  - self.i_hp_cs_ex
        if q_required > 0:  m_required = Q_available/q_required
        else:               m_required = 0
        return m_required
    
    def check_consistency(self):
        """
        Check that the results are consistent.
        """
        if self.s_hp_1x>self.s_hp_1:  self.error=True
        if self.s_hp_2x>self.s_hp_2:  self.error=True
        if self.s_hp_3x<self.s_hp_3:  self.error=True
        if self.s_hp_4 <self.s_hp_3:  self.error=True
        if self.s_hp_3x>self.s_hp_2x: self.error=True
        if self.s_hp_4x>self.s_hp_1x: self.error=True
        if self.x_hp_4 <0:            self.error=True
        if self.eta_hp_cyclen < 0:    self.error=True
        if self.T_hp_4>self.T_hp_cs_ex-0.95*self.parameters['dT_hp_ev_pp']:
            self.error=True
        if self.T_hp_1>self.T_hp_cs_su-0.95*self.parameters['dT_hp_ev_pp']:
            self.error=True
        if self.T_hp_2<self.T_hp_hs_ex+0.95*self.parameters['dT_hp_cd_pp']:
            self.error=True
        if self.T_hp_3<self.T_hp_hs_su+0.95*self.parameters['dT_hp_cd_pp']:
            self.error=True
        if self.T_wi  <self.T_hi      +0.95*self.parameters['dT_hp_cd_pp']:
            self.error = True   
        if abs(self.resi) > 1:              self.error = True
        if self.p_hp_2x < self.p_hp_2x_min: self.error = True
    
    def retrieve_kpi(self):
        """
        Retrieve the KPI's of the heat pump.
        """
        
        if self.options['exergy']: self.compute_exergy()
        
        # --- Basic quantities -------------------------------------------------
        self.w_in_hp_en = self.i_hp_2-self.i_hp_1
        self.q_in_hp_en = self.i_hp_1-self.i_hp_4
        self.w_out_hp_en= self.i_hp_3-self.i_hp_4 
        self.q_out_hp_en= self.i_hp_2-self.i_hp_3
        self.f_hp_hs    = self.q_out_hp_en/(self.i_hp_hs_ex-self.i_hp_hs_su)
        self.f_hp_cs    = self.q_in_hp_en /(self.i_hp_cs_su-self.i_hp_cs_ex)
        
        if self.options['exergy']: 
            self.w_in_hp_ex = self.e_hp_2-self.e_hp_1 
            self.q_in_hp_ex = self.e_hp_1-self.e_hp_4
            self.w_out_hp_ex= self.e_hp_3-self.e_hp_4 
            self.q_out_hp_ex= self.e_hp_2-self.e_hp_3 
        
        # --- Losses ----------------------------------------------------------
            self.loss_hp_compex = self.w_in_hp_en-self.w_in_hp_ex
            self.loss_hp_valvex = self.w_out_hp_ex-self.w_out_hp_en
            self.loss_hp_evapex = self.f_hp_cs*(self.e_hp_cs_su\
                                - self.e_hp_cs_ex)-self.q_in_hp_ex
            self.loss_hp_srcen  = self.f_hp_cs*(self.i_hp_cs_ex\
                                - self.i_hp_cs_ref)
            self.loss_hp_srcex  = self.f_hp_cs*(self.e_hp_cs_ex\
                                - self.e_hp_cs_ref)
            self.loss_hp_condex = self.q_out_hp_ex-self.f_hp_hs\
                                *(self.e_hp_hs_ex-self.e_hp_hs_su)
            self.pow_hp_sinken  = self.f_hp_hs*(self.i_hp_hs_ex\
                                - self.i_hp_hs_su)
            self.pow_hp_sinkex  = self.f_hp_hs*(self.e_hp_hs_ex\
                                - self.e_hp_hs_su)
            self.pow_hp_supplen = self.f_hp_cs*(self.i_hp_cs_su\
                                - self.i_hp_cs_ref)+self.w_in_hp_en
            self.pow_hp_supplex = self.f_hp_cs*(self.e_hp_cs_su\
                                - self.e_hp_cs_ref)+self.w_in_hp_en
        
        # if not exergy mode, still need supplex for exergy efficiency
        if not self.options['exergy']:
            
            self.compute_ref_states()
            self.compute_exergy_boundaries()
              
            self.pow_hp_supplex = self.f_hp_cs*(self.e_hp_cs_su\
                                - self.e_hp_cs_ref)+self.w_in_hp_en
        
        # --- Efficiencies ----------------------------------------------------
        self.eta_hp_cyclen = self.q_out_hp_en/self.w_in_hp_en
        
        if self.options['exergy']:
            self.eta_hp_cyclex = self.q_out_hp_ex/self.w_in_hp_en
            self.eta_hp_compex = self.w_in_hp_ex /self.w_in_hp_en
            self.eta_hp_valvex = self.w_out_hp_en/self.w_out_hp_ex
            self.eta_hp_evapen = 1
            self.eta_hp_evapex = self.f_hp_cs**(-1)*self.q_in_hp_ex\
                               /(self.e_hp_cs_su-self.e_hp_cs_ex)
            if self.i_hp_cs_su == self.i_hp_cs_ref: 
                self.eta_hp_srcen  = float('inf')
            else:                                   
                self.eta_hp_srcen  = (self.i_hp_cs_su-self.i_hp_cs_ex)\
                                    /(self.i_hp_cs_su-self.i_hp_cs_ref)
            if self.e_hp_cs_su == self.e_hp_cs_ref: 
                self.eta_hp_srcex  = float('inf')
            else:                                   
                self.eta_hp_srcex  = (self.e_hp_cs_su-self.e_hp_cs_ex)\
                                    /(self.e_hp_cs_su-self.e_hp_cs_ref)
            self.eta_hp_genven = self.eta_hp_evapen*self.eta_hp_srcen
            self.eta_hp_genvex = self.eta_hp_evapex*self.eta_hp_srcex  
            self.eta_hp_condex = self.f_hp_hs*(self.e_hp_hs_ex\
                                - self.e_hp_hs_su)/self.q_out_hp_ex
            self.eta_hp_toten  = self.pow_hp_sinken/self.pow_hp_supplen
            self.eta_hp_totex  = self.pow_hp_sinkex/self.pow_hp_supplex
        
        # --- Carnot and Lorentz efficiencies ---------------------------------
        T_cs_logAv =(self.T_hp_cs_su-self.T_hp_cs_ex)\
             /np.log(self.T_hp_cs_su/self.T_hp_cs_ex)
        T_hs_logAv =(self.T_hp_hs_ex-self.T_hp_hs_su)\
             /np.log(self.T_hp_hs_ex/self.T_hp_hs_su)
        COP_lorenz =      T_hs_logAv/(     T_hs_logAv-     T_cs_logAv)
        COP_carnot = self.T_hp_hs_ex/(self.T_hp_hs_ex-self.T_hp_cs_su)
        self.cop_hp_lorenz = COP_lorenz
        self.cop_hp_carnot = COP_carnot
        self.psi_hp_lorenz = self.eta_hp_cyclen/COP_lorenz
        self.psi_hp_carnot = self.eta_hp_cyclen/COP_carnot
        
        # --- Preliminary design information ----------------------------------
        self.vol_coef     = self.v_hp_2/self.w_in_hp_en
        self.comp_ratio   = self.p_hp_2/self.p_hp_1
        self.vol_ratio    = self.v_hp_1/self.v_hp_2
        if self.parameters['version'] == 'operational_light':
            # --- Multiply the specific quantities by the mass flow -----------
            self.Q_in_hp_en = self.m_hp*self.q_in_hp_en
            self.W_in_hp_en = self.m_hp*self.w_in_hp_en
            self.Q_out_hp_en= self.m_hp*self.q_out_hp_en
            
    def export_states(self):
        """
        Export the thermodynamic states of the SBVCHP necessary to draw
        the T-s and p-h diagrams.
        """
        self.p_hp = self.p_hp_1,self.p_hp_1x,self.p_hp_2,self.p_hp_2x,\
                    self.p_hp_3,self.p_hp_3x,self.p_hp_4,self.p_hp_4x
        self.T_hp = self.T_hp_1,self.T_hp_1x,self.T_hp_2,self.T_hp_2x,\
                    self.T_hp_3,self.T_hp_3x,self.T_hp_4,self.T_hp_4x
        self.i_hp = self.i_hp_1,self.i_hp_1x,self.i_hp_2,self.i_hp_2x,\
                    self.i_hp_3,self.i_hp_3x,self.i_hp_4,self.i_hp_4x
        self.s_hp = self.s_hp_1,self.s_hp_1x,self.s_hp_2,self.s_hp_2x,\
                    self.s_hp_3,self.s_hp_3x,self.s_hp_4,self.s_hp_4x
        self.x_hp = self.x_hp_1,self.x_hp_1x,self.x_hp_2,self.x_hp_2x,\
                    self.x_hp_3,self.x_hp_3x,self.x_hp_4,self.x_hp_4x
        if self.options['exergy']:
            self.e_hp = self.e_hp_1,self.e_hp_1x,self.e_hp_2,self.e_hp_2x,\
                        self.e_hp_3,self.e_hp_3x,self.e_hp_4,self.e_hp_4x
        else:
            self.e_hp = 0,0,0,0,0,0,0,0
        out = self.p_hp,self.T_hp,self.i_hp,self.s_hp,self.x_hp,self.e_hp
        return out
    
    def compute_ref_states(self):
        """
        Compute the reference state.
        """
        self.T_hp_0 = self.parameters['T_0']
        self.p_hp_0 = self.parameters['p_0']
        self.T_ref  = self.parameters['T_ref']
        self.p_ref  = self.parameters['p_ref']
        
        try:    # Saturated liquid conditions must be taken as reference if the 
                # fluid can condense at T_0: the exergy is therefore 0 for any 
                # state of saturated vapour at T_0. This is because the exergy 
                # variation associated with a phase change at T_0 is always 0.
                # Refer to "THERMAL POWER PLANTS: Energetic and exergetic 
                # approaches" pp.13 for more details.
            self.state_hp.update(CoolProp.QT_INPUTS,0,self.T_hp_0)
            self.i_hp_0 = self.state_hp.hmass()
            self.s_hp_0 = self.state_hp.smass()
            if np.isnan(self.i_hp_0)\
            or np.isnan(self.s_hp_0): 
                raise TypeError('Reference state is NaN!')
        except: # else the fluid above the critical point
            self.state_hp.update(CoolProp.PT_INPUTS,self.p_hp_0,self.T_hp_0)
            self.i_hp_0 = self.state_hp.hmass()
            self.s_hp_0 = self.state_hp.smass()
        
        try:    # Saturated liquid conditions must be taken as reference if the 
                # fluid can condense at T_0: the exergy is therefore 0 for any 
                # state of saturated vapour at T_0. This is because the exergy 
                # variation associated with a phase change at T_0 is always 0.
                # Refer to "THERMAL POWER PLANTS: Energetic and exergetic 
                # approaches" pp.13 for more details.
            self.state_cs.update(CoolProp.QT_INPUTS,0,self.T_hp_0)
            self.i_hp_cs_0 = self.state_cs.hmass()
            self.s_hp_cs_0 = self.state_cs.smass()
            if np.isnan(self.i_hp_cs_0)\
            or np.isnan(self.s_hp_cs_0): 
                raise TypeError('Reference state is NaN!')
        except: # else the fluid above the critical point
            self.state_cs.update(CoolProp.PT_INPUTS,self.p_hp_0,self.T_hp_0)
            self.i_hp_cs_0 = self.state_cs.hmass()
            self.s_hp_cs_0 = self.state_cs.smass()
            
        if self.options['exergy']:
            try:# Saturated liquid conditions must be taken as reference if the 
                # fluid can condense at T_0: the exergy is therefore 0 for any 
                # state of saturated vapour at T_0. This is because the exergy 
                # variation associated with a phase change at T_0 is always 0.
                # Refer to "THERMAL POWER PLANTS: Energetic and exergetic 
                # approaches" pp.13 for more details.
                self.state_hs.update(CoolProp.QT_INPUTS,0,self.T_hp_0)
                self.i_hp_hs_0 = self.state_hs.hmass()
                self.s_hp_hs_0 = self.state_hs.smass()
                if np.isnan(self.i_hp_hs_0)\
                or np.isnan(self.s_hp_hs_0): 
                    raise TypeError('Reference state is NaN!')
            except: # else the fluid above the critical point
                self.state_hs.update(CoolProp.PT_INPUTS,self.p_hp_0,
                                                        self.T_hp_0)
                self.i_hp_hs_0 = self.state_hs.hmass()
                self.s_hp_hs_0 = self.state_hs.smass()
    
    def compute_exergy(self):
        """
        Compute the exergy at each state of the SBVCHP.
        """
        self.compute_ref_states()
        
        # - 1 - Exergy of the cycle -------------------------------------------
        
        self.e_hp_1  = (self.i_hp_1 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_1 -self.s_hp_0)
        self.e_hp_1x = (self.i_hp_1x-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_1x-self.s_hp_0)
        self.e_hp_2  = (self.i_hp_2 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_2 -self.s_hp_0)
        self.e_hp_2x = (self.i_hp_2x-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_2x-self.s_hp_0)
        self.e_hp_3  = (self.i_hp_3 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_3 -self.s_hp_0)
        self.e_hp_3x = (self.i_hp_3x-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_3x-self.s_hp_0)
        self.e_hp_4  = (self.i_hp_4 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_4 -self.s_hp_0)
        self.e_hp_4x = (self.i_hp_4x-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_4x-self.s_hp_0)
        
        if round(self.T_hp_1 ,2)==round(self.T_hp_0,2): self.e_hp_1  = 0
        if round(self.T_hp_1x,2)==round(self.T_hp_0,2): self.e_hp_1x = 0
        if round(self.T_hp_2 ,2)==round(self.T_hp_0,2): self.e_hp_2  = 0
        if round(self.T_hp_2x,2)==round(self.T_hp_0,2): self.e_hp_2x = 0
        if round(self.T_hp_3 ,2)==round(self.T_hp_0,2): self.e_hp_3  = 0
        if round(self.T_hp_3x,2)==round(self.T_hp_0,2): self.e_hp_3x = 0
        if round(self.T_hp_4 ,2)==round(self.T_hp_0,2): self.e_hp_4  = 0
        if round(self.T_hp_4x,2)==round(self.T_hp_0,2): self.e_hp_4x = 0
        
        # 2 - Exergy of the boundaries ----------------------------------------
        
        self.compute_exergy_boundaries()
        
    def compute_exergy_boundaries(self):
        """
        Compute the exergy of the boundaries.
        """
        
        # 1 - Exergy of the source boundary -----------------------------------

        self.e_hp_cs_su = (self.i_hp_cs_su -self.i_hp_cs_0)\
             -self.T_hp_0*(self.s_hp_cs_su -self.s_hp_cs_0)
        self.e_hp_cs_ex = (self.i_hp_cs_ex -self.i_hp_cs_0)\
             -self.T_hp_0*(self.s_hp_cs_ex -self.s_hp_cs_0)
        
        if round(self.T_hp_cs_su,2)==round(self.T_hp_0,2): self.e_hp_cs_su = 0
        if round(self.T_hp_cs_ex,2)==round(self.T_hp_0,2): self.e_hp_cs_ex = 0
        
        self.state_cs.update(CoolProp.PT_INPUTS,self.p_ref,self.T_ref)
        self.i_hp_cs_ref= self.state_cs.hmass()
        self.s_hp_cs_ref= self.state_cs.smass()
        self.e_hp_cs_ref= (self.i_hp_cs_ref -self.i_hp_cs_0)\
             -self.T_hp_0*(self.s_hp_cs_ref -self.s_hp_cs_0)
        
        if round(self.T_ref     ,2)==round(self.T_hp_0,2): self.e_hp_cs_ref= 0
        
        # 2 - Exergy of the sink boundary -------------------------------------
        
        if self.options['exergy']:
        
            self.e_hp_hs_su =(self.i_hp_hs_su -self.i_hp_hs_0)\
                -self.T_hp_0*(self.s_hp_hs_su -self.s_hp_hs_0)
            self.e_hp_hs_ex =(self.i_hp_hs_ex -self.i_hp_hs_0)\
                -self.T_hp_0*(self.s_hp_hs_ex -self.s_hp_hs_0)
            
            if round(self.T_hp_hs_su,2)==round(self.T_hp_0,2): 
                self.e_hp_hs_su = 0
            if round(self.T_hp_hs_ex,2)==round(self.T_hp_0,2): 
                self.e_hp_hs_ex = 0
            
            self.state_hs.update(CoolProp.PT_INPUTS,self.p_ref,self.T_ref)
            self.i_hp_hs_ref= self.state_hs.hmass()
            self.s_hp_hs_ref= self.state_hs.smass()
            self.e_hp_hs_ref= (self.i_hp_hs_ref -self.i_hp_hs_0)\
                 -self.T_hp_0*(self.s_hp_hs_ref -self.s_hp_hs_0)
                
            if round(self.T_ref     ,2)==round(self.T_hp_0,2): 
                self.e_hp_hs_ref= 0

#%% MODEL: SRVCHP

class SRVCHP(SBVCHP):
    """
    Class for the simulation of Subcritical Recuperated Vapor Compression Heat 
    Pumps.

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (heat pump boundaries and operational 
                                        parameters).
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
    # ==== HEAT PUMP MODEL ================================================== #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the cold source
            * 'sink':   available heat flow rate at the hot sink
            * 'power':  electrical power supplied to the compressor
        
        """        
        
        # --- Cycle computation -----------------------------------------------
        
        # ..... two-phase .....................................................
        self.T_hp_4x = min(  self.T_hp_cs_su
                           - self.parameters['dT_hp_ev_pp']
                           - self.parameters['dT_hp_ev_sh'],
                             self.T_hp_cs_ex
                           - self.parameters['dT_hp_ev_pp'])
        if self.T_hp_4x > self.state_hp.T_critical():
            raise ValueError('An inconsistency was detected in the cycle!\
                              In: '+os.path.join(os.path.abspath(__file__),
                              'SRVCHP','evaluate_cycle:\
                    the evaporation temperature is above the critical point!'))
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.QT_INPUTS,0.5,self.T_hp_4x)
        self.p_hp_4x = self.state_hp.p()
        self.state_hp.unspecify_phase()
        # ..... saturated (gas) ...............................................
        self.p_hp_1x = self.p_hp_4x*(1-self.parameters['dp_hp_ev']) 
        self.x_hp_1x = 1.0
        self.state_hp.update(CoolProp.PQ_INPUTS,self.p_hp_1x,self.x_hp_1x)
        self.T_hp_1x = self.state_hp.T()
        self.i_hp_1x = self.state_hp.hmass()
        self.s_hp_1x = self.state_hp.smass()
        # ..... gas ...........................................................
        self.T_hp_1  = self.T_hp_1x+self.parameters['dT_hp_ev_sh']
        self.p_hp_1  = self.p_hp_1x
        if  self.p_hp_1 < self.state_hp.p_critical()\
        and self.T_hp_1 < self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_gas)
        if  self.p_hp_1 < self.state_hp.p_critical()\
        and self.T_hp_1 > self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_hp.update(CoolProp.PT_INPUTS,self.p_hp_1,self.T_hp_1)
        self.i_hp_1  = self.state_hp.hmass()
        self.x_hp_1  = self.state_hp.Q()
        self.s_hp_1  = self.state_hp.smass()
        self.v_hp_1  = self.state_hp.rhomass()**(-1)
        self.state_hp.unspecify_phase()
        # ..... saturated (gas) ...............................................
        self.p_hp_2x = self.find_p()
        self.x_hp_2x = 1.0
        self.state_hp.update(CoolProp.PQ_INPUTS,self.p_hp_2x,self.x_hp_2x)
        self.T_hp_2x = self.state_hp.T()
        self.i_hp_2x = self.state_hp.hmass()       
        self.s_hp_2x = self.state_hp.smass()
        # ..... saturated (liquid) ............................................
        self.p_hp_3x = self.p_hp_2x*(1-self.parameters['dp_hp_cd'])
        self.x_hp_3x = 0.0
        self.state_hp.update(CoolProp.PQ_INPUTS,self.p_hp_3x,self.x_hp_3x)
        self.T_hp_3x = self.state_hp.T()
        self.i_hp_3x = self.state_hp.hmass()
        self.s_hp_3x = self.state_hp.smass()
        # ..... liquid ........................................................
        self.T_hp_3  = self.T_hp_3x-self.parameters['dT_hp_cd_sc']
        self.p_hp_3  = self.p_hp_3x
        self.state_hp.specify_phase(CoolProp.iphase_liquid)
        self.state_hp.update(CoolProp.PT_INPUTS,self.p_hp_3,self.T_hp_3)
        self.i_hp_3  = self.state_hp.hmass()
        self.s_hp_3  = self.state_hp.smass()
        self.x_hp_3  = self.state_hp.Q()
        self.state_hp.unspecify_phase()
        # ..... gas ...........................................................
        self.T_hp_1r  = self.T_hp_1\
                      +(self.T_hp_3-self.T_hp_1)*self.parameters['epsilon']
        self.p_hp_1r  = self.p_hp_1*(1-self.parameters['dp_hp_rg_vp'])
        self.state_hp.update(CoolProp.PT_INPUTS,self.p_hp_1r,self.T_hp_1r)
        self.s_hp_1r = self.state_hp.smass()
        self.i_hp_1r = self.state_hp.hmass()
        self.x_hp_1r = self.state_hp.Q()
        # ..... gas ...........................................................
        self.p_hp_2  = self.p_hp_2x
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_hp.update(CoolProp.PSmass_INPUTS,
                                    self.p_hp_2,self.s_hp_1r)
            self.i_hp_2s = self.state_hp.hmass()
            self.i_hp_2  = self.i_hp_1r \
                         +(self.i_hp_2s-self.i_hp_1r)\
                         / self.parameters['eta_max_cp']
            self.m_hp    = 0
            self.P_cp    = 0
            self.eta_is_cp = self.parameters['eta_max_cp']
            self.cp_flag = 'na'
            self.state_hp.update(CoolProp.HmassP_INPUTS,
                                    self.i_hp_2,self.p_hp_2)
            self.T_hp_2  = self.state_hp.T()
        
        self.s_hp_2  = self.state_hp.smass()
        self.x_hp_2  = self.state_hp.Q()
        self.v_hp_2  = self.state_hp.rhomass()**(-1)
        # ..... liquid ........................................................
        self.i_hp_3r = self.i_hp_3-(self.i_hp_1r-self.i_hp_1)
        self.p_hp_3r = self.p_hp_3*(1-self.parameters['dp_hp_rg_lq'])
        self.state_hp.specify_phase(CoolProp.iphase_liquid)
        self.state_hp.update(CoolProp.HmassP_INPUTS,self.i_hp_3r,self.p_hp_3r)
        self.T_hp_3r = self.state_hp.T()
        self.s_hp_3r = self.state_hp.smass()
        self.x_hp_3r = self.state_hp.Q()
        self.state_hp.unspecify_phase()
        # ..... two-phase .....................................................
        self.p_hp_4  = self.p_hp_4x 
        self.i_hp_4  = self.i_hp_3r
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.HmassP_INPUTS,self.i_hp_4,self.p_hp_4)
        self.T_hp_4  = self.state_hp.T()
        self.s_hp_4  = self.state_hp.smass()
        self.x_hp_4  = self.state_hp.Q()
        self.state_hp.unspecify_phase()
        self.i_hp_4x = self.i_hp_4
        self.T_hp_4x = self.T_hp_4
        self.s_hp_4x = self.s_hp_4
        self.x_hp_4x = self.x_hp_4
        
        # --- Evaluate the operational conditions -----------------------------
        
        if self.parameters['version'] == 'operational_light':
            self.m_hp    = self.find_mass_flow()
            self.P_cp    = self.m_hp    * (self.i_hp_2 - self.i_hp_1r)
            self.m_hp_cs = self.find_flow_secondary_cs()
            self.m_hp_hs = self.find_flow_secondary_hs()
        
        # --- Evaluate the Key Performance Indicators -------------------------
        
        self.retrieve_kpi()
        
    # ======================================================================= #
    # ==== SUB-FUNCTIONS FOR HEAT PUMP MODEL ================================ #
    # ======================================================================= #
    
    def resi_p(self,iter_var):
        """
        Residuals of pinch-point method in condenser.  
        Used to find the condensing pressure in the SRVCHP.
        
        Cells discretization: 0.25K per cell (for primary fluid)
        """
        # --- Guess values ----------------------------------------------------
        p_hp_2x  = iter_var
        # --- Computation -----------------------------------------------------
        p_hp_2   = p_hp_2x
        self.state_hp.update(CoolProp.PQ_INPUTS,p_hp_2x,1.0)
        h_hp_2x  = self.state_hp.hmass()
        T_hp_2x  = self.state_hp.T()
        p_hp_3x  = p_hp_2x*(1-self.parameters['dp_hp_cd'])
        self.state_hp.update(CoolProp.PQ_INPUTS,p_hp_3x,0.0)
        T_hp_3x  = self.state_hp.T()
        h_hp_3x  = self.state_hp.hmass()
        T_hp_3   = T_hp_3x-self.parameters['dT_hp_cd_sc']
        p_hp_3   = p_hp_3x
        if  p_hp_3 < self.state_hp.p_critical()\
        and T_hp_3 < self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_liquid)
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_3,T_hp_3)
        h_hp_3   = self.state_hp.hmass()
        self.state_hp.unspecify_phase()
        T_hp_1r  = self.T_hp_1+(T_hp_3-self.T_hp_1)*self.parameters['epsilon']
        p_hp_1r  = self.p_hp_1*(1-self.parameters['dp_hp_rg_vp'])
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_1r,T_hp_1r)
        s_hp_1r  = self.state_hp.smass()
        h_hp_1r  = self.state_hp.hmass()

        if self.parameters['version'] == 'thermodynamic_full' \
        or self.parameters['version'] == 'operational_light':
            self.state_hp.update(CoolProp.PSmass_INPUTS,p_hp_2,s_hp_1r)
            h_hp_2s  = self.state_hp.hmass()
            h_hp_2   =  h_hp_1r \
                     + (h_hp_2s-h_hp_1r)/self.parameters['eta_max_cp']
            self.state_hp.update(CoolProp.HmassP_INPUTS,h_hp_2,p_hp_2)
            T_hp_2   = self.state_hp.T()
        
        # --> Energy balance on condenser -------------------------------------
        frac_hs = (h_hp_2-h_hp_3)/(self.i_hp_hs_ex-self.i_hp_hs_su)
        T_hi    = T_hp_2x - self.parameters['dT_hp_cd_pp']
        self.state_hs.update(CoolProp.PT_INPUTS,self.p_hp_hs_ex,T_hi)
        h_hi    = self.state_hs.hmass()
        d_hi    = frac_hs*(h_hi-self.i_hp_hs_su)-(h_hp_2x-h_hp_3)
        T_wi    = T_hp_2x
        # .......... check if pinch violated ..................................
        cell_resolution = 1.0 # K
        n_elem = max(int(1/cell_resolution*(self.T_hp_hs_ex-T_hi)+1),3)
        p_hp_wf_vec  = np.linspace(p_hp_2x,        p_hp_2,         n_elem)
        h_hp_wf_vec  = np.linspace(h_hp_2x,        h_hp_2,         n_elem)
        p_hp_hs_vec  = np.linspace(self.p_hp_hs_ex,self.p_hp_hs_ex,n_elem)
        h_hp_hs_vec  = np.linspace(h_hi,           self.i_hp_hs_ex,n_elem)
        T_hp_wf_vec  = np.zeros(n_elem)
        T_hp_hs_vec  = np.zeros(n_elem)
        for i, T in enumerate(T_hp_wf_vec):
            self.state_hp.update(CoolProp.HmassP_INPUTS,h_hp_wf_vec[i],
                                                        p_hp_wf_vec[i])
            self.state_hs.update(CoolProp.HmassP_INPUTS,h_hp_hs_vec[i],
                                                        p_hp_hs_vec[i])
            T_hp_wf_vec[i] = self.state_hp.T()
            T_hp_hs_vec[i] = self.state_hs.T()
        real_pp = round(min(T_hp_wf_vec-T_hp_hs_vec),3)
        trgt_pp = self.parameters['dT_hp_cd_pp']
        chck_pp = real_pp < trgt_pp
        # .......... discretized model ........................................
        if chck_pp:
            self.m_rat      = frac_hs
            cell_resolution = 0.25 # K
            # ======> NEW <====================================================
            n_el_lq         = int(1/cell_resolution*abs((T_hp_3 -T_hp_3x))+2)
            n_el_st         = int(1/cell_resolution*abs((T_hp_3x-T_hp_2x))+2)
            n_el_vp         = int(1/cell_resolution*abs((T_hp_2x-T_hp_2 ))+2)  
            p_hp_wf_lq_vec  = np.linspace(p_hp_3, p_hp_3x,n_el_lq)
            p_hp_wf_st_vec  = np.linspace(p_hp_3x,p_hp_2x,n_el_st)
            p_hp_wf_vp_vec  = np.linspace(p_hp_2x,p_hp_2, n_el_vp)
            p_hp_wf_vec     = np.concatenate((p_hp_wf_lq_vec,
                                              p_hp_wf_st_vec,
                                              p_hp_wf_vp_vec),axis=None)
            h_hp_wf_lq_vec  = np.linspace(h_hp_3, h_hp_3x,n_el_lq)
            h_hp_wf_st_vec  = np.linspace(h_hp_3x,h_hp_2x,n_el_st)
            h_hp_wf_vp_vec  = np.linspace(h_hp_2x,h_hp_2, n_el_vp)
            h_hp_wf_vec     = np.concatenate((h_hp_wf_lq_vec,
                                              h_hp_wf_st_vec,
                                              h_hp_wf_vp_vec),axis=None)
            n_elem          = len(p_hp_wf_vec)
            # --- working fluid -----------------------------------------------
            T_hp_wf_vec     = np.zeros(n_elem)
            for i, h in enumerate(h_hp_wf_vec):
                self.state_hp.update(CoolProp.HmassP_INPUTS,h,p_hp_wf_vec[i])
                T_hp_wf_vec[i] = self.state_hp.T()
            # --- secondary fluid ---------------------------------------------
            p_hp_hs_vec   = np.linspace(self.p_hp_hs_su,self.p_hp_hs_ex,n_elem)
            T_hp_hs_vec   = np.zeros(n_elem)
            T_hp_hs_vec[0]= self.T_hp_hs_su
            h_hp_hs_vec   = np.zeros(n_elem)
            h_hp_hs_vec[0]= self.i_hp_hs_su
            for i in range(len(h_hp_hs_vec)-1):
                h_hp_hs_vec[i+1] = h_hp_hs_vec[i]\
                                  +(h_hp_wf_vec[i+1]-h_hp_wf_vec[i])/self.m_rat
                self.state_hs.update(CoolProp.HmassP_INPUTS,h_hp_hs_vec[i+1],
                                                            p_hp_hs_vec[i+1])
                T_hp_hs_vec[i+1] = self.state_hs.T()
            # --- find pinch --------------------------------------------------
            imin = np.argmin(T_hp_wf_vec-T_hp_hs_vec)
            T_hi = T_hp_hs_vec[imin]
            T_wi = T_hp_wf_vec[imin]
            d_hi = T_hp_wf_vec[imin]-T_hp_hs_vec[imin]\
                                    -self.parameters['dT_hp_cd_pp']
        # --- Final balance ---------------------------------------------------
        out     = d_hi
        self.solutions_p.append(p_hp_2x.copy())
        self.residuals_p.append(d_hi)
        self.T_wi = T_wi
        self.T_hi = T_hi
        return out
    
    def resi_mass_flow(self,iter_var):
        """
        Residuals of different parameters to control the heat pump.
        """
        # --- Evaluate the massflow -------------------------------------------
        cheat_m_hp       = iter_var
        # --- Check if match the demanded power input -------------------------
        if  self.parameters['mode'] == 'power':
            resi         = self.P_hp - cheat_m_hp*(self.i_hp_2-self.i_hp_1r)
        # --- Check if match the maximum heat input ---------------------------
        if  self.parameters['mode'] == 'source':
            q_available  = self.i_hp_cs_su  - self.i_hp_cs_ex
            Q_available  = self.m_hp_cs     * q_available
            q_required   = self.i_hp_1      - self.i_hp_4
            m_required   = Q_available/q_required
            resi         = cheat_m_hp - m_required
        # --- Check if match the maximum heat output --------------------------
        if  self.parameters['mode'] == 'sink':
            q_available  = self.i_hp_hs_ex  - self.i_hp_hs_su
            Q_available  = self.m_hp_hs     * q_available
            q_required   = self.i_hp_2      - self.i_hp_3
            m_required   = Q_available/q_required
            resi         = cheat_m_hp - m_required
        return resi
    
    def check_consistency(self):
        """
        Check that the results are consistent.
        """
        if self.s_hp_1x>self.s_hp_1:  self.error=True
        if self.s_hp_1 >self.s_hp_1r: self.error=True
        if self.s_hp_2x>self.s_hp_2:  self.error=True
        if self.s_hp_3x<self.s_hp_3:  self.error=True
        if self.s_hp_3 <self.s_hp_3r: self.error=True
        if self.s_hp_4 <self.s_hp_3r: self.error=True
        if self.s_hp_3x>self.s_hp_2x: self.error=True
        if self.s_hp_4x>self.s_hp_1x: self.error=True
        if self.x_hp_4 <0:            self.error=True
        if self.eta_hp_cyclen < 0:    self.error=True
        if self.T_hp_4>self.T_hp_cs_ex-0.95*self.parameters['dT_hp_ev_pp']:
            self.error=True
        if self.T_hp_1>self.T_hp_cs_su-0.95*self.parameters['dT_hp_ev_pp']:
            self.error=True
        if self.T_hp_2<self.T_hp_hs_ex+0.95*self.parameters['dT_hp_cd_pp']:
            self.error=True
        if self.T_hp_3<self.T_hp_hs_su+0.95*self.parameters['dT_hp_cd_pp']:
            self.error=True
        if self.T_wi  <self.T_hi      +0.95*self.parameters['dT_hp_cd_pp']:
            self.error = True  
        if self.T_hp_1r>self.T_hp_3   -0.95*self.parameters['dT_hp_cd_pp']:#!!!
            self.error = True
            
        if abs(self.resi) > 1:              self.error = True
        if self.p_hp_2x < self.p_hp_2x_min: self.error = True 
    
    def retrieve_kpi(self):
        """
        Retrieve the KPI's of the heat pump.
        """
        
        if self.options['exergy']: self.compute_exergy()
        
        # --- Basic quantities ------------------------------------------------
        self.w_in_hp_en = self.i_hp_2 -self.i_hp_1r
        self.q_in_hp_en = self.i_hp_1 -self.i_hp_4
        self.w_out_hp_en= self.i_hp_3r-self.i_hp_4
        self.q_out_hp_en= self.i_hp_2 -self.i_hp_3
        self.f_hp_hs    = self.q_out_hp_en/(self.i_hp_hs_ex-self.i_hp_hs_su)
        self.f_hp_cs    = self.q_in_hp_en /(self.i_hp_cs_su-self.i_hp_cs_ex)
        
        if self.options['exergy']: 
            self.w_in_hp_ex = self.e_hp_2 -self.e_hp_1r
            self.q_in_hp_ex = self.e_hp_1 -self.e_hp_4
            self.w_out_hp_ex= self.e_hp_3r-self.e_hp_4 
            self.q_out_hp_ex= self.e_hp_2 -self.e_hp_3 
        
        # --- Losses ----------------------------------------------------------
        if self.options['exergy']:
            self.loss_hp_compex = self.w_in_hp_en-self.w_in_hp_ex
            self.loss_hp_valvex = self.w_out_hp_ex-self.w_out_hp_en
            self.loss_hp_regex  =(self.e_hp_3 -self.e_hp_3r)\
                                -(self.e_hp_1r-self.e_hp_1 )
            self.loss_hp_evapex = self.f_hp_cs*(self.e_hp_cs_su\
                                - self.e_hp_cs_ex)-self.q_in_hp_ex
            self.loss_hp_srcen  = self.f_hp_cs*(self.i_hp_cs_ex\
                                - self.i_hp_cs_ref)
            self.loss_hp_srcex  = self.f_hp_cs*(self.e_hp_cs_ex\
                                - self.e_hp_cs_ref)
            self.loss_hp_condex = self.q_out_hp_ex-self.f_hp_hs\
                                *(self.e_hp_hs_ex-self.e_hp_hs_su)
            self.pow_hp_sinken  = self.f_hp_hs*(self.i_hp_hs_ex\
                                - self.i_hp_hs_su)
            self.pow_hp_sinkex  = self.f_hp_hs*(self.e_hp_hs_ex\
                                - self.e_hp_hs_su)
            self.pow_hp_supplen = self.f_hp_cs*(self.i_hp_cs_su\
                                - self.i_hp_cs_ref)+self.w_in_hp_en
            self.pow_hp_supplex = self.f_hp_cs*(self.e_hp_cs_su\
                                - self.e_hp_cs_ref)+self.w_in_hp_en
        
        # if not exergy mode, still need supplex for exergy efficiency
        if not self.options['exergy']:
            self.compute_ref_states()
            self.compute_exergy_boundaries()
                
            self.pow_hp_supplex = self.f_hp_cs*(self.e_hp_cs_su\
                                - self.e_hp_cs_ref)+self.w_in_hp_en
        
        # --- Efficiencies ----------------------------------------------------
        self.eta_hp_cyclen = self.q_out_hp_en/self.w_in_hp_en
        
        if self.options['exergy']:
            self.eta_hp_cyclex = self.q_out_hp_ex/self.w_in_hp_en
            self.eta_hp_compex = self.w_in_hp_ex /self.w_in_hp_en
            self.eta_hp_valvex = self.w_out_hp_en/self.w_out_hp_ex
            self.eta_hp_regex  =(self.e_hp_1r-self.e_hp_1 )\
                               /(self.e_hp_3 -self.e_hp_3r)
            self.eta_hp_evapen = 1
            self.eta_hp_evapex = self.f_hp_cs**(-1)*self.q_in_hp_ex\
                               /(self.e_hp_cs_su-self.e_hp_cs_ex)
            if self.i_hp_cs_su == self.i_hp_cs_ref: 
                self.eta_hp_srcen  = float('inf')
            else:                                   
                self.eta_hp_srcen  = (self.i_hp_cs_su-self.i_hp_cs_ex)\
                                    /(self.i_hp_cs_su-self.i_hp_cs_ref)
            if self.e_hp_cs_su == self.e_hp_cs_ref: 
                self.eta_hp_srcex  = float('inf')
            else:                                   
                self.eta_hp_srcex  = (self.e_hp_cs_su-self.e_hp_cs_ex)\
                                    /(self.e_hp_cs_su-self.e_hp_cs_ref)
            self.eta_hp_genven = self.eta_hp_evapen*self.eta_hp_srcen
            self.eta_hp_genvex = self.eta_hp_evapex*self.eta_hp_srcex  
            self.eta_hp_condex = self.f_hp_hs*(self.e_hp_hs_ex\
                                - self.e_hp_hs_su)/self.q_out_hp_ex
            self.eta_hp_toten  = self.pow_hp_sinken/self.pow_hp_supplen
            self.eta_hp_totex  = self.pow_hp_sinkex/self.pow_hp_supplex
        
        # --- Carnot and Lorentz efficiencies ---------------------------------
        T_cs_logAv =(self.T_hp_cs_su-self.T_hp_cs_ex)\
             /np.log(self.T_hp_cs_su/self.T_hp_cs_ex)
        T_hs_logAv =(self.T_hp_hs_ex-self.T_hp_hs_su)\
             /np.log(self.T_hp_hs_ex/self.T_hp_hs_su)
        COP_lorenz =      T_hs_logAv/(     T_hs_logAv-     T_cs_logAv)
        COP_carnot = self.T_hp_hs_ex/(self.T_hp_hs_ex-self.T_hp_cs_su)
        self.cop_hp_lorenz = COP_lorenz
        self.cop_hp_carnot = COP_carnot
        self.psi_hp_lorenz = self.eta_hp_cyclen/COP_lorenz
        self.psi_hp_carnot = self.eta_hp_cyclen/COP_carnot
        
        # --- Preliminary design information ----------------------------------
        self.vol_coef     = self.v_hp_2/self.w_in_hp_en
        self.comp_ratio   = self.p_hp_2/self.p_hp_1
        self.vol_ratio    = self.v_hp_1/self.v_hp_2
        if self.parameters['version'] == 'operational_light':
            # --- Multiply the specific quantities by the mass flow -----------
            self.Q_in_hp_en = self.m_hp*self.q_in_hp_en
            self.W_in_hp_en = self.m_hp*self.w_in_hp_en
            self.Q_out_hp_en= self.m_hp*self.q_out_hp_en
    
    def export_states(self):
        """
        Export the thermodynamic states of the SRVCHP necessary to draw
        the T-s and p-h diagrams.
        """
        self.p_hp = self.p_hp_1, self.p_hp_1x,self.p_hp_1r,self.p_hp_2, \
                    self.p_hp_2x,self.p_hp_3, self.p_hp_3x,self.p_hp_3r,\
                    self.p_hp_4, self.p_hp_4x
        self.T_hp = self.T_hp_1, self.T_hp_1x,self.T_hp_1r,self.T_hp_2, \
                    self.T_hp_2x,self.T_hp_3, self.T_hp_3x,self.T_hp_3r,\
                    self.T_hp_4, self.T_hp_4x
        self.i_hp = self.i_hp_1, self.i_hp_1x,self.i_hp_1r,self.i_hp_2, \
                    self.i_hp_2x,self.i_hp_3, self.i_hp_3x,self.i_hp_3r,\
                    self.i_hp_4, self.i_hp_4x
        self.s_hp = self.s_hp_1, self.s_hp_1x,self.s_hp_1r,self.s_hp_2, \
                    self.s_hp_2x,self.s_hp_3, self.s_hp_3x,self.s_hp_3r,\
                    self.s_hp_4, self.s_hp_4x
        self.x_hp = self.x_hp_1, self.x_hp_1x,self.x_hp_1r,self.x_hp_2, \
                    self.x_hp_2x,self.x_hp_3, self.x_hp_3x,self.x_hp_3r,\
                    self.x_hp_4, self.x_hp_4x
        if self.options['exergy']:
                    self.e_hp = self.e_hp_1, self.e_hp_1x,self.e_hp_1r,\
                                self.e_hp_2, self.e_hp_2x,self.e_hp_3, \
                                self.e_hp_3x,self.e_hp_3r,self.e_hp_4, \
                                self.e_hp_4x
        else:
            self.e_hp = 0,0,0,0,0,0,0,0,0,0
        out = self.p_hp,self.T_hp,self.i_hp,self.s_hp,self.x_hp,self.e_hp
        return out
    
    def compute_exergy(self):
        """
        Compute the exergy at each state of the SRVCHP.
        """
        self.compute_ref_states()
        
        # - 1 - Exergy of the cycle -------------------------------------------
        
        self.e_hp_1  = (self.i_hp_1 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_1 -self.s_hp_0)
        self.e_hp_1x = (self.i_hp_1x-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_1x-self.s_hp_0)
        self.e_hp_1r = (self.i_hp_1r-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_1r-self.s_hp_0)
        self.e_hp_2  = (self.i_hp_2 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_2 -self.s_hp_0)
        self.e_hp_2x = (self.i_hp_2x-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_2x-self.s_hp_0)
        self.e_hp_3  = (self.i_hp_3 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_3 -self.s_hp_0)
        self.e_hp_3x = (self.i_hp_3x-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_3x-self.s_hp_0)
        self.e_hp_3r = (self.i_hp_3r-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_3r-self.s_hp_0)
        self.e_hp_4  = (self.i_hp_4 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_4 -self.s_hp_0)
        self.e_hp_4x = (self.i_hp_4x-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_4x-self.s_hp_0)
        
        if round(self.T_hp_1 ,2)==round(self.T_hp_0,2): self.e_hp_1  = 0
        if round(self.T_hp_1x,2)==round(self.T_hp_0,2): self.e_hp_1x = 0
        if round(self.T_hp_1r,2)==round(self.T_hp_0,2): self.e_hp_1r = 0
        if round(self.T_hp_2 ,2)==round(self.T_hp_0,2): self.e_hp_2  = 0
        if round(self.T_hp_2x,2)==round(self.T_hp_0,2): self.e_hp_2x = 0
        if round(self.T_hp_3 ,2)==round(self.T_hp_0,2): self.e_hp_3  = 0
        if round(self.T_hp_3x,2)==round(self.T_hp_0,2): self.e_hp_3x = 0
        if round(self.T_hp_3r,2)==round(self.T_hp_0,2): self.e_hp_3r = 0
        if round(self.T_hp_4 ,2)==round(self.T_hp_0,2): self.e_hp_4  = 0
        if round(self.T_hp_4x,2)==round(self.T_hp_0,2): self.e_hp_4x = 0
        
        # 2 - Exergy of the boundaries ----------------------------------------
        
        self.compute_exergy_boundaries()

#%% MODEL: TBVCHP

class TBVCHP(SBVCHP):
    """
    Class for the simulation of Transcritical Basic Vapor Compression Heat 
    Pumps.

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (heat pump boundaries and operational 
                                        parameters).
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
    # ==== HEAT PUMP MODEL ================================================== #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the cold source
            * 'sink':   available heat flow rate at the hot sink
            * 'power':  electrical power supplied to the compressor
        
        """        
        if self.T_hp_hs_ex-self.parameters['dT_hp_cd_pp']\
            <self.state_hp.T_critical():
            raise ValueError('The heat source is below critical temperature!\
                             In: '+os.path.join(os.path.abspath(__file__),
                             'TBVCHP','evaluate_cycle'))
        
        # --- Cycle computation -----------------------------------------------
        if self.parameters['m_rat'] == 0: self.m_rat = self.opt_eta()
        else:                             self.m_rat = self.parameters['m_rat']

        # ..... two-phase .....................................................
        self.T_hp_4x = min(  self.T_hp_cs_su
                           - self.parameters['dT_hp_ev_pp']
                           - self.parameters['dT_hp_ev_sh'],
                             self.T_hp_cs_ex
                           - self.parameters['dT_hp_ev_pp'])
        if self.T_hp_4x > self.state_hp.T_critical():
            raise ValueError('An inconsistency was detected in the cycle!\
                              In: '+os.path.join(os.path.abspath(__file__),
                              'TBVCHP','evaluate_cycle:\
                    the evaporation temperature is above the critical point!'))
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.QT_INPUTS,0.5,self.T_hp_4x)
        self.p_hp_4x = self.state_hp.p()
        self.state_hp.unspecify_phase()
        # ..... saturated (gas) ...............................................
        self.p_hp_1x = self.p_hp_4x*(1-self.parameters['dp_hp_ev']) 
        self.x_hp_1x = 1.0
        self.state_hp.update(CoolProp.PQ_INPUTS,self.p_hp_1x,self.x_hp_1x)
        self.T_hp_1x = self.state_hp.T()
        self.i_hp_1x = self.state_hp.hmass()
        self.s_hp_1x = self.state_hp.smass()
        # ..... gas ...........................................................
        self.T_hp_1  = self.T_hp_1x+self.parameters['dT_hp_ev_sh']
        self.p_hp_1  = self.p_hp_1x
        if  self.p_hp_1 < self.state_hp.p_critical()\
        and self.T_hp_1 < self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_gas)
        if  self.p_hp_1 < self.state_hp.p_critical()\
        and self.T_hp_1 > self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_hp.update(CoolProp.PT_INPUTS,self.p_hp_1,self.T_hp_1)
        self.i_hp_1  = self.state_hp.hmass()
        self.x_hp_1  = self.state_hp.Q()
        self.s_hp_1  = self.state_hp.smass()
        self.v_hp_1  = self.state_hp.rhomass()**(-1)
        self.state_hp.unspecify_phase()
        # ..... gas (super-critical) ..........................................
        self.p_hp_2  = self.find_p()
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_hp.update(CoolProp.PSmass_INPUTS,
                                    self.p_hp_2,self.s_hp_1)
            self.i_hp_2s = self.state_hp.hmass()
            self.i_hp_2  = self.i_hp_1 \
                         +(self.i_hp_2s-self.i_hp_1)\
                         / self.parameters['eta_max_cp']
            self.m_hp    = 0
            self.P_cp    = 0
            self.eta_is_cp = self.parameters['eta_max_cp']
            self.cp_flag = 'na'
            self.state_hp.update(CoolProp.HmassP_INPUTS,
                                    self.i_hp_2,self.p_hp_2)
            self.T_hp_2  = self.state_hp.T()
        self.s_hp_2  = self.state_hp.smass()
        self.x_hp_2  = self.state_hp.Q()
        self.v_hp_2  = self.state_hp.rhomass()**(-1)
        # ..... liquid (super-critical) .......................................
        self.p_hp_3  = self.p_hp_2*(1-self.parameters['dp_hp_cd'])
        self.i_hp_3  = self.i_hp_2-self.m_rat*(self.i_hp_hs_ex-self.i_hp_hs_su)
        self.state_hp.specify_phase(CoolProp.iphase_supercritical)
        self.state_hp.update(CoolProp.HmassP_INPUTS,self.i_hp_3,self.p_hp_3)
        self.T_hp_3  = self.state_hp.T()   
        self.x_hp_3  = self.state_hp.Q()
        self.s_hp_3  = self.state_hp.smass()
        self.v_hp_3  = self.state_hp.rhomass()**(-1)
        self.state_hp.unspecify_phase()
        # ..... two-phase .....................................................
        self.p_hp_4  = self.p_hp_4x 
        self.i_hp_4  = self.i_hp_3
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.HmassP_INPUTS,self.i_hp_4,self.p_hp_4)
        self.T_hp_4  = self.state_hp.T()
        self.s_hp_4  = self.state_hp.smass()
        self.x_hp_4  = self.state_hp.Q()
        self.state_hp.unspecify_phase()
        self.i_hp_4x = self.i_hp_4
        self.T_hp_4x = self.T_hp_4
        self.s_hp_4x = self.s_hp_4
        self.x_hp_4x = self.x_hp_4
        
        # --- Evaluate the operational conditions -----------------------------
        
        if self.parameters['version'] == 'operational_light':
            self.m_hp    = self.find_mass_flow()
            self.P_cp    = self.m_hp    * (self.i_hp_2 - self.i_hp_1)
            self.m_hp_cs = self.find_flow_secondary_cs()
            self.m_hp_hs = self.find_flow_secondary_hs()
        
        # --- Evaluate the Key Performance Indicators -------------------------
        
        self.retrieve_kpi()
    
    # ======================================================================= #
    # ==== SUB-FUNCTIONS FOR HEAT PUMP MODEL ================================ #
    # ======================================================================= #
    
    def find_p(self):        
        """
        Find the condensing pressure in the TBVCHP based
        on the pinch-point and superheat, and for the selected m_rat.
        Note that phase change is not permitted for the secondary fluid.
        """
        # --- Minimum 'condenser' pressure (firm) -----------------------------
        self.p_hp_2_min     = self.state_hp.p_critical()+1
        # --- Maximum 'condenser' pressure (indicative) -----------------------
        self.p_hp_2_max     = self.state_hp.pmax()-1  
        # --- Find the pressure levels ----------------------------------------
        guess          = self.p_hp_2_min
        try:
            p_sol = scipy.optimize.fsolve(self.resi_p, guess, full_output=True)
        except:
            guess = self.p_hp_2_max
            res   = scipy.optimize.least_squares(self.resi_p, guess, 
                     bounds = (self.p_hp_2_min,self.p_hp_2_max),
                  ftol=1e-12,xtol=1e-12)
            p_sol = res.x
            import warnings
            warnings.warn('In: TSBVCHP.find_p, fsolve failed.'+
                          ' least_squares was tested instead. Message: '+
                          str(res.message))
            
        p_hp_2 = p_sol[0]
        try:    p_hp_2 = p_hp_2[0]
        except: p_hp_2 = p_hp_2
        self.resi   = self.resi_p(p_hp_2)
        return p_hp_2
         
    def resi_p(self,iter_var):
        """
        Residuals of pinch-point method in 'condenser'.  
        Used to find the condensing pressure in the TBVCHP.
        
        Cells discretization: 0.25K per cell (for secondary fluid)
        """
        # --- Guess values ----------------------------------------------------
        p_hp_2 = iter_var
        # --- Computation -----------------------------------------------------
        T_hp_4x = min(  self.T_hp_cs_su
                      - self.parameters['dT_hp_ev_pp']
                      - self.parameters['dT_hp_ev_sh'],
                        self.T_hp_cs_ex
                      - self.parameters['dT_hp_ev_pp'])
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.QT_INPUTS,0.25,T_hp_4x)
        p_hp_4x = self.state_hp.p()
        self.state_hp.unspecify_phase()
        p_hp_1x = p_hp_4x*(1-self.parameters['dp_hp_ev']) 
        x_hp_1x = 1.0
        self.state_hp.update(CoolProp.PQ_INPUTS,p_hp_1x,x_hp_1x)
        T_hp_1x = self.state_hp.T()
        self.state_hp.unspecify_phase()
        T_hp_1  = T_hp_1x+self.parameters['dT_hp_ev_sh']
        p_hp_1  = p_hp_1x
        if  p_hp_1 < self.state_hp.p_critical()\
        and T_hp_1 < self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_gas)
        if  p_hp_1 < self.state_hp.p_critical()\
        and T_hp_1 > self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_1,T_hp_1)
        h_hp_1  = self.state_hp.hmass()
        s_hp_1  = self.state_hp.smass()
        self.state_hp.unspecify_phase()
        if self.parameters['version'] == 'thermodynamic_full' \
        or self.parameters['version'] == 'operational_light':
            self.state_hp.update(CoolProp.PSmass_INPUTS,p_hp_2,s_hp_1)
            h_hp_2s  = self.state_hp.hmass()
            h_hp_2   = h_hp_1 \
                     + (h_hp_2s-h_hp_1)/self.parameters['eta_max_cp']
        
        # --- Energy balance on condenser -------------------------------------
        cell_resolution = 0.25 # K
        n_elem = int(1/cell_resolution*(self.T_hp_hs_ex-self.T_hp_hs_su)+1)
        # --- secondary fluid -------------------------------------------------
        T_hp_hs_vec       = np.linspace(self.T_hp_hs_su,self.T_hp_hs_ex,n_elem)
        p_hp_hs_vec       = np.linspace(self.p_hp_hs_su,self.p_hp_hs_ex,n_elem)
        h_hp_hs_vec  = np.zeros(n_elem)
        for i, T in enumerate(h_hp_hs_vec):
            self.state_hs.update(CoolProp.PT_INPUTS,p_hp_hs_vec[i],
                                                    T_hp_hs_vec[i])
            h_hp_hs_vec[i] = self.state_hs.hmass()
        T_hp_hs_vec       = np.flip(T_hp_hs_vec)
        h_hp_hs_vec       = np.flip(h_hp_hs_vec)
        # --- working fluid ---------------------------------------------------
        h_hp_wf_vec       = np.zeros(n_elem)
        h_hp_wf_vec[0]    = h_hp_2
        for i in range(len(h_hp_hs_vec)-1):
            h_hp_wf_vec[i+1] = h_hp_wf_vec[i]\
                              - self.m_rat*(h_hp_hs_vec[i]-h_hp_hs_vec[i+1])
        p_hp_3 = p_hp_2*(1-self.parameters['dp_hp_cd'])
        p_hp_wf_vec  = np.linspace(p_hp_2,p_hp_3,n_elem)
        T_hp_wf_vec  = np.zeros(n_elem)
        for i, T in enumerate(T_hp_wf_vec):
            self.state_hp.update(CoolProp.HmassP_INPUTS,h_hp_wf_vec[i],
                                                        p_hp_wf_vec[i])
            T_hp_wf_vec[i] = self.state_hp.T()
        # --- find pinch ------------------------------------------------------
        imin = np.argmin(T_hp_wf_vec-T_hp_hs_vec)
        T_hi = T_hp_hs_vec[imin]
        T_wi = T_hp_wf_vec[imin]
        d_hi = T_hp_wf_vec[imin]-T_hp_hs_vec[imin]\
                                -self.parameters['dT_hp_cd_pp']
        # --- Final balance ---------------------------------------------------
        out     = d_hi
        self.T_hi = T_hi
        self.T_wi = T_wi
        return out
     
    def opt_eta(self):
        """
        Find the optimum m_rat in the TBVCHP based on the pinch-point and for 
        the selected superheating.
        """
        T_hp_3_min = self.T_hp_hs_su+self.parameters['dT_hp_cd_pp']
        T_hp_3_max = self.T_hp_hs_ex+self.parameters['dT_hp_cd_pp']
        
        p_hp_3_min = self.state_hp.p_critical()+1               # firm
        p_hp_2_min = p_hp_3_min/(1-self.parameters['dp_hp_cd']) # firm
        p_hp_2_max = self.state_hp.pmax()-1                     # indicative
        p_hp_3_max = p_hp_2_max*(1-self.parameters['dp_hp_cd']) # indicative
        
        # --- Minimum 'condenser' pressure (firm) -----------------------------
        # ..... evaporator ....................................................
        T_hp_4x = min(  self.T_hp_cs_su
                      - self.parameters['dT_hp_ev_pp']
                      - self.parameters['dT_hp_ev_sh'],
                        self.T_hp_cs_ex
                      - self.parameters['dT_hp_ev_pp'])
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.QT_INPUTS,0.5,T_hp_4x)
        p_hp_4x = self.state_hp.p()
        self.state_hp.unspecify_phase()
        p_hp_1x = p_hp_4x*(1-self.parameters['dp_hp_ev']) 
        x_hp_1x = 1.0
        self.state_hp.update(CoolProp.PQ_INPUTS,p_hp_1x,x_hp_1x)
        T_hp_1x = self.state_hp.T()
        T_hp_1  = T_hp_1x+self.parameters['dT_hp_ev_sh']
        p_hp_1  = p_hp_1x
        if  p_hp_1 < self.state_hp.p_critical()\
        and T_hp_1 < self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_gas)
        if  p_hp_1 < self.state_hp.p_critical()\
        and T_hp_1 > self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_1,T_hp_1)
        h_hp_1  = self.state_hp.hmass()
        s_hp_1  = self.state_hp.smass()
        self.state_hp.unspecify_phase()
        # ..... compressor ....................................................
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            # min p2
            self.state_hp.update(CoolProp.PSmass_INPUTS,p_hp_2_min,s_hp_1)
            h_hp_2s      = self.state_hp.hmass()
            h_hp_2_minp2 = h_hp_1+(h_hp_2s-h_hp_1)\
                                                /self.parameters['eta_max_cp']
            # max p2
            try:
                self.state_hp.update(CoolProp.PSmass_INPUTS,p_hp_2_max,s_hp_1)
                h_hp_2s  = self.state_hp.hmass()
            except:
                h_hp_2s  = CP.PropsSI('H','S',s_hp_1,'P',p_hp_2_max,
                                      self.parameters['fluid_hp'])
            h_hp_2_maxp2 = h_hp_1+(h_hp_2s-h_hp_1)\
                                                /self.parameters['eta_max_cp']
        
        # ..... 'condenser' ...................................................
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_3_min,T_hp_3_min)  
        h_hp_3_minp2_minT3 = self.state_hp.hmass()
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_3_min,T_hp_3_max)  
        h_hp_3_minp2_maxT3 = self.state_hp.hmass()
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_3_max,T_hp_3_min)  
        h_hp_3_maxp2_minT3 = self.state_hp.hmass()
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_3_max,T_hp_3_max)  
        h_hp_3_maxp2_maxT3 = self.state_hp.hmass()
        den = (self.i_hp_hs_ex-self.i_hp_hs_su)
        r_hp_minp2_minT3 = (h_hp_2_minp2-h_hp_3_minp2_minT3)/den
        r_hp_minp2_maxT3 = (h_hp_2_minp2-h_hp_3_minp2_maxT3)/den
        r_hp_maxp2_minT3 = (h_hp_2_maxp2-h_hp_3_maxp2_minT3)/den
        r_hp_maxp2_maxT3 = (h_hp_2_maxp2-h_hp_3_maxp2_maxT3)/den
        # ... find the bounds and guess value for the optimisation problem ....
        r_hp_min = min([r_hp_minp2_minT3,r_hp_minp2_maxT3,
                        r_hp_maxp2_minT3,r_hp_maxp2_maxT3])
        r_hp_max = max([r_hp_minp2_minT3,r_hp_minp2_maxT3,
                        r_hp_maxp2_minT3,r_hp_maxp2_maxT3])
        r_hp_gue = r_hp_minp2_minT3#np.mean([r_hp_min,r_hp_max])
        # --- Solve the optimisation problem ----------------------------------
        try:
            sol = scipy.optimize.minimize(self.resi_eta,x0=r_hp_gue,
                                          bounds=[(max([0,   r_hp_min]),
                                                   min([1e+5,r_hp_max]))])
            res = sol.x[0]
            if res == r_hp_gue\
            or round(self.find_p()) == round(self.state_hp.pmax()-1):
                import warnings
                warnings.warn('In: TBVCHP.opt_eta, sol = guess or p = p_max. '+
                              'Full computation was done instead.')
                r_hp_gue = np.arange(r_hp_min,r_hp_max,0.025)
                out      = np.zeros(len(r_hp_gue))
                for i, r in enumerate(r_hp_gue):
                    try:    out[i] = self.resi_eta([r])
                    except: out[i] = 1e+5
                imax = np.argmax(1/out)
                # refine the results
                res = r_hp_gue[imax]
                r_hp_gue = np.arange(0.96*res,1.04*res,0.001)
                out      = np.zeros(len(r_hp_gue))
                for i, r in enumerate(r_hp_gue):
                    try:    out[i] = self.resi_eta([r])
                    except: out[i] = 1e+5
                imax = np.argmax(1/out)
                res = r_hp_gue[imax]
        except:
            import warnings
            warnings.warn('In: TBVCHP.opt_eta, minimize failed. '+
                          'Full computation was done instead.')
            r_hp_gue = np.arange(r_hp_min,r_hp_max,0.025)
            out      = np.zeros(len(r_hp_gue))
            for i, r in enumerate(r_hp_gue):
                try:    out[i] = self.resi_eta([r])
                except: out[i] = 1e+5
            imax = np.argmax(1/out)
            # refine the results
            res = r_hp_gue[imax]
            r_hp_gue = np.arange(0.96*res,1.04*res,0.001)
            out      = np.zeros(len(r_hp_gue))
            for i, r in enumerate(r_hp_gue):
                try:    out[i] = self.resi_eta([r])
                except: out[i] = 1e+5
            imax = np.argmax(1/out)
            res = r_hp_gue[imax]
        return res
    
    def resi_eta(self,iter_var):
        """
        Inverse of eta_cyclen (function to minimise).  
        Used to find the optimum m_rat in the TBVCHP.
        """
        # --- Guess values ----------------------------------------------------
        self.m_rat = iter_var[0]
        p_hp_2 = self.find_p()
        # --- Computation -----------------------------------------------------
        T_hp_4x = min(  self.T_hp_cs_su
                      - self.parameters['dT_hp_ev_pp']
                      - self.parameters['dT_hp_ev_sh'],
                        self.T_hp_cs_ex
                      - self.parameters['dT_hp_ev_pp'])
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.QT_INPUTS,0.5,T_hp_4x)
        p_hp_4x = self.state_hp.p()
        self.state_hp.unspecify_phase()
        p_hp_1x = p_hp_4x*(1-self.parameters['dp_hp_ev']) 
        x_hp_1x = 1.0
        self.state_hp.update(CoolProp.PQ_INPUTS,p_hp_1x,x_hp_1x)
        T_hp_1x = self.state_hp.T()
        self.state_hp.unspecify_phase()
        T_hp_1  = T_hp_1x+self.parameters['dT_hp_ev_sh']
        p_hp_1  = p_hp_1x
        if  p_hp_1 < self.state_hp.p_critical()\
        and T_hp_1 < self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_gas)
        if  p_hp_1 < self.state_hp.p_critical()\
        and T_hp_1 > self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_1,T_hp_1)
        h_hp_1  = self.state_hp.hmass()
        s_hp_1  = self.state_hp.smass()
        self.state_hp.unspecify_phase()
        if self.parameters['version'] == 'thermodynamic_full' \
        or self.parameters['version'] == 'operational_light':
            self.state_hp.update(CoolProp.PSmass_INPUTS,
                                 p_hp_2,s_hp_1)
            h_hp_2s  = self.state_hp.hmass()
            h_hp_2   = h_hp_1 \
                     + (h_hp_2s-h_hp_1)/self.parameters['eta_max_cp']
        
        h_hp_3  = h_hp_2-self.m_rat*(self.i_hp_hs_ex-self.i_hp_hs_su)
        h_hp_4  = h_hp_3
        # --- Proceed to efficiency -------------------------------------------
        w_in_en     = h_hp_2 -h_hp_1 
        q_out_en    = h_hp_2 -h_hp_3
        w_out_en    = h_hp_3 -h_hp_4 
        w_net_en    = w_in_en-w_out_en
        eta_cyclen  = q_out_en/w_net_en
        return 1/eta_cyclen
    
    def check_consistency(self):
        """
        Check that the results are consistent.
        """
        if self.s_hp_1x>self.s_hp_1:  self.error=True
        if self.s_hp_1 >self.s_hp_2:  self.error=True
        if self.s_hp_2 <self.s_hp_3:  self.error=True
        if self.s_hp_4 <self.s_hp_3:  self.error=True
        if self.s_hp_4x>self.s_hp_1x: self.error=True
        if self.x_hp_4 <0:            self.error=True
        if self.eta_hp_cyclen < 0:    self.error=True
        if self.T_hp_4>self.T_hp_cs_ex-0.95*self.parameters['dT_hp_ev_pp']:
            self.error=True
        if self.T_hp_1>self.T_hp_cs_su-0.95*self.parameters['dT_hp_ev_pp']:
            self.error=True
        if self.T_hp_2<self.T_hp_hs_ex+0.95*self.parameters['dT_hp_cd_pp']:
            self.error=True
        if self.T_hp_3<self.T_hp_hs_su+0.95*self.parameters['dT_hp_cd_pp']:
            self.error=True
            
        if abs(self.resi) > 1:            self.error = True
        if self.p_hp_2 < self.p_hp_2_min: self.error = True
    
    def export_states(self):
        """
        Export the thermodynamic states of the TBVCHP necessary to draw
        the T-s and p-h diagrams.
        """
        self.p_hp = self.p_hp_1,self.p_hp_1x,self.p_hp_2,\
                    self.p_hp_3,self.p_hp_4,self.p_hp_4x
        self.T_hp = self.T_hp_1,self.T_hp_1x,self.T_hp_2,\
                    self.T_hp_3,self.T_hp_4,self.T_hp_4x
        self.i_hp = self.i_hp_1,self.i_hp_1x,self.i_hp_2,\
                    self.i_hp_3,self.i_hp_4,self.i_hp_4x
        self.s_hp = self.s_hp_1,self.s_hp_1x,self.s_hp_2,\
                    self.s_hp_3,self.s_hp_4,self.s_hp_4x
        self.x_hp = self.x_hp_1,self.x_hp_1x,self.x_hp_2,\
                    self.x_hp_3,self.x_hp_4,self.x_hp_4x
        if self.options['exergy']:
            self.e_hp = self.e_hp_1,self.e_hp_1x,self.e_hp_2,\
                        self.e_hp_3,self.e_hp_4,self.e_hp_4x
        else:
            self.e_hp = 0,0,0,0,0,0,0,0
        out = self.p_hp,self.T_hp,self.i_hp,self.s_hp,self.x_hp,self.e_hp
        return out
    
    def compute_exergy(self):
        """
        Compute the exergy at each state of the TBVCHP.
        """
        self.compute_ref_states()
        
        # - 1 - Exergy of the cycle -------------------------------------------
        
        self.e_hp_1  = (self.i_hp_1 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_1 -self.s_hp_0)
        self.e_hp_1x = (self.i_hp_1x-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_1x-self.s_hp_0)
        self.e_hp_2  = (self.i_hp_2 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_2 -self.s_hp_0)
        self.e_hp_3  = (self.i_hp_3 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_3 -self.s_hp_0)
        self.e_hp_4  = (self.i_hp_4 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_4 -self.s_hp_0)
        self.e_hp_4x = (self.i_hp_4x-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_4x-self.s_hp_0)
        
        if round(self.T_hp_1 ,2)==round(self.T_hp_0,2): self.e_hp_1  = 0
        if round(self.T_hp_1x,2)==round(self.T_hp_0,2): self.e_hp_1x = 0
        if round(self.T_hp_2 ,2)==round(self.T_hp_0,2): self.e_hp_2  = 0
        if round(self.T_hp_3 ,2)==round(self.T_hp_0,2): self.e_hp_3  = 0
        if round(self.T_hp_4 ,2)==round(self.T_hp_0,2): self.e_hp_4  = 0
        if round(self.T_hp_4x,2)==round(self.T_hp_0,2): self.e_hp_4x = 0
        
        # 2 - Exergy of the boundaries ----------------------------------------
        
        self.compute_exergy_boundaries()

#%% MODEL: TRVCHP

class TRVCHP(TBVCHP):
    """
    Class for the simulation of Transcritical Recuperated Vapor Compression 
    Heat Pumps.

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (heat pump boundaries and operational 
                                        parameters).
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
    # ==== HEAT PUMP MODEL ================================================== #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the cold source
            * 'sink':   available heat flow rate at the hot sink
            * 'power':  electrical power supplied to the compressor
        
        """        
        if self.T_hp_hs_ex-self.parameters['dT_hp_cd_pp']\
            <self.state_hp.T_critical():
            raise ValueError('The heat source is below critical temperature!\
                             In: '+os.path.join(os.path.abspath(__file__),
                             'TRVCHP','evaluate_cycle'))
        # --- Cycle computation -----------------------------------------------
        if self.parameters['m_rat'] == 0: self.m_rat = self.opt_eta()
        else:                             self.m_rat = self.parameters['m_rat']

        # ..... two-phase .....................................................
        self.T_hp_4x = min(  self.T_hp_cs_su
                           - self.parameters['dT_hp_ev_pp']
                           - self.parameters['dT_hp_ev_sh'],
                             self.T_hp_cs_ex
                           - self.parameters['dT_hp_ev_pp'])
        if self.T_hp_4x > self.state_hp.T_critical():
            raise ValueError('An inconsistency was detected in the cycle!\
                              In: '+os.path.join(os.path.abspath(__file__),
                              'TRVCHP','evaluate_cycle:\
                    the evaporation temperature is above the critical point!'))
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.QT_INPUTS,0.5,self.T_hp_4x)
        self.p_hp_4x = self.state_hp.p()
        self.state_hp.unspecify_phase()
        # ..... saturated (gas) ...............................................
        self.p_hp_1x = self.p_hp_4x*(1-self.parameters['dp_hp_ev']) 
        self.x_hp_1x = 1.0
        self.state_hp.update(CoolProp.PQ_INPUTS,self.p_hp_1x,self.x_hp_1x)
        self.T_hp_1x = self.state_hp.T()
        self.i_hp_1x = self.state_hp.hmass()
        self.s_hp_1x = self.state_hp.smass()
        # ..... gas ...........................................................
        self.T_hp_1  = self.T_hp_1x+self.parameters['dT_hp_ev_sh']
        self.p_hp_1  = self.p_hp_1x
        if  self.p_hp_1 < self.state_hp.p_critical()\
        and self.T_hp_1 < self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_gas)
        if  self.p_hp_1 < self.state_hp.p_critical()\
        and self.T_hp_1 > self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_hp.update(CoolProp.PT_INPUTS,self.p_hp_1,self.T_hp_1)
        self.i_hp_1  = self.state_hp.hmass()
        self.x_hp_1  = self.state_hp.Q()
        self.s_hp_1  = self.state_hp.smass()
        self.v_hp_1  = self.state_hp.rhomass()**(-1)
        self.state_hp.unspecify_phase()
        # ..... gas (super-critical) ..........................................
        self.p_hp_2  = self.find_p()
        self.p_hp_3  = self.p_hp_2*(1-self.parameters['dp_hp_cd'])
        self.p_hp_1r = self.p_hp_1*(1-self.parameters['dp_hp_rg_vp'])
        self.i_hp_1r = self.find_recup((self.i_hp_1,self.T_hp_1,self.p_hp_1r,
                                        self.p_hp_2,self.p_hp_3))
        self.state_hp.update(CoolProp.HmassP_INPUTS,self.i_hp_1r,self.p_hp_1r)
        self.s_hp_1r = self.state_hp.smass()
        self.T_hp_1r = self.state_hp.T()
        self.x_hp_1r = self.state_hp.Q()
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_hp.update(CoolProp.PSmass_INPUTS,
                                    self.p_hp_2,self.s_hp_1r)
            self.i_hp_2s = self.state_hp.hmass()
            self.i_hp_2  = self.i_hp_1r \
                         +(self.i_hp_2s-self.i_hp_1r)\
                         / self.parameters['eta_max_cp']
            self.m_hp    = 0
            self.P_cp    = 0
            self.eta_is_cp = self.parameters['eta_max_cp']
            self.cp_flag = 'na'
            self.state_hp.update(CoolProp.HmassP_INPUTS,
                                    self.i_hp_2,self.p_hp_2)
            self.T_hp_2  = self.state_hp.T()
        
        self.s_hp_2  = self.state_hp.smass()
        self.x_hp_2  = self.state_hp.Q()
        self.v_hp_2  = self.state_hp.rhomass()**(-1)
        # ..... liquid (super-critical) .......................................
        self.i_hp_3  = self.i_hp_2-self.m_rat*(self.i_hp_hs_ex-self.i_hp_hs_su)
        self.state_hp.specify_phase(CoolProp.iphase_supercritical)
        self.state_hp.update(CoolProp.HmassP_INPUTS,self.i_hp_3,self.p_hp_3)
        self.T_hp_3  = self.state_hp.T()   
        self.x_hp_3  = self.state_hp.Q()
        self.s_hp_3  = self.state_hp.smass()
        self.v_hp_3  = self.state_hp.rhomass()**(-1)
        self.state_hp.unspecify_phase()
        # ..... liquid ........................................................
        self.i_hp_3r = self.i_hp_3-(self.i_hp_1r-self.i_hp_1)
        self.p_hp_3r = self.p_hp_3*(1-self.parameters['dp_hp_rg_lq'])
        self.state_hp.specify_phase(CoolProp.iphase_liquid)
        self.state_hp.update(CoolProp.HmassP_INPUTS,self.i_hp_3r,self.p_hp_3r)
        self.T_hp_3r = self.state_hp.T()
        self.s_hp_3r = self.state_hp.smass()
        self.x_hp_3r = self.state_hp.Q()
        self.state_hp.unspecify_phase()
        # ..... two-phase .....................................................
        self.p_hp_4  = self.p_hp_4x 
        self.i_hp_4  = self.i_hp_3r
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.HmassP_INPUTS,self.i_hp_4,self.p_hp_4)
        self.T_hp_4  = self.state_hp.T()
        self.s_hp_4  = self.state_hp.smass()
        self.x_hp_4  = self.state_hp.Q()
        self.state_hp.unspecify_phase()
        self.i_hp_4x = self.i_hp_4
        self.T_hp_4x = self.T_hp_4
        self.s_hp_4x = self.s_hp_4
        self.x_hp_4x = self.x_hp_4

        # --- Evaluate the operational conditions -----------------------------

        if self.parameters['version'] == 'operational_light':
            self.m_hp    = self.find_mass_flow()
            self.P_cp    = self.m_hp    * (self.i_hp_2 - self.i_hp_1r)
            self.m_hp_cs = self.find_flow_secondary_cs()
            self.m_hp_hs = self.find_flow_secondary_hs()

        # --- Evaluate the Key Performance Indicators -------------------------

        self.retrieve_kpi()
    
    # ======================================================================= #
    # ==== SUB-FUNCTIONS FOR HEAT PUMP MODEL ================================ #
    # ======================================================================= #
         
    def resi_p(self,iter_var):
        """
        Residuals of pinch-point method in 'condenser'.  
        Used to find the condensing pressure in the TRVCHP.
        
        Cells discretization: 0.25K per cell (for secondary fluid)
        """
        # --- Guess values ----------------------------------------------------
        p_hp_2 = iter_var
        # --- Computation -----------------------------------------------------
        T_hp_4x = min(  self.T_hp_cs_su
                      - self.parameters['dT_hp_ev_pp']
                      - self.parameters['dT_hp_ev_sh'],
                        self.T_hp_cs_ex
                      - self.parameters['dT_hp_ev_pp'])
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.QT_INPUTS,0.25,T_hp_4x)
        p_hp_4x = self.state_hp.p()
        self.state_hp.unspecify_phase()
        p_hp_1x = p_hp_4x*(1-self.parameters['dp_hp_ev']) 
        x_hp_1x = 1.0
        self.state_hp.update(CoolProp.PQ_INPUTS,p_hp_1x,x_hp_1x)
        T_hp_1x = self.state_hp.T()
        self.state_hp.unspecify_phase()
        T_hp_1  = T_hp_1x+self.parameters['dT_hp_ev_sh']
        p_hp_1  = p_hp_1x
        if  p_hp_1 < self.state_hp.p_critical()\
        and T_hp_1 < self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_gas)
        if  p_hp_1 < self.state_hp.p_critical()\
        and T_hp_1 > self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_1,T_hp_1)
        h_hp_1  = self.state_hp.hmass()
        self.state_hp.unspecify_phase()
        p_hp_1r = p_hp_1*(1-self.parameters['dp_hp_rg_vp'])
        p_hp_3  = p_hp_2*(1-self.parameters['dp_hp_cd'])
        extra   = h_hp_1,T_hp_1,p_hp_1r,p_hp_2,p_hp_3
        h_hp_1r = self.find_recup(extra)
        self.state_hp.update(CoolProp.HmassP_INPUTS,h_hp_1r,p_hp_1r)
        s_hp_1r = self.state_hp.smass()
        T_hp_1r = self.state_hp.T()
        if self.parameters['version'] == 'thermodynamic_full' \
        or self.parameters['version'] == 'operational_light':
            self.state_hp.update(CoolProp.PSmass_INPUTS,p_hp_2,s_hp_1r)
            h_hp_2s  = self.state_hp.hmass()
            h_hp_2   = h_hp_1r \
                      + (h_hp_2s-h_hp_1r)/self.parameters['eta_max_cp']
        
        # --- Energy balance on condenser -------------------------------------
        cell_resolution   = 0.25 # K
        n_elem = int(1/cell_resolution*(self.T_hp_hs_ex-self.T_hp_hs_su)+1)
        # --- secondary fluid -------------------------------------------------
        T_hp_hs_vec       = np.linspace(self.T_hp_hs_su,self.T_hp_hs_ex,n_elem)
        p_hp_hs_vec       = np.linspace(self.p_hp_hs_su,self.p_hp_hs_ex,n_elem)
        h_hp_hs_vec  = np.zeros(n_elem)
        for i, T in enumerate(h_hp_hs_vec):
            self.state_hs.update(CoolProp.PT_INPUTS,p_hp_hs_vec[i],
                                                    T_hp_hs_vec[i])
            h_hp_hs_vec[i] = self.state_hs.hmass()
        T_hp_hs_vec       = np.flip(T_hp_hs_vec)
        h_hp_hs_vec       = np.flip(h_hp_hs_vec)
        # --- working fluid ---------------------------------------------------
        h_hp_wf_vec       = np.zeros(n_elem)
        h_hp_wf_vec[0]    = h_hp_2
        for i in range(len(h_hp_hs_vec)-1):
            h_hp_wf_vec[i+1] = h_hp_wf_vec[i]\
                              - self.m_rat*(h_hp_hs_vec[i]-h_hp_hs_vec[i+1])
        p_hp_3 = p_hp_2*(1-self.parameters['dp_hp_cd'])
        p_hp_wf_vec  = np.linspace(p_hp_2,p_hp_3,n_elem)
        T_hp_wf_vec  = np.zeros(n_elem)
        for i, T in enumerate(T_hp_wf_vec):
            self.state_hp.update(CoolProp.HmassP_INPUTS,h_hp_wf_vec[i],
                                                        p_hp_wf_vec[i])
            T_hp_wf_vec[i] = self.state_hp.T()
        # --- find pinch ------------------------------------------------------
        imin    = np.argmin(T_hp_wf_vec-T_hp_hs_vec)
        T_hi    = T_hp_hs_vec[imin]
        T_wi    = T_hp_wf_vec[imin]
        d_hi    = T_hp_wf_vec[imin]-T_hp_hs_vec[imin]\
                                   -self.parameters['dT_hp_cd_pp']
        # --- Final balance ---------------------------------------------------
        out     = d_hi
        self.T_hi = T_hi
        self.T_wi = T_wi
        return out
     
    def find_recup(self,extra):
        """
        Find the condensing pressure in the TRVCHP based
        on the pinch-point, and for the selected m_rat.
        """
        # --- Find the pressure levels ----------------------------------------
        h_hp_1,T_hp_1,p_hp_1r,p_hp_2,p_hp_3 = extra

        # --- Find the temperatures -------------------------------------------
        guess = h_hp_1*1.1
        
        try:
            res   = scipy.optimize.fsolve       (self.resi_recup, guess, extra, 
                                          full_output=True)
            h_sol = res[0]
        except:
            res   = scipy.optimize.least_squares(self.resi_recup, guess, extra, 
                                               bounds = (h_hp_1,h_hp_1*1e+3))
            h_sol = res.x
            import warnings
            warnings.warn('In: TRVCHP.find_recup, fsolve failed.'+
                          ' least_squares was tested instead. Message: '+
                          str(res.message))
        try:    h_hp_1r = h_sol[0]
        except: h_hp_1r = h_sol
        self.resi_rec = self.resi_recup(h_hp_1r,*extra)
        return h_hp_1r
    
    def resi_recup(self,iter_var,*extra):
        """
        Residuals of effectiveness method in recuperator.
        Used to find the exit temperatures in the TRVCHP.
        """
        # --- Guess values ----------------------------------------------------
        try:    h_hp_1r_guess = iter_var[0]
        except: h_hp_1r_guess = iter_var
        h_hp_1,T_hp_1,p_hp_1r,p_hp_2,p_hp_3 = extra
        # --- Computation -----------------------------------------------------
        self.state_hp.update(CoolProp.HmassP_INPUTS,h_hp_1r_guess,p_hp_1r)
        T_hp_1r_guess = self.state_hp.T() 
        s_hp_1r_guess = self.state_hp.smass()
        if self.parameters['version'] == 'thermodynamic_full' \
        or self.parameters['version'] == 'operational_light':
            self.state_hp.update(CoolProp.PSmass_INPUTS,p_hp_2,s_hp_1r_guess)
            h_hp_2s  = self.state_hp.hmass()
            h_hp_2   =   h_hp_1r_guess \
                      + (h_hp_2s-h_hp_1r_guess)/self.parameters['eta_max_cp']
        
        h_hp_3 = h_hp_2-self.m_rat*(self.i_hp_hs_ex-self.i_hp_hs_su)
        self.state_hp.update(CoolProp.HmassP_INPUTS,h_hp_3,p_hp_3)
        T_hp_3  = self.state_hp.T()   
        T_hp_1r = T_hp_1+(T_hp_3-T_hp_1)*self.parameters['epsilon']
        # --- Final energy balance --------------------------------------------
        out  = T_hp_1r-T_hp_1r_guess
        return out
        
    def opt_eta(self):
        """
        Find the optimum m_rat in the TRVCHP based on the pinch-point and for 
        the selected superheating.
        """
        T_hp_3_min = self.T_hp_hs_su+self.parameters['dT_hp_cd_pp']
        T_hp_3_max = self.T_hp_hs_ex+self.parameters['dT_hp_cd_pp']
        
        p_hp_3_min = self.state_hp.p_critical()+1               # firm
        p_hp_2_min = p_hp_3_min/(1-self.parameters['dp_hp_cd']) # firm
        p_hp_2_max = self.state_hp.pmax()-1                     # indicative
        p_hp_3_max = p_hp_2_max*(1-self.parameters['dp_hp_cd']) # indicative
        
        # --- Minimum 'condenser' pressure (firm) -----------------------------
        # ..... evaporator ....................................................
        T_hp_4x = min(  self.T_hp_cs_su
                      - self.parameters['dT_hp_ev_pp']
                      - self.parameters['dT_hp_ev_sh'],
                        self.T_hp_cs_ex
                      - self.parameters['dT_hp_ev_pp'])
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.QT_INPUTS,0.5,T_hp_4x)
        p_hp_4x = self.state_hp.p()
        self.state_hp.unspecify_phase()
        p_hp_1x = p_hp_4x*(1-self.parameters['dp_hp_ev']) 
        x_hp_1x = 1.0
        self.state_hp.update(CoolProp.PQ_INPUTS,p_hp_1x,x_hp_1x)
        T_hp_1x = self.state_hp.T()
        T_hp_1  = T_hp_1x+self.parameters['dT_hp_ev_sh']
        p_hp_1  = p_hp_1x
        if  p_hp_1 < self.state_hp.p_critical()\
        and T_hp_1 < self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_gas)
        if  p_hp_1 < self.state_hp.p_critical()\
        and T_hp_1 > self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_1,T_hp_1)
        h_hp_1  = self.state_hp.hmass()
        s_hp_1  = self.state_hp.smass()
        self.state_hp.unspecify_phase()
        
        p_hp_1r     = p_hp_1*(1-self.parameters['dp_hp_rg_vp'])
        h_hp_1r_min = h_hp_1
        s_hp_1r_min = s_hp_1
        T_hp_1r_min = T_hp_1
        h_hp_1r_max = h_hp_1
        s_hp_1r_max = s_hp_1
        T_hp_1r_max = T_hp_1
        
        # ..... compressor ....................................................
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            # min p2
            self.state_hp.update(CoolProp.PSmass_INPUTS,p_hp_2_min,s_hp_1r_min)
            h_hp_2s      = self.state_hp.hmass()
            h_hp_2_minp2 = h_hp_1r_min+(h_hp_2s-h_hp_1r_min)\
                                                /self.parameters['eta_max_cp']
            # max p2
            try:
                self.state_hp.update(CoolProp.PSmass_INPUTS,p_hp_2_max,
                                     s_hp_1r_max)
                h_hp_2s  = self.state_hp.hmass()
            except:
                h_hp_2s  = CP.PropsSI('H','S',s_hp_1r_max,'P',p_hp_2_max,
                                      self.parameters['fluid_hp'])
            h_hp_2_maxp2 = h_hp_1r_max+(h_hp_2s-h_hp_1r_max)\
                                                /self.parameters['eta_max_cp']
        
        # ..... 'condenser' ...................................................
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_3_min,T_hp_3_min)  
        h_hp_3_minp2_minT3 = self.state_hp.hmass()
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_3_min,T_hp_3_max)  
        h_hp_3_minp2_maxT3 = self.state_hp.hmass()
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_3_max,T_hp_3_min)  
        h_hp_3_maxp2_minT3 = self.state_hp.hmass()
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_3_max,T_hp_3_max)  
        h_hp_3_maxp2_maxT3 = self.state_hp.hmass()
        den = (self.i_hp_hs_ex-self.i_hp_hs_su)
        r_hp_minp2_minT3 = (h_hp_2_minp2-h_hp_3_minp2_minT3)/den
        r_hp_minp2_maxT3 = (h_hp_2_minp2-h_hp_3_minp2_maxT3)/den
        r_hp_maxp2_minT3 = (h_hp_2_maxp2-h_hp_3_maxp2_minT3)/den
        r_hp_maxp2_maxT3 = (h_hp_2_maxp2-h_hp_3_maxp2_maxT3)/den
        # ... find the bounds and guess value for the optimisation problem ....
        r_hp_min = min([r_hp_minp2_minT3,r_hp_minp2_maxT3,
                        r_hp_maxp2_minT3,r_hp_maxp2_maxT3])
        r_hp_max = max([r_hp_minp2_minT3,r_hp_minp2_maxT3,
                        r_hp_maxp2_minT3,r_hp_maxp2_maxT3])
        r_hp_gue = r_hp_minp2_minT3#np.mean([r_hp_min,r_hp_max])
        # --- Solve the optimisation problem ----------------------------------
        try:
            sol = scipy.optimize.minimize(self.resi_eta,x0=r_hp_gue,
                                          bounds=[(max([0,   r_hp_min]),
                                                    min([1e+5,r_hp_max]))])
            res = sol.x[0]
            if res == r_hp_gue\
            or round(self.find_p()) == round(self.state_hp.pmax()-1):
                import warnings
                warnings.warn('In: TRVCHP.opt_eta, sol = guess or p = p_max. '+
                              'Full computation was done instead.')
                r_hp_gue = np.arange(r_hp_min,r_hp_max,0.025)
                out      = np.zeros(len(r_hp_gue))
                for i, r in enumerate(r_hp_gue):
                    try:    out[i] = self.resi_eta([r])
                    except: out[i] = 1e+5
                imax = np.argmax(1/out)
                # refine the results
                res = r_hp_gue[imax]
                r_hp_gue = np.arange(0.96*res,1.04*res,0.001)
                out      = np.zeros(len(r_hp_gue))
                for i, r in enumerate(r_hp_gue):
                    try:    out[i] = self.resi_eta([r])
                    except: out[i] = 1e+5
                imax = np.argmax(1/out)
                res = r_hp_gue[imax]
        except:
            import warnings
            warnings.warn('In: TRVCHP.opt_eta, minimize failed. '+
                          'Full computation was done instead.')
            r_hp_gue = np.arange(r_hp_min,r_hp_max,0.025)
            out      = np.zeros(len(r_hp_gue))
            for i, r in enumerate(r_hp_gue):
                try:    out[i] = self.resi_eta([r])
                except: out[i] = 1e+5
            imax = np.argmax(1/out)
            # refine the results
            res = r_hp_gue[imax]
            r_hp_gue = np.arange(0.96*res,1.04*res,0.001)
            out      = np.zeros(len(r_hp_gue))
            for i, r in enumerate(r_hp_gue):
                try:    out[i] = self.resi_eta([r])
                except: out[i] = 1e+5
            imax = np.argmax(1/out)
            res = r_hp_gue[imax]
        return res
    
    def resi_eta(self,iter_var):
        """
        Inverse of eta_cyclen (function to minimise).  
        Used to find the optimum m_rat in the TRVCHP.
        """
        # --- Guess values ----------------------------------------------------
        self.m_rat = iter_var[0]
        p_hp_2 = self.find_p()
        # --- Computation -----------------------------------------------------
        T_hp_4x = min(  self.T_hp_cs_su
                      - self.parameters['dT_hp_ev_pp']
                      - self.parameters['dT_hp_ev_sh'],
                        self.T_hp_cs_ex
                      - self.parameters['dT_hp_ev_pp'])
        self.state_hp.specify_phase(CoolProp.iphase_twophase)
        self.state_hp.update(CoolProp.QT_INPUTS,0.25,T_hp_4x)
        p_hp_4x = self.state_hp.p()
        self.state_hp.unspecify_phase()
        p_hp_1x = p_hp_4x*(1-self.parameters['dp_hp_ev']) 
        x_hp_1x = 1.0
        self.state_hp.update(CoolProp.PQ_INPUTS,p_hp_1x,x_hp_1x)
        T_hp_1x = self.state_hp.T()
        self.state_hp.unspecify_phase()
        T_hp_1  = T_hp_1x+self.parameters['dT_hp_ev_sh']
        p_hp_1  = p_hp_1x
        if  p_hp_1 < self.state_hp.p_critical()\
        and T_hp_1 < self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_gas)
        if  p_hp_1 < self.state_hp.p_critical()\
        and T_hp_1 > self.state_hp.T_critical():
            self.state_hp.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_hp.update(CoolProp.PT_INPUTS,p_hp_1,T_hp_1)
        h_hp_1  = self.state_hp.hmass()
        self.state_hp.unspecify_phase()
        p_hp_1r = p_hp_1*(1-self.parameters['dp_hp_rg_vp'])
        p_hp_3  = p_hp_2*(1-self.parameters['dp_hp_cd'])
        extra   = h_hp_1,T_hp_1,p_hp_1r,p_hp_2,p_hp_3
        h_hp_1r = self.find_recup(extra)
        self.state_hp.update(CoolProp.HmassP_INPUTS,h_hp_1r,p_hp_1r)
        s_hp_1r = self.state_hp.smass()
        T_hp_1r = self.state_hp.T()
        
        if self.parameters['version'] == 'thermodynamic_full' \
        or self.parameters['version'] == 'operational_light':
            self.state_hp.update(CoolProp.PSmass_INPUTS,
                                  p_hp_2,s_hp_1r)
            h_hp_2s  = self.state_hp.hmass()
            h_hp_2   = h_hp_1r \
                      + (h_hp_2s-h_hp_1r)/self.parameters['eta_max_cp']
        
        h_hp_3  = h_hp_2-self.m_rat*(self.i_hp_hs_ex-self.i_hp_hs_su)
        h_hp_3r = h_hp_3-(h_hp_1r-h_hp_1)
        h_hp_4  = h_hp_3r
        # --- Proceed to efficiency -------------------------------------------
        w_in_en     = h_hp_2 -h_hp_1r 
        q_out_en    = h_hp_2 -h_hp_3
        w_out_en    = h_hp_3r-h_hp_4 
        w_net_en    = w_in_en-w_out_en
        eta_cyclen  = q_out_en/w_net_en
        return 1/eta_cyclen
    
    def check_consistency(self):
        """
        Check that the results are consistent.
        """
        
        if self.s_hp_1x>self.s_hp_1:  self.error=True
        if self.s_hp_1 >self.s_hp_1r: self.error=True        
        if self.s_hp_2 <self.s_hp_3:  self.error=True
        if self.s_hp_3 <self.s_hp_3r: self.error=True
        if self.s_hp_4 <self.s_hp_3r: self.error=True
        if self.s_hp_4x>self.s_hp_1x: self.error=True
        if self.x_hp_4 <0:            self.error=True
        if self.eta_hp_cyclen < 0:    self.error=True
        
        if self.T_hp_4>self.T_hp_cs_ex-0.95*self.parameters['dT_hp_ev_pp']:
            self.error=True
        if self.T_hp_1>self.T_hp_cs_su-0.95*self.parameters['dT_hp_ev_pp']:
            self.error=True
        if self.T_hp_2<self.T_hp_hs_ex+0.95*self.parameters['dT_hp_cd_pp']:
            self.error=True
        if self.T_hp_3<self.T_hp_hs_su+0.95*self.parameters['dT_hp_cd_pp']:
            self.error=True
        if self.T_hp_1r>self.T_hp_3   -0.95*self.parameters['dT_hp_cd_pp']:#!!!
            self.error = True
        
        if abs(self.resi) > 1:              self.error = True
        if self.p_hp_2 < self.p_hp_2_min:   self.error = True
    
    def retrieve_kpi(self):
        """
        Retrieve the KPI's of the heat pump.
        """
        
        if self.options['exergy']: self.compute_exergy()
        
        # --- Basic quantities ------------------------------------------------
        self.w_in_hp_en = self.i_hp_2 -self.i_hp_1r
        self.q_in_hp_en = self.i_hp_1 -self.i_hp_4
        self.w_out_hp_en= self.i_hp_3r-self.i_hp_4
        self.q_out_hp_en= self.i_hp_2 -self.i_hp_3
        self.f_hp_hs    = self.q_out_hp_en/(self.i_hp_hs_ex-self.i_hp_hs_su)
        self.f_hp_cs    = self.q_in_hp_en /(self.i_hp_cs_su-self.i_hp_cs_ex)
        
        if self.options['exergy']: 
            self.w_in_hp_ex = self.e_hp_2 -self.e_hp_1r
            self.q_in_hp_ex = self.e_hp_1 -self.e_hp_4
            self.w_out_hp_ex= self.e_hp_3r-self.e_hp_4 
            self.q_out_hp_ex= self.e_hp_2 -self.e_hp_3 
        
        # --- Losses ----------------------------------------------------------
        if self.options['exergy']:
            self.loss_hp_compex = self.w_in_hp_en-self.w_in_hp_ex
            self.loss_hp_valvex = self.w_out_hp_ex-self.w_out_hp_en
            self.loss_hp_regex  =(self.e_hp_3 -self.e_hp_3r)\
                                -(self.e_hp_1r-self.e_hp_1 )
            self.loss_hp_evapex = self.f_hp_cs*(self.e_hp_cs_su\
                                - self.e_hp_cs_ex)-self.q_in_hp_ex
            self.loss_hp_srcen  = self.f_hp_cs*(self.i_hp_cs_ex\
                                - self.i_hp_cs_ref)
            self.loss_hp_srcex  = self.f_hp_cs*(self.e_hp_cs_ex\
                                - self.e_hp_cs_ref)
            self.loss_hp_condex = self.q_out_hp_ex-self.f_hp_hs\
                                *(self.e_hp_hs_ex-self.e_hp_hs_su)
            self.pow_hp_sinken  = self.f_hp_hs*(self.i_hp_hs_ex\
                                - self.i_hp_hs_su)
            self.pow_hp_sinkex  = self.f_hp_hs*(self.e_hp_hs_ex\
                                - self.e_hp_hs_su)
            self.pow_hp_supplen = self.f_hp_cs*(self.i_hp_cs_su\
                                - self.i_hp_cs_ref)+self.w_in_hp_en
            self.pow_hp_supplex = self.f_hp_cs*(self.e_hp_cs_su\
                                - self.e_hp_cs_ref)+self.w_in_hp_en
        
        # if not exergy mode, still need supplex for exergy efficiency
        if not self.options['exergy']:
            self.compute_ref_states()
            self.compute_exergy_boundaries()
                
            self.pow_hp_supplex = self.f_hp_cs*(self.e_hp_cs_su\
                                - self.e_hp_cs_ref)+self.w_in_hp_en
        
        # --- Efficiencies ----------------------------------------------------
        self.eta_hp_cyclen = self.q_out_hp_en/self.w_in_hp_en
        
        if self.options['exergy']:
            self.eta_hp_cyclex = self.q_out_hp_ex/self.w_in_hp_en
            self.eta_hp_compex = self.w_in_hp_ex /self.w_in_hp_en
            self.eta_hp_valvex = self.w_out_hp_en/self.w_out_hp_ex
            self.eta_hp_regex  =(self.e_hp_1r-self.e_hp_1 )\
                               /(self.e_hp_3 -self.e_hp_3r)
            self.eta_hp_evapen = 1
            self.eta_hp_evapex = self.f_hp_cs**(-1)*self.q_in_hp_ex\
                               /(self.e_hp_cs_su-self.e_hp_cs_ex)
            if self.i_hp_cs_su == self.i_hp_cs_ref: 
                self.eta_hp_srcen  = float('inf')
            else:                                   
                self.eta_hp_srcen  = (self.i_hp_cs_su-self.i_hp_cs_ex)\
                                    /(self.i_hp_cs_su-self.i_hp_cs_ref)
            if self.e_hp_cs_su == self.e_hp_cs_ref: 
                self.eta_hp_srcex  = float('inf')
            else:                                   
                self.eta_hp_srcex  = (self.e_hp_cs_su-self.e_hp_cs_ex)\
                                    /(self.e_hp_cs_su-self.e_hp_cs_ref)
            self.eta_hp_genven = self.eta_hp_evapen*self.eta_hp_srcen
            self.eta_hp_genvex = self.eta_hp_evapex*self.eta_hp_srcex  
            self.eta_hp_condex = self.f_hp_hs*(self.e_hp_hs_ex\
                                - self.e_hp_hs_su)/self.q_out_hp_ex
            self.eta_hp_toten  = self.pow_hp_sinken/self.pow_hp_supplen
            self.eta_hp_totex  = self.pow_hp_sinkex/self.pow_hp_supplex
        
        # --- Carnot and Lorentz efficiencies ---------------------------------
        T_cs_logAv =(self.T_hp_cs_su-self.T_hp_cs_ex)\
             /np.log(self.T_hp_cs_su/self.T_hp_cs_ex)
        T_hs_logAv =(self.T_hp_hs_ex-self.T_hp_hs_su)\
             /np.log(self.T_hp_hs_ex/self.T_hp_hs_su)
        COP_lorenz =      T_hs_logAv/(     T_hs_logAv-     T_cs_logAv)
        COP_carnot = self.T_hp_hs_ex/(self.T_hp_hs_ex-self.T_hp_cs_su)
        self.cop_hp_lorenz = COP_lorenz
        self.cop_hp_carnot = COP_carnot
        self.psi_hp_lorenz = self.eta_hp_cyclen/COP_lorenz
        self.psi_hp_carnot = self.eta_hp_cyclen/COP_carnot
        
        # --- Preliminary design information ----------------------------------
        self.vol_coef     = self.v_hp_2/self.w_in_hp_en
        self.comp_ratio   = self.p_hp_2/self.p_hp_1
        self.vol_ratio    = self.v_hp_1/self.v_hp_2
        if self.parameters['version'] == 'operational_light':
            # --- Multiply the specific quantities by the mass flow -----------
            self.Q_in_hp_en = self.m_hp*self.q_in_hp_en
            self.W_in_hp_en = self.m_hp*self.w_in_hp_en
            self.Q_out_hp_en= self.m_hp*self.q_out_hp_en
        
    def export_states(self):
        """
        Export the thermodynamic states of the SRVCHP necessary to draw
        the T-s and p-h diagrams.
        """
        self.p_hp = self.p_hp_1, self.p_hp_1x,self.p_hp_1r,self.p_hp_2, \
                    self.p_hp_3, self.p_hp_3r,self.p_hp_4, self.p_hp_4x
        self.T_hp = self.T_hp_1, self.T_hp_1x,self.T_hp_1r,self.T_hp_2, \
                    self.T_hp_3, self.T_hp_3r,self.T_hp_4, self.T_hp_4x
        self.i_hp = self.i_hp_1, self.i_hp_1x,self.i_hp_1r,self.i_hp_2, \
                    self.i_hp_3, self.i_hp_3r,self.i_hp_4, self.i_hp_4x
        self.s_hp = self.s_hp_1, self.s_hp_1x,self.s_hp_1r,self.s_hp_2, \
                    self.s_hp_3, self.s_hp_3r,self.s_hp_4, self.s_hp_4x
        self.x_hp = self.x_hp_1, self.x_hp_1x,self.x_hp_1r,self.x_hp_2, \
                    self.x_hp_3, self.x_hp_3r,self.x_hp_4, self.x_hp_4x
        if self.options['exergy']:
                    self.e_hp = self.e_hp_1, self.e_hp_1x,self.e_hp_1r,\
                                self.e_hp_2, self.e_hp_3, self.e_hp_3r,\
                                self.e_hp_4, self.e_hp_4x
        else:
            self.e_hp = 0,0,0,0,0,0,0,0
        out = self.p_hp,self.T_hp,self.i_hp,self.s_hp,self.x_hp,self.e_hp
        return out
    
    def compute_exergy(self):
        """
        Compute the exergy at each state of the SRVCHP.
        """
        self.compute_ref_states()
        
        # - 1 - Exergy of the cycle -------------------------------------------
        
        self.e_hp_1  = (self.i_hp_1 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_1 -self.s_hp_0)
        self.e_hp_1x = (self.i_hp_1x-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_1x-self.s_hp_0)
        self.e_hp_1r = (self.i_hp_1r-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_1r-self.s_hp_0)
        self.e_hp_2  = (self.i_hp_2 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_2 -self.s_hp_0)
        self.e_hp_3  = (self.i_hp_3 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_3 -self.s_hp_0)
        self.e_hp_3r = (self.i_hp_3r-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_3r-self.s_hp_0)
        self.e_hp_4  = (self.i_hp_4 -self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_4 -self.s_hp_0)
        self.e_hp_4x = (self.i_hp_4x-self.i_hp_0)\
          -self.T_hp_0*(self.s_hp_4x-self.s_hp_0)
        
        if round(self.T_hp_1 ,2)==round(self.T_hp_0,2): self.e_hp_1  = 0
        if round(self.T_hp_1x,2)==round(self.T_hp_0,2): self.e_hp_1x = 0
        if round(self.T_hp_1r,2)==round(self.T_hp_0,2): self.e_hp_1r = 0
        if round(self.T_hp_2 ,2)==round(self.T_hp_0,2): self.e_hp_2  = 0
        if round(self.T_hp_3 ,2)==round(self.T_hp_0,2): self.e_hp_3  = 0
        if round(self.T_hp_3r,2)==round(self.T_hp_0,2): self.e_hp_3r = 0
        if round(self.T_hp_4 ,2)==round(self.T_hp_0,2): self.e_hp_4  = 0
        if round(self.T_hp_4x,2)==round(self.T_hp_0,2): self.e_hp_4x = 0
        
        # 2 - Exergy of the boundaries ----------------------------------------
        
        self.compute_exergy_boundaries()