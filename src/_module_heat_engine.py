"""
- Date: March 2025
- Author: Antoine Laterre

@ Project: "CBSim - Carnot Battery Simulator" 

-> Module containing the heat engines models

"""

#%% IMPORT PACKAGES

import scipy.optimize
import numpy as np

import CoolProp

from   CoolProp import AbstractState

import os
import sys
SIMULATOR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(SIMULATOR,'COMPONENTS'))

#%% MODEL: SBORC

class SBORC:
    """
    Class for the simulation of Subcritical Basic Organic Rankine Cycles.

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (heat engine boundaries and operational 
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
        Create a SBORC object.
        """
        
        p_he_cs_su, i_he_cs_su, p_he_cs_ex, i_he_cs_ex, m_he_cs, fluid_he_cs, \
        p_he_hs_su, i_he_hs_su, p_he_hs_ex, i_he_hs_ex, m_he_hs, fluid_he_hs, \
        P_he = inputs        
        
        # --- Heat engine -----------------------------------------------------
        self.P_he        = P_he       # [W]   heat engine net power output
        backend          = 'TTSE'     # HEOS for Helm. Eq. Of State (slower)
        self.state_he = AbstractState(backend,parameters['fluid_he'])
        self.state_he.unspecify_phase()
        
        # --- Heat source -----------------------------------------------------
        self.p_he_hs_su  = p_he_hs_su # [Pa]   hot  heat source supply pressure
        self.i_he_hs_su  = i_he_hs_su # [J/kg] hot  heat source supply enthalpy
        self.p_he_hs_ex  = p_he_hs_ex # [Pa]   hot  heat source exit pressure
        self.i_he_hs_ex  = i_he_hs_ex # [J/kg] hot  heat source exit   enthalpy
        self.m_he_hs     = m_he_hs    # [kg/s] hot  heat source mass flow rate
        self.fluid_he_hs = fluid_he_hs# [-]    hot  heat source fluid
        self.state_hs    = AbstractState(backend,self.fluid_he_hs)
        self.state_hs.unspecify_phase()
        
        self.state_hs.update(CoolProp.HmassP_INPUTS,self.i_he_hs_su,
                                                    self.p_he_hs_su)
        self.T_he_hs_su  = self.state_hs.T()
        self.s_he_hs_su  = self.state_hs.smass()
        self.state_hs.update(CoolProp.HmassP_INPUTS,self.i_he_hs_ex,
                                                    self.p_he_hs_ex)
        self.T_he_hs_ex  = self.state_hs.T()
        self.s_he_hs_ex  = self.state_hs.smass()

        if self.T_he_hs_su < self.T_he_hs_ex:
            raise ValueError('The glide on the evaporator must be positive!\
                             In: '+os.path.join(os.path.abspath(__file__),
                             'SBORC','__init__'))
    
        # --- Heat sink -------------------------------------------------------
        self.p_he_cs_su  = p_he_cs_su # [Pa]   cold heat sink   supply pressure
        self.i_he_cs_su  = i_he_cs_su # [J/kg] cold heat sink   supply enthalpy
        self.p_he_cs_ex  = p_he_cs_ex # [Pa]   cold heat sink   exit   pressure
        self.i_he_cs_ex  = i_he_cs_ex # [J/kg] cold heat sink   exit   enthalpy
        self.m_he_cs     = m_he_cs    # [kg/s] cold heat sink   mass flow rate
        self.fluid_he_cs = fluid_he_cs# [-]    cold heat sink   fluid
        self.state_cs    = AbstractState(backend,self.fluid_he_cs)
        self.state_cs.unspecify_phase()
        
        self.state_cs.update(CoolProp.HmassP_INPUTS,self.i_he_cs_su,
                                                    self.p_he_cs_su)
        self.T_he_cs_su  = self.state_cs.T()
        self.s_he_cs_su  = self.state_cs.smass()
        self.state_cs.update(CoolProp.HmassP_INPUTS,self.i_he_cs_ex,
                                                    self.p_he_cs_ex)
        self.T_he_cs_ex  = self.state_cs.T()
        self.s_he_cs_ex  = self.state_cs.smass()
        
        if self.T_he_cs_ex < self.T_he_cs_su:
            raise ValueError('The glide on the evaporator must be positive!\
                             In: '+os.path.join(os.path.abspath(__file__),
                             'SBORC','__init__'))
        
        # --- Parameters and options ------------------------------------------
        self.parameters  = parameters
        self.options     = options
    
    def evaluate(self):
        """
        Evaluate the SBORC cycle for the given boundary conditions.
        """
        
        self.error = True
        
        if self.parameters['version'] not in ['thermodynamic_full',
                                              'operational_light']:
            raise ValueError('An inconsistency was detected!\
                              In: '+os.path.join(os.path.abspath(__file__),
                'SBORC','evaluate: parameters["version"] is not valid'))
        
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
                              'SBORC','check_consistency'))
    
    # ======================================================================= #
    # ==== HEAT ENGINE MODEL ================================================ #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the hot source
            * 'sink':   available heat flow rate at the cold sink
            * 'power':  electrical power produced by the heat engine
        
        """        
        
        # --- Cycle computation -----------------------------------------------
        
        # ..... saturated (liquid) ............................................
        self.p_he_1x, self.p_he_3x = self.find_p()
        self.x_he_1x = 0.0
        self.state_he.update(CoolProp.PQ_INPUTS,self.p_he_1x,self.x_he_1x)
        self.T_he_1x = self.state_he.T()
        self.i_he_1x = self.state_he.hmass()
        self.s_he_1x = self.state_he.smass()
        # ..... liquid ........................................................
        self.p_he_1  = self.p_he_1x
        self.T_he_1  = self.T_he_1x-self.parameters['dT_he_cd_sc']
        self.state_he.specify_phase(CoolProp.iphase_liquid)
        self.state_he.update(CoolProp.PT_INPUTS,self.p_he_1,self.T_he_1)
        self.i_he_1  = self.state_he.hmass()
        self.x_he_1  = self.state_he.Q()
        self.s_he_1  = self.state_he.smass()
        self.state_he.unspecify_phase()
        # ..... saturated (gas) ...............................................
        self.p_he_4x = self.p_he_1x/(1-self.parameters['dp_he_cd'])
        self.x_he_4x = 1.0
        self.state_he.update(CoolProp.PQ_INPUTS,self.p_he_4x,self.x_he_4x)
        self.T_he_4x = self.state_he.T()
        self.i_he_4x = self.state_he.hmass()
        self.s_he_4x = self.state_he.smass()
        # ..... saturated (gas) ...............................................
        self.x_he_3x = 1.0
        self.state_he.update(CoolProp.PQ_INPUTS,self.p_he_3x,self.x_he_3x)
        self.T_he_3x = self.state_he.T()
        self.i_he_3x = self.state_he.hmass()
        self.s_he_3x = self.state_he.smass()
        # ..... gas ...........................................................
        self.T_he_3  = self.T_he_3x+self.parameters['dT_he_ev_sh']        
        self.p_he_3  = self.p_he_3x
        if  self.p_he_3 < self.state_he.p_critical()\
        and self.T_he_3 < self.state_he.T_critical():
            self.state_he.specify_phase(CoolProp.iphase_gas)
        if  self.p_he_3 < self.state_he.p_critical()\
        and self.T_he_3 > self.state_he.T_critical():
            self.state_he.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_he.update(CoolProp.PT_INPUTS,self.p_he_3,self.T_he_3)
        self.i_he_3  = self.state_he.hmass()
        self.x_he_3  = self.state_he.Q()
        self.s_he_3  = self.state_he.smass()
        self.v_he_3  = self.state_he.rhomass()**(-1)
        self.state_he.unspecify_phase()
        # ..... saturated (liquid) ............................................
        self.p_he_2x = self.p_he_3x/(1-self.parameters['dp_he_ev'])
        self.x_he_2x = 0.0
        self.state_he.update(CoolProp.PQ_INPUTS,self.p_he_2x,self.x_he_2x)
        self.T_he_2x = self.state_he.T()
        self.i_he_2x = self.state_he.hmass()
        self.s_he_2x = self.state_he.smass()
        # ..... liquid ........................................................
        self.p_he_2  = self.p_he_2x
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,self.p_he_2,
                                                        self.s_he_1)
            self.i_he_2s = self.state_he.hmass()
            self.i_he_2  = self.i_he_1\
                         +(self.i_he_2s-self.i_he_1)/self.parameters['eta_pm']
            self.state_he.update(CoolProp.HmassP_INPUTS,
                                 self.i_he_2,self.p_he_2)
            self.T_he_2  = self.state_he.T()
            self.P_pm    = 0
            self.eta_is_pm = self.parameters['eta_pm']
        self.x_he_2  = self.state_he.Q()
        self.s_he_2  = self.state_he.smass()
        # ..... gas ...........................................................
        self.p_he_4  = self.p_he_4x
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,
                                 self.p_he_4,self.s_he_3)
            self.i_he_4s = self.state_he.hmass()
            self.i_he_4  = self.i_he_3\
                         -(self.i_he_3-self.i_he_4s)\
                         * self.parameters['eta_max_ex']
            self.m_he    = 0
            self.P_ex    = 0
            self.eta_is_ex = self.parameters['eta_max_ex']
            self.ex_flag = 'na'
            self.state_he.update(CoolProp.HmassP_INPUTS,
                                 self.i_he_4,self.p_he_4)
            self.T_he_4  = self.state_he.T()
        self.s_he_4  = self.state_he.smass()
        self.x_he_4  = self.state_he.Q()
        self.v_he_4  = self.state_he.rhomass()**(-1)
        
        # --- Evaluate the operational conditions -----------------------------
        
        if self.parameters['version'] == 'operational_light':
            self.m_he    = self.find_mass_flow()
            self.P_ex    = self.m_he    * (self.i_he_3 - self.i_he_4)
            self.m_he_cs = self.find_flow_secondary_cs()
            self.m_he_hs = self.find_flow_secondary_hs()
        
        # --- Evaluate the Key Performance Indicators -------------------------
        
        self.retrieve_kpi()
    
    # ======================================================================= #
    # ==== SUB-FUNCTIONS FOR HEAT ENGINE MODEL ============================== #
    # ======================================================================= #

    def find_p(self):
        """
        Find the condensing and evaporating pressures in the SBORC based
        on the pinch-point and for the selected superheating and subcooling.
        Note that phase change is not permitted for the secondary fluid.
        """
        # --- Minimum condenser pressure (firm) -------------------------------
        T_he_1x_min         = self.T_he_cs_su \
                            + self.parameters['dT_he_cd_pp'] \
                            + self.parameters['dT_he_cd_sc']
        if T_he_1x_min > self.state_he.T_critical():
            raise ValueError('An inconsistency was detected in the cycle!\
                              In: '+os.path.join(os.path.abspath(__file__),
                              'SBORC','find_p:\
            the min. condensation temperature is above the critical point!'))
        self.state_he.update(CoolProp.QT_INPUTS,0.0,T_he_1x_min)  
        self.p_he_1x_min    = self.state_he.p()
        # --- Maximum condenser pressure (indicative) -------------------------
        T_he_1x_max         = max(self.T_he_cs_ex\
                                 +self.parameters['dT_he_cd_pp'],
                                  T_he_1x_min)
        if T_he_1x_max > self.state_he.T_critical():
            T_he_1x_max = self.state_he.T_critical() - 1
        self.state_he.update(CoolProp.QT_INPUTS,0.0,T_he_1x_max)
        self.p_he_1x_max    = self.state_he.p() 
        # --- Maximum evaporator pressure (firm) ------------------------------
        T_he_3x_max         = self.T_he_hs_su\
                            - self.parameters['dT_he_ev_pp']\
                            - self.parameters['dT_he_ev_sh']
        if T_he_3x_max >= self.state_he.T_critical() - 1: 
            T_he_3x_max = self.state_he.T_critical() - 1
        self.state_he.update(CoolProp.QT_INPUTS,1.0,T_he_3x_max)
        self.p_he_3x_max    = min(self.state_he.p(),
                                  self.state_he.p_critical()-1)
        # --- Minimum evaporator pressure (indicative) ------------------------
        T_he_3x_min = min(self.T_he_hs_ex-self.parameters['dT_he_ev_pp'],
                          T_he_3x_max)
        if T_he_3x_min >= self.state_he.T_critical() - 1:
            T_he_3x_min = self.state_he.T_critical() - 1
        self.state_he.update(CoolProp.QT_INPUTS,1.0,T_he_3x_min)
        self.p_he_3x_min = self.state_he.p()  
        # --- Find the pressure levels ----------------------------------------
        guess = self.p_he_1x_max, self.p_he_3x_min 
        try:
            p_sol = scipy.optimize.fsolve(self.resi_p, guess)
        except:
            guess = self.p_he_1x_min, self.p_he_3x_max 
            res = scipy.optimize.least_squares(self.resi_p, guess, 
                  bounds = ((self.p_he_1x_min,self.p_he_1x_min),
                            (self.p_he_3x_max,self.p_he_3x_max)),
                  ftol=1e-12,xtol=1e-12)
            p_sol = res.x
            import warnings
            warnings.warn('In: SBORC.find_p, fsolve failed.'+
                          ' least_squares was tested instead. Message: '+
                          str(res.message))
        self.resi_1, self.resi_2 = self.resi_p((p_sol[0],p_sol[1]))
        if abs(self.resi_1) > 1 or abs(self.resi_2) > 1\
        or p_sol[1] > self.p_he_3x_max:
            guess = (self.p_he_1x_min+self.p_he_1x_max)/2,\
                    (self.p_he_3x_min+self.p_he_3x_max)/2
            # print(f'p_he_1x_min:   {self.p_he_1x_min*1e-5}')
            # print(f'p_he_1x_max:   {self.p_he_1x_max*1e-5}')
            # print(f'p_he_1x_guess: {guess[0]*1e-5}')
            # print(f'p_he_3x_min:   {self.p_he_3x_min*1e-5}')
            # print(f'p_he_3x_max:   {self.p_he_3x_max*1e-5}')
            # print(f'p_he_3x_guess: {guess[1]*1e-5}')
            res = scipy.optimize.least_squares(self.resi_p, guess, 
                  bounds = ((self.p_he_1x_min,self.p_he_1x_min),
                            (self.p_he_3x_max,self.p_he_3x_max)),
                  ftol=1e-12,xtol=1e-12)
            p_sol = res.x
            import warnings
            warnings.warn('In: SBORC.find_p, fsolve did not converge.'+
                          ' least_squares was tested instead. Message: '+
                          str(res.message))
        self.resi_1, self.resi_2 = self.resi_p((p_sol[0],p_sol[1]))
        p_he_1x  = p_sol[0]
        p_he_3x  = p_sol[1]
        # print('p_1x,min: %0.3f' % (self.p_he_1x_min*1e-5))
        # print('p_1x:     %0.3f' % (     p_he_1x    *1e-5))
        # print(f'p_3x:         {p_he_3x    *1e-5}')
        # print('resi_1:   %0.3f' % (self.resi_1))
        # print('resi_2:   %0.3f' % (self.resi_2))
        return p_he_1x, p_he_3x        

    def resi_p(self,iter_var):
        """
        Residuals of pinch-point method in condenser and evaporator.  
        Used to find the condensing and evaporating pressures in the SBORC.
        
        Cells discretization: 0.25K per cell (for primary fluid)
        """
        # --- Guess values ----------------------------------------------------
        p_he_1x, p_he_3x = iter_var
        # --- Condenser -------------------------------------------------------
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_1x,0.0)
        T_he_1x = self.state_he.T()
        p_he_1  = p_he_1x
        T_he_1  = T_he_1x-self.parameters['dT_he_cd_sc']
        self.state_he.specify_phase(CoolProp.iphase_liquid)
        self.state_he.update(CoolProp.PT_INPUTS,p_he_1,T_he_1)
        h_he_1  = self.state_he.hmass()
        s_he_1  = self.state_he.smass()
        self.state_he.unspecify_phase()
        p_he_4x = p_he_1x/(1-self.parameters['dp_he_cd'])
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_4x,1.0)
        T_he_4x = self.state_he.T()
        h_he_4x = self.state_he.hmass()
        # --- Evaporator ------------------------------------------------------
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_3x,1.0)
        T_he_3x = self.state_he.T()
        h_he_3x = self.state_he.hmass()
        p_he_3  = p_he_3x
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            T_he_3  = T_he_3x+self.parameters['dT_he_ev_sh']
            # T_he_3  = self.T_he_hs_su-self.parameters['dT_he_ev_pp']
        if  p_he_3 < self.state_he.p_critical()\
        and T_he_3 < self.state_he.T_critical():
            self.state_he.specify_phase(CoolProp.iphase_gas)
        if  p_he_3 < self.state_he.p_critical()\
        and T_he_3 > self.state_he.T_critical():
            self.state_he.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_he.update(CoolProp.PT_INPUTS,p_he_3,T_he_3)
        h_he_3  = self.state_he.hmass()
        s_he_3  = self.state_he.smass()
        self.state_he.unspecify_phase()
        p_he_2x = p_he_3x/(1-self.parameters['dp_he_ev'])
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_2x,0.0)
        T_he_2x = self.state_he.T()
        h_he_2x = self.state_he.hmass()
        # --- Pump ------------------------------------------------------------
        p_he_2  = p_he_2x
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,p_he_2,s_he_1)
            h_he_2s = self.state_he.hmass()
            h_he_2  = h_he_1+(h_he_2s-h_he_1)/self.parameters['eta_pm']
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_2,p_he_2)
            T_he_2  = self.state_he.T() 
        # --- Expander --------------------------------------------------------
        p_he_4  = p_he_4x
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,p_he_4,s_he_3)
            h_he_4s = self.state_he.hmass()
            h_he_4  = h_he_3-(h_he_3-h_he_4s)*self.parameters['eta_max_ex']
        
        # --> Energy balance on evaporator ------------------------------------
        frac_hs = (h_he_3-h_he_2)/(self.i_he_hs_su-self.i_he_hs_ex)
        T_hi    = T_he_2x + self.parameters['dT_he_ev_pp']
        self.state_hs.update(CoolProp.PT_INPUTS,self.p_he_hs_ex,T_hi)
        h_hi    = self.state_hs.hmass()  
        d_hi    = frac_hs*(self.i_he_hs_su-h_hi)-(h_he_3-h_he_2x)
        T_wh    = T_he_2x
        # .......... check if pinch violated ..................................
        cell_resolution = 1.0 # K
        n_elem = max(int(1/cell_resolution*(T_he_2x-T_he_2)+1),3)
        p_he_wf_vec  = np.linspace(p_he_2,         p_he_2x,        n_elem)
        h_he_wf_vec  = np.linspace(h_he_2,         h_he_2x,        n_elem)
        p_he_hs_vec  = np.linspace(self.p_he_hs_ex,self.p_he_hs_ex,n_elem)
        h_he_hs_vec  = np.linspace(self.i_he_hs_ex,           h_hi,n_elem)
        T_he_wf_vec  = np.zeros(n_elem)
        T_he_hs_vec  = np.zeros(n_elem)
        for i, T in enumerate(T_he_wf_vec):
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_wf_vec[i],
                                                        p_he_wf_vec[i])
            self.state_hs.update(CoolProp.HmassP_INPUTS,h_he_hs_vec[i],
                                                        p_he_hs_vec[i])
            T_he_wf_vec[i] = self.state_he.T()
            T_he_hs_vec[i] = self.state_hs.T()
        real_pp = round(min(T_he_hs_vec-T_he_wf_vec),3)
        trgt_pp = self.parameters['dT_he_ev_pp']
        chck_pp = real_pp < trgt_pp
        # .......... discretized model ........................................
        if chck_pp:
            self.m_rat      = frac_hs
            cell_resolution = 0.25 # K
            # ======> NEW <====================================================
            n_el_lq         = int(1/cell_resolution*abs((T_he_2 -T_he_2x))+2)
            n_el_st         = int(1/cell_resolution*abs((T_he_2x-T_he_3x))+2)
            n_el_vp         = int(1/cell_resolution*abs((T_he_3x-T_he_3 ))+2)  
            p_he_wf_lq_vec  = np.linspace(p_he_2, p_he_2x,n_el_lq)
            p_he_wf_st_vec  = np.linspace(p_he_2x,p_he_3x,n_el_st)
            p_he_wf_vp_vec  = np.linspace(p_he_3x,p_he_3, n_el_vp)
            p_he_wf_vec     = np.concatenate((p_he_wf_lq_vec,
                                              p_he_wf_st_vec,
                                              p_he_wf_vp_vec),axis=None)
            h_he_wf_lq_vec  = np.linspace(h_he_2, h_he_2x,n_el_lq)
            h_he_wf_st_vec  = np.linspace(h_he_2x,h_he_3x,n_el_st)
            h_he_wf_vp_vec  = np.linspace(h_he_3x,h_he_3, n_el_vp)
            h_he_wf_vec     = np.concatenate((h_he_wf_lq_vec,
                                              h_he_wf_st_vec,
                                              h_he_wf_vp_vec),axis=None)
            n_elem          = len(p_he_wf_vec)
            # --- working fluid -----------------------------------------------
            T_he_wf_vec     = np.zeros(n_elem)
            for i, h in enumerate(h_he_wf_vec):
                self.state_he.update(CoolProp.HmassP_INPUTS,h,p_he_wf_vec[i])
                T_he_wf_vec[i] = self.state_he.T()
            # --- secondary fluid ---------------------------------------------
            p_he_hs_vec   = np.linspace(self.p_he_hs_ex,self.p_he_hs_su,n_elem)
            T_he_hs_vec   = np.zeros(n_elem)
            T_he_hs_vec[0]= self.T_he_hs_ex
            h_he_hs_vec   = np.zeros(n_elem)
            h_he_hs_vec[0]= self.i_he_hs_ex
            for i in range(len(h_he_hs_vec)-1):
                h_he_hs_vec[i+1] = h_he_hs_vec[i]\
                                 +(h_he_wf_vec[i+1]-h_he_wf_vec[i])/self.m_rat
                self.state_hs.update(CoolProp.HmassP_INPUTS,h_he_hs_vec[i+1],
                                                            p_he_hs_vec[i+1])
                T_he_hs_vec[i+1] = self.state_hs.T()
            # --- find pinch --------------------------------------------------
            imin = np.argmin(T_he_hs_vec-T_he_wf_vec)
            T_hi = T_he_hs_vec[imin]
            T_wh = T_he_wf_vec[imin]
            d_hi = T_he_hs_vec[imin]-T_he_wf_vec[imin]\
                                    -self.parameters['dT_he_ev_pp']
        # --- Energy balance on condenser -------------------------------------
        frac_cs = (h_he_4-h_he_1)/(self.i_he_cs_ex-self.i_he_cs_su)
        T_ci    = T_he_4x-self.parameters['dT_he_cd_pp']
        self.state_cs.update(CoolProp.PT_INPUTS,self.p_he_cs_ex,T_ci)
        h_ci    = self.state_cs.hmass()  
        d_ci    = frac_cs*(h_ci-self.i_he_cs_su)-(h_he_4x-h_he_1)
        if p_he_1x < self.p_he_1x_min: # pinch at the sink inlet
            T_ci    = T_he_1-self.parameters['v']
            self.state_cs.update(CoolProp.PT_INPUTS,self.p_he_cs_su,T_ci)
            h_ci    = self.state_cs.hmass()  
            d_ci    = h_ci-self.i_he_cs_su
        # --- Final balance ---------------------------------------------------
        out     = d_ci, d_hi
        self.T_ci = T_ci
        self.T_hi = T_hi
        self.T_wh = T_wh
        return out

    def find_mass_flow(self):
        """
        Find the mass flow rate for the selected operational conditions.
        """
        gue     = 1
        m_he    = scipy.optimize.fsolve(self.resi_mass_flow, gue)[0]
        if m_he < 0:    m_he = 0
        return m_he

    def resi_mass_flow(self,iter_var):
        """
        Residuals of different parameters to control the heat engine.
        """
        # --- Evaluate the massflow -------------------------------------------
        cheat_m_he     = iter_var
        # --- Check if match the demanded power input -------------------------
        if  self.parameters['mode'] == 'power':
            w_in_he_en   = self.i_he_2   - self.i_he_1
            W_in_he_en   = cheat_m_he    * w_in_he_en
            w_out_he_en  = self.i_he_3   - self.i_he_4
            W_out_he_en  = cheat_m_he    * w_out_he_en
            cheat_P_he   = W_out_he_en   - W_in_he_en
            resi         = self.P_he     - cheat_P_he
        # --- Check if match the maximum heat output --------------------------
        if  self.parameters['mode'] == 'sink':
            q_available  = self.i_he_cs_ex  - self.i_he_cs_su
            Q_available  = self.m_he_cs     * q_available
            q_required   = self.i_he_4      - self.i_he_1
            m_required   = Q_available/q_required
            resi         = cheat_m_he       - m_required
        # --- Check if match the maximum heat input ---------------------------
        if  self.parameters['mode'] == 'source':
            q_available  = self.i_he_hs_su  - self.i_he_hs_ex
            Q_available  = self.m_he_hs     * q_available
            q_required   = self.i_he_3      - self.i_he_2
            m_required   = Q_available/q_required
            resi         = cheat_m_he       - m_required
        return resi

    def find_flow_secondary_hs(self):
        """
        Finds the mass flow rate of the hot secondary fluid to match the 
        specified temperature profile based on the refrigerant flow rate.
        """
        q_available  = self.i_he_3      - self.i_he_2
        Q_available  = self.m_he        * q_available
        q_required   = self.i_he_hs_su  - self.i_he_hs_ex
        if q_required > 0:  m_required = Q_available/q_required
        else:               m_required = 0
        return m_required
    
    def find_flow_secondary_cs(self):
        """
        Finds the mass flow rate of the cold secondary fluid to match the 
        specified temperature profile based on the refrigerant flow rate.
        """
        q_available  = self.i_he_4     - self.i_he_1
        Q_available  = self.m_he       * q_available
        q_required   = self.i_he_cs_ex - self.i_he_cs_su
        if q_required > 0:  m_required = Q_available/q_required
        else:               m_required = 0
        return m_required

    
    def check_consistency(self):
        """
        Check that the results are consistent.
        """
        if self.s_he_4x > self.s_he_4:                  self.error = True
        if self.s_he_1x < self.s_he_1:                  self.error = True
        if self.s_he_2x < self.s_he_2:                  self.error = True
        if self.s_he_3x > self.s_he_3:                  self.error = True
        if self.s_he_4  < self.s_he_3:                  self.error = True
        if self.s_he_3x < self.s_he_2x:                 self.error = True
        if self.s_he_4x < self.s_he_1x:                 self.error = True
        if self.T_he_2  > self.T_he_hs_ex-0.95*self.parameters['dT_he_ev_pp']:
            self.error = True
        if self.T_he_3  > self.T_he_hs_su-0.95*self.parameters['dT_he_ev_pp']:
            self.error = True
        if self.T_he_4  < self.T_he_cs_ex+0.95*self.parameters['dT_he_cd_pp']: 
            self.error = True
        if self.T_he_1  < self.T_he_cs_su+0.95*self.parameters['dT_he_cd_pp']: 
            self.error = True
        if self.eta_he_cyclen < 0:                      self.error = True
        if self.T_hi-self.T_wh    < 0.95*self.parameters['dT_he_ev_pp']:  
            self.error = True
        if self.T_he_4x-self.T_ci < 0.95*self.parameters['dT_he_cd_pp']:    
            self.error = True  
        
        if self.p_he_1x < 0.99*self.p_he_1x_min: self.error = True
        if self.p_he_3x > 1.00*self.p_he_3x_max: self.error = True
        
        if abs(self.resi_1) > 1e-2: self.error = True
        if abs(self.resi_2) > 1e-2: self.error = True
        
        if not bool(self.parameters['wet_ex']):
            p_vec = np.linspace(self.p_he_3,self.p_he_4)
            for i, p in enumerate(p_vec):
                self.state_he.update(CoolProp.PSmass_INPUTS,p,self.s_he_3) 
                h_is = self.state_he.hmass()
                h    = self.i_he_3-self.eta_is_ex*(self.i_he_3-h_is)
                self.state_he.update(CoolProp.HmassP_INPUTS,h,p)
                if 0 < self.state_he.Q() < 1:
                    self.error = True
                if  self.error: break
    
    def retrieve_kpi(self):
        """
        Retrieve the KPI's of the SBORC
        """        
        
        if self.options['exergy']: self.compute_exergy()
        
        # --- Basic quantities ------------------------------------------------
        self.w_in_he_en = self.i_he_2 -self.i_he_1 
        self.q_in_he_en = self.i_he_3 -self.i_he_2
        self.w_out_he_en= self.i_he_3 -self.i_he_4 
        self.q_out_he_en= self.i_he_4 -self.i_he_1 
        self.w_net_he_en= self.w_out_he_en-self.w_in_he_en
        self.f_he_hs    = self.q_in_he_en /(self.i_he_hs_su-self.i_he_hs_ex)
        self.f_he_cs    = self.q_out_he_en/(self.i_he_cs_ex-self.i_he_cs_su)
        
        if self.options['exergy']: 
            self.w_in_he_ex = self.e_he_2 -self.e_he_1 
            self.q_in_he_ex = self.e_he_3 -self.e_he_2
            self.w_out_he_ex= self.e_he_3 -self.e_he_4 
            self.q_out_he_ex= self.e_he_4 -self.e_he_1 
            self.w_net_he_ex= self.w_out_he_ex-self.w_in_he_ex
        
        # --- Efficiencies ----------------------------------------------------
        self.eta_he_cyclen = self.w_net_he_en/self.q_in_he_en
        
        if self.options['exergy']: 
            self.eta_he_evapen = self.f_he_hs**(-1)*self.q_in_he_en\
                               / (self.i_he_hs_su-self.i_he_hs_ex)
            self.eta_he_srcen  = (self.i_he_hs_su-self.i_he_hs_ex) \
                               / (self.i_he_hs_su-self.i_he_hs_ref)
            self.eta_he_genven = self.eta_he_evapen*self.eta_he_srcen
            self.eta_he_toten  = self.eta_he_cyclen*self.eta_he_genven
            self.eta_he_cyclex = self.w_net_he_en/self.q_in_he_ex
            self.eta_he_pumpex = self.w_in_he_ex /self.w_in_he_en
            self.eta_he_expex  = self.w_out_he_en/self.w_out_he_ex
            self.eta_he_evapex = self.f_he_hs**(-1) * self.q_in_he_ex\
                               / (self.e_he_hs_su-self.e_he_hs_ex)
            self.eta_he_srcex  = (self.e_he_hs_su-self.e_he_hs_ex) \
                               / (self.e_he_hs_su-self.e_he_hs_ref)
            self.eta_he_genvex = self.eta_he_evapex*self.eta_he_srcex  
            self.eta_he_condex = self.f_he_cs*(self.e_he_cs_ex\
                               - self.e_he_cs_su)/self.q_out_he_ex
            self.eta_he_totex  = self.eta_he_cyclex*self.eta_he_genvex
        
        # --- Losses ----------------------------------------------------------
        if self.options['exergy']: 
            self.loss_he_srcen  = self.f_he_hs \
                                * (self.i_he_hs_ex-self.i_he_hs_ref)
            self.loss_he_sinken = self.f_he_cs \
                                * (self.i_he_cs_ex-self.i_he_cs_su)
            self.pow_he_supplen = self.f_he_hs \
                                * (self.i_he_hs_su-self.i_he_hs_ref)
            self.loss_he_pumpex = self.w_in_he_en-self.w_in_he_ex
            self.loss_he_expex  = self.w_out_he_ex-self.w_out_he_en
            self.loss_he_evapex = self.f_he_hs*(self.e_he_hs_su \
                                - self.e_he_hs_ex) - self.q_in_he_ex
            self.loss_he_srcex  = self.f_he_hs*(self.e_he_hs_ex \
                                - self.e_he_hs_ref)
            self.loss_he_condex = self.q_out_he_ex - self.f_he_cs \
                                * (self.e_he_cs_ex - self.e_he_cs_su)
            self.loss_he_sinkex = self.f_he_cs*(self.e_he_cs_ex \
                                - self.e_he_cs_su)
            self.pow_he_supplex = self.f_he_hs*(self.e_he_hs_su \
                                - self.e_he_hs_ref)
        
        # --- Carnot and Lorentz efficiencies ---------------------------------
        T_cs_logAv =(self.T_he_cs_ex-self.T_he_cs_su)\
             /np.log(self.T_he_cs_ex/self.T_he_cs_su)
        T_hs_logAv =(self.T_he_hs_su-self.T_he_hs_ex)\
             /np.log(self.T_he_hs_su/self.T_he_hs_ex)
        ETA_lorenz =(     T_hs_logAv-     T_cs_logAv)/     T_hs_logAv
        ETA_carnot =(self.T_he_hs_su-self.T_he_cs_su)/self.T_he_hs_su
        self.eta_he_lorenz = ETA_lorenz
        self.eta_he_carnot = ETA_carnot
        self.psi_he_lorenz = self.eta_he_cyclen/ETA_lorenz
        self.psi_he_carnot = self.eta_he_cyclen/ETA_carnot
        
        # --- Preliminary design information ----------------------------------
        self.vol_coef       = self.v_he_3/self.w_out_he_en
        self.exp_ratio      = self.p_he_3/self.p_he_4
        self.vol_ratio      = self.v_he_4/self.v_he_3
        if self.parameters['version'] == 'operational_light':
            # --- Multiply the specific quantities by the mass flow -----------
            self.Q_in_he_en = self.m_he*self.q_in_he_en
            self.W_in_he_en = self.m_he*self.w_in_he_en
            self.Q_out_he_en= self.m_he*self.q_out_he_en
            self.W_out_he_en= self.m_he*self.w_out_he_en
            self.W_net_he_en= self.m_he*self.w_net_he_en  
    
    def export_states(self):
        """
        Export the thermodynamic states of the SBORC necessary to draw
        the T-s and p-h diagrams.
        """
        self.p_he = self.p_he_1,self.p_he_1x,self.p_he_2,self.p_he_2x,\
                    self.p_he_3,self.p_he_3x,self.p_he_4,self.p_he_4x
        self.T_he = self.T_he_1,self.T_he_1x,self.T_he_2,self.T_he_2x,\
                    self.T_he_3,self.T_he_3x,self.T_he_4,self.T_he_4x
        self.i_he = self.i_he_1,self.i_he_1x,self.i_he_2,self.i_he_2x,\
                    self.i_he_3,self.i_he_3x,self.i_he_4,self.i_he_4x
        self.s_he = self.s_he_1,self.s_he_1x,self.s_he_2,self.s_he_2x,\
                    self.s_he_3,self.s_he_3x,self.s_he_4,self.s_he_4x
        self.x_he = self.x_he_1,self.x_he_1x,self.x_he_2,self.x_he_2x,\
                    self.x_he_3,self.x_he_3x,self.x_he_4,self.x_he_4x
        if self.options['exergy']: 
            self.e_he = self.e_he_1,self.e_he_1x,self.e_he_2,self.e_he_2x,\
                        self.e_he_3,self.e_he_3x,self.e_he_4,self.e_he_4x
        else:
            self.e_he = 0,0,0,0,0,0,0,0
        out = self.p_he,self.T_he,self.i_he,self.s_he,self.x_he,self.e_he
        return out
    
    def compute_ref_states(self):
        """
        Compute the reference state.
        """
        self.T_he_0 = self.parameters['T_0']
        self.p_he_0 = self.parameters['p_0']
        self.T_ref  = self.parameters['T_ref']
        self.p_ref  = self.parameters['p_ref']
        
        # Saturated liquid conditions must be taken as reference if the 
        # fluid can condense at T_0: the exergy is therefore 0 for any 
        # state of saturated vapour at T_0. This is because the exergy 
        # variation associated with a phase change at T_0 is always 0.
        # Refer to "THERMAL POWER PLANTS: Energetic and exergetic 
        # approaches" pp.13 for more details.
        self.state_he.update(CoolProp.QT_INPUTS,0,self.T_he_0)
        self.i_he_0 = self.state_he.hmass()
        self.s_he_0 = self.state_he.smass()
        if np.isnan(self.i_he_0): # the fluid above the critical point
            self.state_he.update(CoolProp.PT_INPUTS,self.p_he_0,self.T_he_0)
            self.i_he_0 = self.state_he.hmass()
            self.s_he_0 = self.state_he.smass()
        
        # Saturated liquid conditions must be taken as reference if the 
        # fluid can condense at T_0: the exergy is therefore 0 for any 
        # state of saturated vapour at T_0. This is because the exergy 
        # variation associated with a phase change at T_0 is always 0.
        # Refer to "THERMAL POWER PLANTS: Energetic and exergetic 
        # approaches" pp.13 for more details.
        self.state_cs.update(CoolProp.QT_INPUTS,0,self.T_he_0)
        self.i_he_cs_0 = self.state_cs.hmass()
        self.s_he_cs_0 = self.state_cs.smass()   
        if np.isnan(self.i_he_cs_0): # the fluid above the critical point
            self.state_cs.update(CoolProp.PT_INPUTS,self.p_he_0,self.T_he_0)
            self.i_he_cs_0 = self.state_cs.hmass()
            self.s_he_cs_0 = self.state_cs.smass()
        
        if self.options['exergy']:
            # Saturated liquid conditions must be taken as reference if the 
            # fluid can condense at T_0: the exergy is therefore 0 for any 
            # state of saturated vapour at T_0. This is because the exergy 
            # variation associated with a phase change at T_0 is always 0.
            # Refer to "THERMAL POWER PLANTS: Energetic and exergetic 
            # approaches" pp.13 for more details.
            self.state_hs.update(CoolProp.QT_INPUTS,0,self.T_he_0)
            self.i_he_hs_0 = self.state_hs.hmass()
            self.s_he_hs_0 = self.state_hs.smass() 
            if np.isnan(self.i_he_hs_0): # the fluid above the critical point  
                self.state_hs.update(CoolProp.PT_INPUTS,self.p_he_0,self.T_he_0)
                self.i_he_hs_0 = self.state_hs.hmass()
                self.s_he_hs_0 = self.state_hs.smass() 
        
    def compute_exergy(self):
        """
        Compute the exergy at each state of the SBORC.
        """
        self.compute_ref_states()
        
        # 1 - Exergy of the cycle ---------------------------------------------
        
        self.e_he_1  = (self.i_he_1 -self.i_he_0)\
          -self.T_he_0*(self.s_he_1 -self.s_he_0)
        self.e_he_1x = (self.i_he_1x-self.i_he_0)\
          -self.T_he_0*(self.s_he_1x-self.s_he_0)
        self.e_he_2  = (self.i_he_2 -self.i_he_0)\
          -self.T_he_0*(self.s_he_2 -self.s_he_0)
        self.e_he_2x = (self.i_he_2x-self.i_he_0)\
          -self.T_he_0*(self.s_he_2x-self.s_he_0)
        self.e_he_3  = (self.i_he_3 -self.i_he_0)\
          -self.T_he_0*(self.s_he_3 -self.s_he_0)
        self.e_he_3x = (self.i_he_3x-self.i_he_0)\
          -self.T_he_0*(self.s_he_3x-self.s_he_0)
        self.e_he_4  = (self.i_he_4 -self.i_he_0)\
          -self.T_he_0*(self.s_he_4 -self.s_he_0)
        self.e_he_4x = (self.i_he_4x-self.i_he_0)\
          -self.T_he_0*(self.s_he_4x-self.s_he_0)
        
        if round(self.T_he_1 ,2)==round(self.T_he_0,2): self.e_he_1  = 0
        if round(self.T_he_1x,2)==round(self.T_he_0,2): self.e_he_1x = 0
        if round(self.T_he_2 ,2)==round(self.T_he_0,2): self.e_he_2  = 0
        if round(self.T_he_2x,2)==round(self.T_he_0,2): self.e_he_2x = 0
        if round(self.T_he_3 ,2)==round(self.T_he_0,2): self.e_he_3  = 0
        if round(self.T_he_3x,2)==round(self.T_he_0,2): self.e_he_3x = 0
        if round(self.T_he_4 ,2)==round(self.T_he_0,2): self.e_he_4  = 0
        if round(self.T_he_4x,2)==round(self.T_he_0,2): self.e_he_4x = 0
        
        # 2 - Exergy of the boundaries ----------------------------------------
        
        self.compute_exergy_boundaries()

    def compute_exergy_boundaries(self):
        """
        Compute the exergy of the boundaries.
        """
        
        # 1 - Exergy of the sink boundary -------------------------------------
        
        self.e_he_cs_su = (self.i_he_cs_su -self.i_he_cs_0)\
             -self.T_he_0*(self.s_he_cs_su -self.s_he_cs_0)
        self.e_he_cs_ex = (self.i_he_cs_ex -self.i_he_cs_0)\
             -self.T_he_0*(self.s_he_cs_ex -self.s_he_cs_0)

        if round(self.T_he_cs_su,2)==round(self.T_he_0,2): self.e_he_cs_su = 0
        if round(self.T_he_cs_ex,2)==round(self.T_he_0,2): self.e_he_cs_ex = 0
        
        self.state_cs.update(CoolProp.PT_INPUTS,self.p_ref,self.T_ref)
        self.i_he_cs_ref  = self.state_cs.hmass()
        self.s_he_cs_ref  = self.state_cs.smass()
        self.e_he_cs_ref  = (self.i_he_cs_ref -self.i_he_cs_0)\
               -self.T_he_0*(self.s_he_cs_ref -self.s_he_cs_0)

        if round(self.T_ref     ,2)==round(self.T_he_0,2): self.e_he_cs_ref= 0
        
        # 2 - Exergy of the source boundary -----------------------------------

        if self.options['exergy']:
            self.e_he_hs_su = (self.i_he_hs_su -self.i_he_hs_0)\
                 -self.T_he_0*(self.s_he_hs_su -self.s_he_hs_0)
            self.e_he_hs_ex = (self.i_he_hs_ex -self.i_he_hs_0)\
                 -self.T_he_0*(self.s_he_hs_ex -self.s_he_hs_0)
    
            if round(self.T_he_hs_su,2)==round(self.T_he_0,2): 
                self.e_he_hs_su = 0
            if round(self.T_he_hs_ex,2)==round(self.T_he_0,2): 
                self.e_he_hs_ex = 0
    
            self.state_hs.update(CoolProp.PT_INPUTS,self.p_ref,self.T_ref)
            self.i_he_hs_ref  = self.state_hs.hmass()
            self.s_he_hs_ref  = self.state_hs.smass()
            self.e_he_hs_ref  = (self.i_he_hs_ref -self.i_he_hs_0)\
                   -self.T_he_0*(self.s_he_hs_ref -self.s_he_hs_0)
            
            if round(self.T_ref     ,2)==round(self.T_he_0,2): 
                self.e_he_hs_ref= 0

#%% MODEL: SRORC

class SRORC(SBORC):
    """
    Class for the simulation of Subcritical Recuperated Organic Rankine Cycles.

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (heat engine boundaries and operational 
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
    # ==== HEAT ENGINE MODEL ================================================ #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the hot source
            * 'sink':   available heat flow rate at the cold sink
            * 'power':  electrical power produced by the heat engine
            
        """        
        
        # --- Cycle computation -----------------------------------------------
        
        # ..... saturated (liquid) ............................................
        self.p_he_1x, self.p_he_3x = self.find_p()
        self.x_he_1x = 0.0
        self.state_he.update(CoolProp.PQ_INPUTS,self.p_he_1x,self.x_he_1x)
        self.T_he_1x = self.state_he.T()
        self.i_he_1x = self.state_he.hmass()
        self.s_he_1x = self.state_he.smass()
        # ..... liquid ........................................................
        self.p_he_1  = self.p_he_1x
        self.T_he_1  = self.T_he_1x-self.parameters['dT_he_cd_sc']
        self.state_he.specify_phase(CoolProp.iphase_liquid)
        self.state_he.update(CoolProp.PT_INPUTS,self.p_he_1,self.T_he_1)
        self.i_he_1  = self.state_he.hmass()       
        self.x_he_1  = self.state_he.Q()
        self.s_he_1  = self.state_he.smass()
        self.state_he.unspecify_phase()
        # ..... saturated (gas) ...............................................
        self.p_he_4x = self.p_he_1x/(1-self.parameters['dp_he_cd'])
        self.x_he_4x = 1.0
        self.state_he.update(CoolProp.PQ_INPUTS,self.p_he_4x,self.x_he_4x)
        self.T_he_4x = self.state_he.T()
        self.i_he_4x = self.state_he.hmass()
        self.s_he_4x = self.state_he.smass()
        # ..... saturated (gas) ...............................................
        self.x_he_3x = 1.0
        self.state_he.update(CoolProp.PQ_INPUTS,self.p_he_3x,self.x_he_3x)
        self.T_he_3x = self.state_he.T()
        self.i_he_3x = self.state_he.hmass()
        self.s_he_3x = self.state_he.smass()        
        # ..... gas ...........................................................
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.T_he_3  = self.T_he_3x+self.parameters['dT_he_ev_sh']
        self.p_he_3  = self.p_he_3x
        if  self.p_he_3 < self.state_he.p_critical()\
        and self.T_he_3 < self.state_he.T_critical():
            self.state_he.specify_phase(CoolProp.iphase_gas)
        if  self.p_he_3 < self.state_he.p_critical()\
        and self.T_he_3 > self.state_he.T_critical():
            self.state_he.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_he.update(CoolProp.PT_INPUTS,self.p_he_3,self.T_he_3)
        self.i_he_3  = self.state_he.hmass()
        self.x_he_3  = self.state_he.Q()
        self.s_he_3  = self.state_he.smass()
        self.v_he_3  = self.state_he.rhomass()**(-1)
        self.state_he.unspecify_phase()
        # ..... saturated (liquid) ............................................
        self.p_he_2x = self.p_he_3x/(1-self.parameters['dp_he_ev'])
        self.x_he_2x = 0.0
        self.state_he.update(CoolProp.PQ_INPUTS,self.p_he_2x,self.x_he_2x)
        self.T_he_2x = self.state_he.T()
        self.i_he_2x = self.state_he.hmass()
        self.s_he_2x = self.state_he.smass()
        # ..... liquid ........................................................
        self.p_he_2  = self.p_he_2x/(1-self.parameters['dp_he_rg_lq'])
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,self.p_he_2,
                                                        self.s_he_1)
            self.i_he_2s = self.state_he.hmass()
            self.i_he_2  = self.i_he_1\
                         +(self.i_he_2s-self.i_he_1)/self.parameters['eta_pm']
            self.state_he.update(CoolProp.HmassP_INPUTS,
                                 self.i_he_2,self.p_he_2)
            self.T_he_2  = self.state_he.T()
            self.P_pm    = 0
            self.eta_is_pm = self.parameters['eta_pm']
        self.x_he_2  = self.state_he.Q()
        self.s_he_2  = self.state_he.smass()
        # ..... gas ...........................................................
        self.p_he_4  = self.p_he_4x/(1-self.parameters['dp_he_rg_vp'])
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,
                                 self.p_he_4,self.s_he_3)
            self.i_he_4s = self.state_he.hmass()
            self.i_he_4  = self.i_he_3\
                         -(self.i_he_3-self.i_he_4s)\
                         * self.parameters['eta_max_ex']
            self.m_he    = 0
            self.P_ex    = 0
            self.eta_is_ex = self.parameters['eta_max_ex']
            self.ex_flag = 'na'
            self.state_he.update(CoolProp.HmassP_INPUTS,
                                 self.i_he_4,self.p_he_4)
            self.T_he_4  = self.state_he.T()
        self.s_he_4  = self.state_he.smass()
        self.x_he_4  = self.state_he.Q()
        self.v_he_4  = self.state_he.rhomass()**(-1)
        # ..... gas (potentially saturated) ...................................
        self.p_he_4r = self.p_he_4x
        self.T_he_4r = self.T_he_4\
                     -(self.T_he_4-self.T_he_2)*self.parameters['epsilon']
        if self.T_he_4r == self.T_he_4x:           
            self.state_he.update(CoolProp.QT_INPUTS,1.0,self.T_he_4r)
        else:
            if  self.p_he_4r < self.state_he.p_critical()\
            and self.T_he_4r < self.state_he.T_critical():
                self.state_he.specify_phase(CoolProp.iphase_gas)
            if  self.p_he_4r < self.state_he.p_critical()\
            and self.T_he_4r > self.state_he.T_critical():
                self.state_he.specify_phase(CoolProp.iphase_supercritical_gas)
            self.state_he.update(CoolProp.PT_INPUTS,self.p_he_4r,self.T_he_4r)
            self.state_he.unspecify_phase()
        self.i_he_4r = self.state_he.hmass() 
        self.s_he_4r = self.state_he.smass() 
        self.x_he_4r = self.state_he.Q()
        # ..... liquid (potentially saturated) ................................
        self.i_he_2r = self.i_he_2+(self.i_he_4-self.i_he_4r)
        self.p_he_2r = self.p_he_2x
        self.state_he.update(CoolProp.HmassP_INPUTS,self.i_he_2r,self.p_he_2r)
        self.T_he_2r = self.state_he.T()
        self.s_he_2r = self.state_he.smass()
        self.x_he_2r = self.state_he.Q()
        
        # --- Evaluate the operational conditions -----------------------------
        
        if self.parameters['version'] == 'operational_light':
            
            self.m_he    = self.find_mass_flow()
            self.P_ex    = self.m_he    * (self.i_he_3 - self.i_he_4)
            self.m_he_cs = self.find_flow_secondary_cs()
            self.m_he_hs = self.find_flow_secondary_hs()
        
        # --- Evaluate the Key Performance Indicators -------------------------
        
        self.retrieve_kpi()
    
    # ======================================================================= #
    # ==== SUB-FUNCTIONS FOR HEAT ENGINE MODEL ============================== #
    # ======================================================================= #

    def resi_p(self,iter_var):
        """
        Residuals of pinch-point method in condenser and evaporator.  
        Used to find the condensing and evaporating pressures in the SRORC.
        
        Cells discretization: 0.25K per cell (for primary fluid)
        """
        # --- Guess values ----------------------------------------------------
        p_he_1x, p_he_3x = iter_var
        # --- Condenser -------------------------------------------------------        
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_1x,0.0)
        T_he_1x = self.state_he.T()
        p_he_1  = p_he_1x
        T_he_1  = T_he_1x-self.parameters['dT_he_cd_sc']
        self.state_he.specify_phase(CoolProp.iphase_liquid)
        self.state_he.update(CoolProp.PT_INPUTS,p_he_1,T_he_1)
        h_he_1  = self.state_he.hmass()
        s_he_1  = self.state_he.smass()
        self.state_he.unspecify_phase()
        p_he_4x = p_he_1x/(1-self.parameters['dp_he_cd'])
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_4x,1.0)
        T_he_4x = self.state_he.T()
        h_he_4x = self.state_he.hmass()
        # --- Evaporator ------------------------------------------------------
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_3x,1.0)
        T_he_3x = self.state_he.T()
        h_he_3x = self.state_he.hmass()
        p_he_3  = p_he_3x
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            T_he_3  = T_he_3x+self.parameters['dT_he_ev_sh']
        if  p_he_3 < self.state_he.p_critical()\
        and T_he_3 < self.state_he.T_critical():
            self.state_he.specify_phase(CoolProp.iphase_gas)
        if  p_he_3 < self.state_he.p_critical()\
        and T_he_3 > self.state_he.T_critical():
            self.state_he.specify_phase(CoolProp.iphase_supercritical_gas)
        self.state_he.update(CoolProp.PT_INPUTS,p_he_3,T_he_3)
        h_he_3  = self.state_he.hmass()
        s_he_3  = self.state_he.smass()
        self.state_he.unspecify_phase()
        p_he_2x = p_he_3x/(1-self.parameters['dp_he_ev'])
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_2x,0.0)
        T_he_2x = self.state_he.T()
        h_he_2x = self.state_he.hmass()
        # --- Pump ------------------------------------------------------------                
        p_he_2  = p_he_2x/(1-self.parameters['dp_he_rg_lq'])
        if self.parameters['version'] == 'thermodynamic_full' \
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,p_he_2,s_he_1)
            h_he_2s = self.state_he.hmass()
            h_he_2  = h_he_1+(h_he_2s-h_he_1)/self.parameters['eta_pm']
            self.state_he.specify_phase(CoolProp.iphase_liquid)
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_2,p_he_2)
            T_he_2  = self.state_he.T()
            self.state_he.unspecify_phase()
        # --- Expander --------------------------------------------------------
        p_he_4  = p_he_4x/(1-self.parameters['dp_he_rg_vp'])
        if self.parameters['version'] == 'thermodynamic_full' \
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,p_he_4,s_he_3)
            h_he_4s = self.state_he.hmass()
            h_he_4  = h_he_3-(h_he_3-h_he_4s)*self.parameters['eta_max_ex']
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_4,p_he_4)
            T_he_4  = self.state_he.T()
        # --- Regenerator -----------------------------------------------------
        p_he_4r = p_he_4x
        T_he_4r = T_he_4-(T_he_4-T_he_2)*self.parameters['epsilon']
        if T_he_4r == T_he_4x:           
            self.state_he.update(CoolProp.QT_INPUTS,1.0,T_he_4r)
        else:
            self.state_he.update(CoolProp.PT_INPUTS,p_he_4r,T_he_4r)
        h_he_4r = self.state_he.hmass()       
        h_he_2r = h_he_2 + (h_he_4-h_he_4r)
        p_he_2r = p_he_2x
        self.state_he.update(CoolProp.HmassP_INPUTS,h_he_2r,p_he_2r)
        T_he_2r = self.state_he.T()
        if h_he_2r > h_he_2x:
            # only for energy balance!
            h_he_2x = h_he_2r
        
        # --> Energy balance on evaporator ------------------------------------
        frac_hs = (h_he_3-h_he_2r)/(self.i_he_hs_su-self.i_he_hs_ex)
        T_hi    = T_he_2x + self.parameters['dT_he_ev_pp']
        self.state_hs.update(CoolProp.PT_INPUTS,self.p_he_hs_ex,T_hi)
        h_hi = self.state_hs.hmass() 
        d_hi    = frac_hs*(self.i_he_hs_su-h_hi)-(h_he_3-h_he_2x)
        T_wh    = T_he_2x
        # .......... check if pinch violated ..................................
        cell_resolution = 1.0 # K
        n_elem = max(int(1/cell_resolution*(T_he_2x-T_he_2r)+1),3)
        p_he_wf_vec  = np.linspace(p_he_2r,        p_he_2x,        n_elem)
        h_he_wf_vec  = np.linspace(h_he_2r,        h_he_2x,        n_elem)
        p_he_hs_vec  = np.linspace(self.p_he_hs_ex,self.p_he_hs_ex,n_elem)
        h_he_hs_vec  = np.linspace(self.i_he_hs_ex,           h_hi,n_elem)
        T_he_wf_vec  = np.zeros(n_elem)
        T_he_hs_vec  = np.zeros(n_elem)
        for i, T in enumerate(T_he_wf_vec):
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_wf_vec[i],
                                                        p_he_wf_vec[i])
            self.state_hs.update(CoolProp.HmassP_INPUTS,h_he_hs_vec[i],
                                                        p_he_hs_vec[i])
            T_he_wf_vec[i] = self.state_he.T()
            T_he_hs_vec[i] = self.state_hs.T()
        real_pp = round(min(T_he_hs_vec-T_he_wf_vec),3)
        trgt_pp = self.parameters['dT_he_ev_pp']
        chck_pp = real_pp < trgt_pp
        # .......... discretized model ........................................
        if chck_pp:
            self.m_rat      = frac_hs
            cell_resolution = 0.25 # K
            # ======> NEW <====================================================
            n_el_lq         = int(1/cell_resolution*abs((T_he_2r-T_he_2x))+2)
            n_el_st         = int(1/cell_resolution*abs((T_he_2x-T_he_3x))+2)
            n_el_vp         = int(1/cell_resolution*abs((T_he_3x-T_he_3 ))+2)  
            p_he_wf_lq_vec  = np.linspace(p_he_2r,p_he_2x,n_el_lq)
            p_he_wf_st_vec  = np.linspace(p_he_2x,p_he_3x,n_el_st)
            p_he_wf_vp_vec  = np.linspace(p_he_3x,p_he_3, n_el_vp)
            p_he_wf_vec     = np.concatenate((p_he_wf_lq_vec,
                                              p_he_wf_st_vec,
                                              p_he_wf_vp_vec),axis=None)
            h_he_wf_lq_vec  = np.linspace(h_he_2r,h_he_2x,n_el_lq)
            h_he_wf_st_vec  = np.linspace(h_he_2x,h_he_3x,n_el_st)
            h_he_wf_vp_vec  = np.linspace(h_he_3x,h_he_3, n_el_vp)
            h_he_wf_vec     = np.concatenate((h_he_wf_lq_vec,
                                              h_he_wf_st_vec,
                                              h_he_wf_vp_vec),axis=None)
            n_elem          = len(p_he_wf_vec)
            # --- working fluid -----------------------------------------------
            T_he_wf_vec     = np.zeros(n_elem)
            for i, h in enumerate(h_he_wf_vec):
                self.state_he.update(CoolProp.HmassP_INPUTS,h,p_he_wf_vec[i])
                T_he_wf_vec[i] = self.state_he.T()
            # --- secondary fluid ---------------------------------------------
            p_he_hs_vec   = np.linspace(self.p_he_hs_ex,self.p_he_hs_su,n_elem)
            T_he_hs_vec   = np.zeros(n_elem)
            T_he_hs_vec[0]= self.T_he_hs_ex
            h_he_hs_vec   = np.zeros(n_elem)
            h_he_hs_vec[0]= self.i_he_hs_ex
            for i in range(len(h_he_hs_vec)-1):
                h_he_hs_vec[i+1] = h_he_hs_vec[i]\
                                  +(h_he_wf_vec[i+1]-h_he_wf_vec[i])/self.m_rat
                self.state_hs.update(CoolProp.HmassP_INPUTS,h_he_hs_vec[i+1],
                                                            p_he_hs_vec[i+1])
                T_he_hs_vec[i+1] = self.state_hs.T()
            # --- find pinch --------------------------------------------------
            imin = np.argmin(T_he_hs_vec-T_he_wf_vec)
            T_hi = T_he_hs_vec[imin]
            T_wh = T_he_wf_vec[imin]
            d_hi = T_he_hs_vec[imin]-T_he_wf_vec[imin]\
                                    -self.parameters['dT_he_ev_pp']
        
        # --> Energy balance on condenser -------------------------------------
        frac_cs = (h_he_4r-h_he_1)/(self.i_he_cs_ex-self.i_he_cs_su)
        T_ci    = T_he_4x-self.parameters['dT_he_cd_pp']
        self.state_cs.update(CoolProp.PT_INPUTS,self.p_he_cs_ex,T_ci)
        h_ci    = self.state_cs.hmass()  
        d_ci    = frac_cs*(h_ci-self.i_he_cs_su)-(h_he_4x-h_he_1)
        if p_he_1x < self.p_he_1x_min: # pinch at the sink inlet
            T_ci    = T_he_1-self.parameters['dT_he_cd_pp']
            self.state_cs.update(CoolProp.PT_INPUTS,self.p_he_cs_su,T_ci)
            h_ci    = self.state_cs.hmass()  
            d_ci    = h_ci-self.i_he_cs_su
        # --- Final balance ---------------------------------------------------
        out     = d_ci, d_hi
        self.T_ci = T_ci
        self.T_hi = T_hi
        self.T_wh = T_wh
        return out

    def resi_mass_flow(self,iter_var):
        """
        Residuals of different parameters to control the heat engine.
        """
        # --- Evaluate the massflow -------------------------------------------
        cheat_m_he     = iter_var
        # --- Check if match the demanded power input -------------------------
        if  self.parameters['mode'] == 'power':
            w_in_he_en   = self.i_he_2   - self.i_he_1
            W_in_he_en   = cheat_m_he    * w_in_he_en
            w_out_he_en  = self.i_he_3   - self.i_he_4
            W_out_he_en  = cheat_m_he    * w_out_he_en
            cheat_P_he   = W_out_he_en   - W_in_he_en
            resi         = self.P_he     - cheat_P_he
        # --- Check if match the maximum heat output --------------------------
        if  self.parameters['mode'] == 'sink':
            q_available  = self.i_he_cs_ex  - self.i_he_cs_su
            Q_available  = self.m_he_cs     * q_available
            q_required   = self.i_he_4r     - self.i_he_1
            m_required   = Q_available/q_required
            resi         = cheat_m_he       - m_required
        # --- Check if match the maximum heat input ---------------------------
        if  self.parameters['mode'] == 'source':
            q_available  = self.i_he_hs_su  - self.i_he_hs_ex
            Q_available  = self.m_he_hs     * q_available
            q_required   = self.i_he_3      - self.i_he_2r
            m_required   = Q_available/q_required
            resi         = cheat_m_he       - m_required
        return resi
    
    def find_flow_secondary_hs(self):
        """
        Finds the mass flow rate of the hot secondary fluid to match the 
        specified temperature profile based on the refrigerant flow rate.
        """
        q_available  = self.i_he_3      - self.i_he_2r
        Q_available  = self.m_he        * q_available
        q_required   = self.i_he_hs_su  - self.i_he_hs_ex
        if q_required > 0:  m_required = Q_available/q_required
        else:               m_required = 0
        return m_required
    
    def find_flow_secondary_cs(self):
        """
        Finds the mass flow rate of the cold secondary fluid to match the 
        specified temperature profile based on the refrigerant flow rate.
        """
        q_available  = self.i_he_4r    - self.i_he_1
        Q_available  = self.m_he       * q_available
        q_required   = self.i_he_cs_ex - self.i_he_cs_su
        if q_required > 0:  m_required = Q_available/q_required
        else:               m_required = 0
        return m_required

    def check_consistency(self):
        """
        Check that the results are consistent.
        """
        if self.s_he_4x > self.s_he_4:                  self.error = True
        if self.s_he_4r > self.s_he_4:                  self.error = True
        if self.s_he_4x > self.s_he_4r:                 self.error = True
        if self.s_he_1x < self.s_he_1:                  self.error = True
        if self.s_he_2x < self.s_he_2:                  self.error = True
        if self.s_he_3x > self.s_he_3:                  self.error = True
        if self.s_he_4  < self.s_he_3:                  self.error = True
        if self.s_he_3x < self.s_he_2x:                 self.error = True
        if self.s_he_4x < self.s_he_1x:                 self.error = True
        if self.T_he_2r > self.T_he_hs_ex-0.95*self.parameters['dT_he_ev_pp']:
            self.error = True
        if self.T_he_3  > self.T_he_hs_su-0.95*self.parameters['dT_he_ev_pp']:
            self.error = True
        if self.T_he_4r < self.T_he_cs_ex+0.95*self.parameters['dT_he_cd_pp']: 
            self.error = True
        if self.T_he_1  < self.T_he_cs_su+0.95*self.parameters['dT_he_cd_pp']:
            self.error = True
        if self.eta_he_cyclen < 0:                      self.error = True
        if self.T_hi-self.T_wh    < 0.95*self.parameters['dT_he_ev_pp']:  
            self.error = True
        if self.T_he_4x-self.T_ci < 0.95*self.parameters['dT_he_cd_pp']:    
            self.error = True 
        
        if self.p_he_1x < 0.99*self.p_he_1x_min: self.error = True
        if self.p_he_3x > 1.00*self.p_he_3x_max: self.error = True

        if abs(self.resi_1) > 1e-2: self.error = True
        if abs(self.resi_2) > 1e-2: self.error = True
        
        if not bool(self.parameters['wet_ex']):
            p_vec = np.linspace(self.p_he_3,self.p_he_4)
            for i, p in enumerate(p_vec):
                self.state_he.update(CoolProp.PSmass_INPUTS,p,self.s_he_3) 
                h_is = self.state_he.hmass()
                h    = self.i_he_3-self.eta_is_ex*(self.i_he_3-h_is)
                self.state_he.update(CoolProp.HmassP_INPUTS,h,p)
                if 0 < self.state_he.Q() < 1:
                    self.error = True
                if  self.error: break
            
    def retrieve_kpi(self):
        """
        Retrieve the KPI's of the SRORC
        """
        
        if self.options['exergy']: self.compute_exergy()
        
        # --- Basic quantities ------------------------------------------------
        self.w_in_he_en = self.i_he_2 -self.i_he_1 
        self.q_in_he_en = self.i_he_3 -self.i_he_2r
        self.w_out_he_en= self.i_he_3 -self.i_he_4 
        self.q_out_he_en= self.i_he_4r-self.i_he_1
        self.w_net_he_en= self.w_out_he_en-self.w_in_he_en
        self.f_he_hs    = self.q_in_he_en /(self.i_he_hs_su-self.i_he_hs_ex)
        self.f_he_cs    = self.q_out_he_en/(self.i_he_cs_ex-self.i_he_cs_su)
        
        if self.options['exergy']:
            self.w_in_he_ex = self.e_he_2 -self.e_he_1 
            self.q_in_he_ex = self.e_he_3 -self.e_he_2r
            self.w_out_he_ex= self.e_he_3 -self.e_he_4 
            self.q_out_he_ex= self.e_he_4r-self.e_he_1 
            self.w_net_he_ex= self.w_out_he_ex-self.w_in_he_ex
        
        # --- Efficiencies ----------------------------------------------------
        self.eta_he_cyclen = self.w_net_he_en/self.q_in_he_en
        
        if self.options['exergy']:
            self.eta_he_cyclex = self.w_net_he_en/self.q_in_he_ex
            self.eta_he_pumpex = self.w_in_he_ex /self.w_in_he_en
            self.eta_he_expex  = self.w_out_he_en/self.w_out_he_ex
            self.eta_he_regex  = (self.e_he_2r-self.e_he_2) \
                               / (self.e_he_4-self.e_he_4r)
            self.eta_he_evapen = self.f_he_hs**(-1) * self.q_in_he_en\
                               / (self.i_he_hs_su-self.i_he_hs_ex)
            self.eta_he_evapex = self.f_he_hs**(-1) * self.q_in_he_ex\
                               / (self.e_he_hs_su-self.e_he_hs_ex)
            self.eta_he_srcen  = (self.i_he_hs_su-self.i_he_hs_ex) \
                               / (self.i_he_hs_su-self.i_he_hs_ref)
            self.eta_he_srcex  = (self.e_he_hs_su-self.e_he_hs_ex) \
                               / (self.e_he_hs_su-self.e_he_hs_ref)
            self.eta_he_genven = self.eta_he_evapen*self.eta_he_srcen
            self.eta_he_genvex = self.eta_he_evapex*self.eta_he_srcex  
            self.eta_he_condex = self.f_he_cs*(self.e_he_cs_ex\
                               - self.e_he_cs_su) / self.q_out_he_ex
            self.eta_he_toten  = self.eta_he_cyclen*self.eta_he_genven
            self.eta_he_totex  = self.eta_he_cyclex*self.eta_he_genvex
        
        # --- Losses ----------------------------------------------------------
        if self.options['exergy']:
            self.loss_he_pumpex = self.w_in_he_en-self.w_in_he_ex
            self.loss_he_expex  = self.w_out_he_ex-self.w_out_he_en
            self.loss_he_regex  = (self.e_he_4-self.e_he_4r) \
                                - (self.e_he_2r-self.e_he_2)
            self.loss_he_evapex = self.f_he_hs*(self.e_he_hs_su\
                                - self.e_he_hs_ex) - self.q_in_he_ex
            self.loss_he_srcen  = self.f_he_hs*(self.i_he_hs_ex\
                                - self.i_he_hs_ref)
            self.loss_he_srcex  = self.f_he_hs*(self.e_he_hs_ex\
                                - self.e_he_hs_ref)
            self.loss_he_condex = self.q_out_he_ex-self.f_he_cs \
                                * (self.e_he_cs_ex-self.e_he_cs_su)
            self.loss_he_sinken = self.f_he_cs*(self.i_he_cs_ex\
                                                -self.i_he_cs_su)
            self.loss_he_sinkex = self.f_he_cs*(self.e_he_cs_ex\
                                                -self.e_he_cs_su)
            self.pow_he_supplen = self.f_he_hs*(self.i_he_hs_su\
                                                -self.i_he_hs_ref)
            self.pow_he_supplex = self.f_he_hs*(self.e_he_hs_su\
                                                -self.e_he_hs_ref)
        
        # --- Carnot and Lorentz efficiencies ---------------------------------
        T_cs_logAv =(self.T_he_cs_ex-self.T_he_cs_su)\
             /np.log(self.T_he_cs_ex/self.T_he_cs_su)
        T_hs_logAv =(self.T_he_hs_su-self.T_he_hs_ex)\
             /np.log(self.T_he_hs_su/self.T_he_hs_ex)
        ETA_lorenz =(     T_hs_logAv-     T_cs_logAv)/     T_hs_logAv
        ETA_carnot =(self.T_he_hs_su-self.T_he_cs_su)/self.T_he_hs_su
        self.eta_he_lorenz = ETA_lorenz
        self.eta_he_carnot = ETA_carnot
        self.psi_he_lorenz = self.eta_he_cyclen/ETA_lorenz
        self.psi_he_carnot = self.eta_he_cyclen/ETA_carnot
        
        # --- Preliminary design information ----------------------------------
        self.vol_coef       = self.v_he_3/self.w_out_he_en
        self.exp_ratio      = self.p_he_3/self.p_he_4
        self.vol_ratio      = self.v_he_4/self.v_he_3
        
        if self.parameters['version'] == 'operational_light':
            # --- Multiply the specific quantities by the mass flow -----------
            self.Q_in_he_en = self.m_he*self.q_in_he_en
            self.W_in_he_en = self.m_he*self.w_in_he_en
            self.Q_out_he_en= self.m_he*self.q_out_he_en
            self.W_out_he_en= self.m_he*self.w_out_he_en
            self.W_net_he_en= self.m_he*self.w_net_he_en 
       
    def export_states(self):
        """
        Export the thermodynamic states of the SRORC necessary to draw
        the T-s and p-h diagrams.
        """
        self.p_he = self.p_he_1, self.p_he_1x,self.p_he_2, self.p_he_2x,\
                    self.p_he_2r,self.p_he_3, self.p_he_3x,self.p_he_4,\
                    self.p_he_4x,self.p_he_4r
        self.T_he = self.T_he_1, self.T_he_1x,self.T_he_2, self.T_he_2x,\
                    self.T_he_2r,self.T_he_3, self.T_he_3x,self.T_he_4,\
                    self.T_he_4x,self.T_he_4r
        self.i_he = self.i_he_1, self.i_he_1x,self.i_he_2, self.i_he_2x,\
                    self.i_he_2r,self.i_he_3, self.i_he_3x,self.i_he_4,\
                    self.i_he_4x,self.i_he_4r
        self.s_he = self.s_he_1, self.s_he_1x,self.s_he_2, self.s_he_2x,\
                    self.s_he_2r,self.s_he_3, self.s_he_3x,self.s_he_4,\
                    self.s_he_4x,self.s_he_4r
        self.x_he = self.x_he_1, self.x_he_1x,self.x_he_2, self.x_he_2x,\
                    self.x_he_2r,self.x_he_3, self.x_he_3x,self.x_he_4,\
                    self.x_he_4x,self.x_he_4r
        if self.options['exergy']: 
            self.e_he = self.e_he_1, self.e_he_1x,self.e_he_2, self.e_he_2x,\
                        self.e_he_2r,self.e_he_3, self.e_he_3x,self.e_he_4,\
                        self.e_he_4x,self.e_he_4r
        else:
            self.e_he = 0,0,0,0,0,0,0,0,0,0
        out = self.p_he,self.T_he,self.i_he,self.s_he,self.x_he,self.e_he
        return out    
        
    def compute_exergy(self):
        """
        Compute the exergy at each state of the SRORC.
        """
        self.compute_ref_states()
        
        # 1 - Exergy of the cycle ---------------------------------------------
        
        self.e_he_1  = (self.i_he_1 -self.i_he_0)\
          -self.T_he_0*(self.s_he_1 -self.s_he_0)
        self.e_he_1x = (self.i_he_1x-self.i_he_0)\
          -self.T_he_0*(self.s_he_1x-self.s_he_0)
        self.e_he_2  = (self.i_he_2 -self.i_he_0)\
          -self.T_he_0*(self.s_he_2 -self.s_he_0)
        self.e_he_2x = (self.i_he_2x-self.i_he_0)\
          -self.T_he_0*(self.s_he_2x-self.s_he_0)
        self.e_he_2r = (self.i_he_2r-self.i_he_0)\
          -self.T_he_0*(self.s_he_2r-self.s_he_0)
        self.e_he_3  = (self.i_he_3 -self.i_he_0)\
          -self.T_he_0*(self.s_he_3 -self.s_he_0)
        self.e_he_3x = (self.i_he_3x-self.i_he_0)\
          -self.T_he_0*(self.s_he_3x-self.s_he_0)
        self.e_he_4  = (self.i_he_4 -self.i_he_0)\
          -self.T_he_0*(self.s_he_4 -self.s_he_0)
        self.e_he_4x = (self.i_he_4x-self.i_he_0)\
          -self.T_he_0*(self.s_he_4x-self.s_he_0)
        self.e_he_4r = (self.i_he_4r-self.i_he_0)\
          -self.T_he_0*(self.s_he_4r-self.s_he_0)
        
        if round(self.T_he_1 ,2)==round(self.T_he_0,2): self.e_he_1  = 0
        if round(self.T_he_1x,2)==round(self.T_he_0,2): self.e_he_1x = 0
        if round(self.T_he_2 ,2)==round(self.T_he_0,2): self.e_he_2  = 0
        if round(self.T_he_2x,2)==round(self.T_he_0,2): self.e_he_2x = 0
        if round(self.T_he_2r,2)==round(self.T_he_0,2): self.e_he_2r = 0
        if round(self.T_he_3 ,2)==round(self.T_he_0,2): self.e_he_3  = 0
        if round(self.T_he_3x,2)==round(self.T_he_0,2): self.e_he_3x = 0
        if round(self.T_he_4 ,2)==round(self.T_he_0,2): self.e_he_4  = 0
        if round(self.T_he_4x,2)==round(self.T_he_0,2): self.e_he_4x = 0
        if round(self.T_he_4r,2)==round(self.T_he_0,2): self.e_he_4r = 0
        
        # 2 - Exergy of the boundaries ----------------------------------------
        
        self.compute_exergy_boundaries()

#%% MODEL: TBORC

class TBORC(SBORC):
    """
    Class for the simulation of Transcritical Basic Organic Rankine Cycles.

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (heat engine boundaries and operational 
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
    # ==== HEAT ENGINE MODEL ================================================ #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the hot source
            * 'sink':   available heat flow rate at the cold sink
            * 'power':  electrical power produced by the heat engine
        
        """        

        if self.T_he_hs_su-self.parameters['dT_he_ev_pp'] < \
           self.state_he.T_critical():
            raise ValueError('The heat source is below critical temperature!\
                             In: '+os.path.join(os.path.abspath(__file__),
                             'TBORC','evaluate_cycle'))
                
        # --- Cycle computation -----------------------------------------------
        if self.parameters['m_rat'] == 0: self.m_rat = self.opt_eta()
        else:                             self.m_rat = self.parameters['m_rat']
        
        # ..... saturated (liquid) ............................................
        self.p_he_1x, self.p_he_3 = self.find_p()
        self.x_he_1x = 0.0
        self.state_he.update(CoolProp.PQ_INPUTS,self.p_he_1x,self.x_he_1x)
        self.T_he_1x = self.state_he.T()
        self.i_he_1x = self.state_he.hmass()
        self.s_he_1x = self.state_he.smass()
        # ..... liquid ........................................................
        self.p_he_1  = self.p_he_1x
        self.T_he_1  = self.T_he_1x-self.parameters['dT_he_cd_sc']
        self.state_he.specify_phase(CoolProp.iphase_liquid)
        self.state_he.update(CoolProp.PT_INPUTS,self.p_he_1,self.T_he_1)
        self.i_he_1  = self.state_he.hmass()     
        self.x_he_1  = self.state_he.Q()
        self.s_he_1  = self.state_he.smass()
        self.state_he.unspecify_phase()
        # ..... saturated (gas) ...............................................
        self.p_he_4x = self.p_he_1x/(1-self.parameters['dp_he_cd'])
        self.x_he_4x = 1.0
        self.state_he.update(CoolProp.PQ_INPUTS,self.p_he_4x,self.x_he_4x)
        self.T_he_4x = self.state_he.T()
        self.i_he_4x = self.state_he.hmass()
        self.s_he_4x = self.state_he.smass()
        # ..... liquid (super-critical) .......................................
        self.p_he_2  = self.p_he_3/(1-self.parameters['dp_he_ev'])
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,self.p_he_2,
                                                        self.s_he_1)
            self.i_he_2s = self.state_he.hmass()
            self.i_he_2  = self.i_he_1\
                         +(self.i_he_2s-self.i_he_1)/self.parameters['eta_pm']
            self.state_he.update(CoolProp.HmassP_INPUTS,
                                 self.i_he_2,self.p_he_2)
            self.T_he_2  = self.state_he.T()
            self.P_pm    = 0
            self.eta_is_pm = self.parameters['eta_pm']
        self.x_he_2  = self.state_he.Q()
        self.s_he_2  = self.state_he.smass()
        # ..... gas (super-critical) ..........................................
        self.i_he_3  = self.m_rat*(self.i_he_hs_su-self.i_he_hs_ex)+self.i_he_2
        self.state_he.specify_phase(CoolProp.iphase_supercritical)
        self.state_he.update(CoolProp.HmassP_INPUTS,self.i_he_3,self.p_he_3)
        self.T_he_3  = self.state_he.T()   
        self.x_he_3  = self.state_he.Q()
        self.s_he_3  = self.state_he.smass()
        self.v_he_3  = self.state_he.rhomass()**(-1)
        self.state_he.unspecify_phase()
        # ..... gas ...........................................................
        self.p_he_4  = self.p_he_4x
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,
                                 self.p_he_4,self.s_he_3)
            self.i_he_4s = self.state_he.hmass()
            self.i_he_4  = self.i_he_3\
                         -(self.i_he_3-self.i_he_4s)\
                         * self.parameters['eta_max_ex']
            self.m_he    = 0
            self.P_ex    = 0
            self.eta_is_ex = self.parameters['eta_max_ex']
            self.ex_flag = 'na'
            self.state_he.update(CoolProp.HmassP_INPUTS,
                                 self.i_he_4,self.p_he_4)
            self.T_he_4  = self.state_he.T()
        self.s_he_4  = self.state_he.smass()
        self.x_he_4  = self.state_he.Q()
        self.v_he_4  = self.state_he.rhomass()**(-1)
        
        # --- Evaluate the operational conditions -----------------------------
        
        if self.parameters['version'] == 'operational_light':
            self.m_he    = self.find_mass_flow()
            self.P_ex    = self.m_he    * (self.i_he_3 - self.i_he_4)
            self.m_he_cs = self.find_flow_secondary_cs()
            self.m_he_hs = self.find_flow_secondary_hs()
        
        # --- Evaluate the Key Performance Indicators -------------------------
        
        self.retrieve_kpi()
        
    # ======================================================================= #
    # ==== SUB-FUNCTIONS FOR HEAT ENGINE MODEL ============================== #
    # ======================================================================= #

    def find_p(self):
        """
        Find the condensing and evaporating pressures in the TBORC based
        on the pinch-point and subcooling, and for the selected m_rat.
        Note that phase change is not permitted for the secondary fluid.
        """
        # --- Minimum condenser pressure (firm) -------------------------------
        T_he_1x_min         = self.T_he_cs_su \
                            + self.parameters['dT_he_cd_pp'] \
                            + self.parameters['dT_he_cd_sc']
        if T_he_1x_min > self.state_he.T_critical():
            raise ValueError('An inconsistency was detected in the cycle!\
                              In: '+os.path.join(os.path.abspath(__file__),
                              'TBORC','find_p:\
            the min. condensation temperature is above the critical point!'))
        self.state_he.update(CoolProp.QT_INPUTS,0.0,T_he_1x_min)  
        self.p_he_1x_min    = self.state_he.p()
        # --- Maximum condenser pressure (indicative) -------------------------
        T_he_1x_max         = max(self.T_he_cs_ex\
                                 +self.parameters['dT_he_cd_pp'],
                                  T_he_1x_min)
        if T_he_1x_max > self.state_he.T_critical():
            T_he_1x_max = self.state_he.T_critical() - 1
        self.state_he.update(CoolProp.QT_INPUTS,0.0,T_he_1x_max)
        self.p_he_1x_max    = self.state_he.p() 
        # --- Minimum 'evaporator' pressure (firm) ----------------------------
        self.p_he_3_min     = self.state_he.p_critical()+1
        # --- Maximum 'evaporator' pressure (indicative) ----------------------
        self.p_he_3_max     = self.state_he.pmax()-1
        # --- Find the pressure levels ----------------------------------------
        guess               = self.p_he_1x_min, self.p_he_3_min
        try:
            res   = scipy.optimize.fsolve(self.resi_p, guess, full_output=True)
            p_sol = res[0]
        except:
            guess = self.p_he_1x_min, self.p_he_3_min 
            res   = scipy.optimize.least_squares(self.resi_p, guess, 
                     bounds = ((self.p_he_1x_min,self.p_he_1x_min),
                               (self.p_he_3_max ,self.p_he_3_max )),
                  ftol=1e-12,xtol=1e-12)
            p_sol = res.x
            import warnings
            warnings.warn('In: TBORC.find_p, fsolve failed.'+
                          ' least_squares was tested instead. Message: '+
                          str(res.message))
        self.resi_1, self.resi_2 = self.resi_p((p_sol[0],p_sol[1]))
        if abs(self.resi_1) > 1 or abs(self.resi_2) > 1:
            guess = (self.p_he_1x_min+self.p_he_1x_max)/2,\
                    (self.p_he_3_min +self.p_he_3_max) /2
            res = scipy.optimize.least_squares(self.resi_p, guess, 
                  bounds = ((self.p_he_1x_min,self.p_he_1x_min),
                            (self.p_he_3_max, self.p_he_3_max )),
                  ftol=1e-12,xtol=1e-12)
            p_sol = res.x
            import warnings
            warnings.warn('In: TBORC.find_p, fsolve did not converge.'+
                          ' least_squares was tested instead. Message: '+
                          str(res.message))
        self.resi_1, self.resi_2 = self.resi_p((p_sol[0],p_sol[1]))
        p_he_1x  = p_sol[0]
        p_he_3   = p_sol[1]
        return p_he_1x, p_he_3        

    def resi_p(self,iter_var):
        """
        Residuals of pinch-point method in condenser and 'evaporator'.
        Used to find the condensing and 'evaporating' pressures in the TBORC.
        
        Cells discretization: 0.25K per cell (for secondary fluid)
        """
        # --- Guess values ----------------------------------------------------
        p_he_1x_guess, p_he_3_guess = iter_var
        # --- Condenser -------------------------------------------------------
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_1x_guess,0.0)
        T_he_1x = self.state_he.T()
        p_he_1  = p_he_1x_guess
        T_he_1  = T_he_1x-self.parameters['dT_he_cd_sc']
        self.state_he.specify_phase(CoolProp.iphase_liquid)
        self.state_he.update(CoolProp.PT_INPUTS,p_he_1,T_he_1)
        h_he_1  = self.state_he.hmass()
        s_he_1  = self.state_he.smass()
        self.state_he.unspecify_phase()
        p_he_4x = p_he_1x_guess/(1-self.parameters['dp_he_cd'])
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_4x,1.0)
        T_he_4x = self.state_he.T()
        h_he_4x = self.state_he.hmass()
        # --- Pump ------------------------------------------------------------
        p_he_2  = p_he_3_guess/(1-self.parameters['dp_he_ev'])
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,p_he_2,s_he_1)
            h_he_2s = self.state_he.hmass()
            h_he_2  = h_he_1+(h_he_2s-h_he_1)/self.parameters['eta_pm']
        # --- 'Evaporator' ----------------------------------------------------
        cell_resolution   = 0.05 # K
        n_elem = int(1/cell_resolution*(self.T_he_hs_su-self.T_he_hs_ex)+1)
        T_he_hs_vec       = np.linspace(self.T_he_hs_su,self.T_he_hs_ex,n_elem)
        p_he_hs_vec       = np.linspace(self.p_he_hs_su,self.p_he_hs_ex,n_elem)
        h_he_hs_vec  = np.zeros(n_elem)
        for i, T in enumerate(h_he_hs_vec):
            self.state_hs.update(CoolProp.PT_INPUTS,p_he_hs_vec[i],
                                                    T_he_hs_vec[i])
            h_he_hs_vec[i] = self.state_hs.hmass()
        T_he_hs_vec       = np.flip(T_he_hs_vec)
        h_he_hs_vec       = np.flip(h_he_hs_vec)
        h_he_wf_vec       = np.zeros(n_elem)
        h_he_wf_vec[0]    = h_he_2
        for i in range(len(h_he_hs_vec)-1):
            h_he_wf_vec[i+1] = self.m_rat*(h_he_hs_vec[i+1]-h_he_hs_vec[i])\
                             + h_he_wf_vec[i]
        p_he_wf_vec  = np.linspace(p_he_2         ,p_he_3_guess   ,n_elem)
        T_he_wf_vec  = np.zeros(n_elem)
        for i, T in enumerate(T_he_wf_vec):
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_wf_vec[i],
                                                        p_he_wf_vec[i])
            T_he_wf_vec[i] = self.state_he.T()
        p_he_3       = p_he_wf_vec[-1]
        T_he_3       = T_he_wf_vec[-1]
        self.state_he.specify_phase(CoolProp.iphase_supercritical)
        self.state_he.update(CoolProp.PT_INPUTS,p_he_3,T_he_3)
        h_he_3  = self.state_he.hmass()
        s_he_3  = self.state_he.smass()
        self.state_he.unspecify_phase()
        # --- Expander --------------------------------------------------------
        p_he_4  = p_he_4x
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,p_he_4,s_he_3)
            h_he_4s = self.state_he.hmass()
            h_he_4  = h_he_3-(h_he_3-h_he_4s)*self.parameters['eta_max_ex']
        # --- Energy balance on evaporator ------------------------------------
        imin    = np.argmin(T_he_hs_vec-T_he_wf_vec)
        T_hi    = T_he_hs_vec[imin]
        T_wi    = T_he_wf_vec[imin]
        d_hi    = T_he_hs_vec[imin]-T_he_wf_vec[imin]\
                                   -self.parameters['dT_he_ev_pp']
        # --- Energy balance on condenser -------------------------------------
        frac_cs = (h_he_4-h_he_1)/(self.i_he_cs_ex-self.i_he_cs_su)
        T_ci    = T_he_4x-self.parameters['dT_he_cd_pp']
        self.state_cs.update(CoolProp.PT_INPUTS,self.p_he_cs_ex,T_ci)
        h_ci    = self.state_cs.hmass()  
        d_ci    = frac_cs*(h_ci-self.i_he_cs_su)-(h_he_4x-h_he_1)
        if p_he_1x_guess < self.p_he_1x_min: # pinch at the sink inlet
            T_ci    = T_he_1-self.parameters['dT_he_cd_pp']
            self.state_cs.update(CoolProp.PT_INPUTS,self.p_he_cs_su,T_ci)
            h_ci    = self.state_cs.hmass()  
            d_ci    = h_ci-self.i_he_cs_su
        # --- Final balance ---------------------------------------------------
        out     = d_hi, d_ci
        self.T_ci = T_ci
        self.T_hi = T_hi
        self.T_wi = T_wi
        return out
    
    def opt_eta(self):
        """
        Find the optimum m_rat in the TBORC based on the pinch-point and for 
        the selected subcooling.
        """
        T_he_3_min = self.state_he.T_critical()+0.1
        T_he_3_max = self.T_he_hs_su-self.parameters['dT_he_ev_pp']
        p_he_2_min = self.state_he.p_critical()+1
        p_he_3_min = p_he_2_min*(1-self.parameters['dp_he_ev'])
        
        # --- Minimum 'evaporator' pressure (firm) ----------------------------
        # ... condenser .......................................................
        T_he_1x_min= self.T_he_cs_su+self.parameters['dT_he_cd_pp']\
                                    +self.parameters['dT_he_cd_sc']
        self.state_he.update(CoolProp.QT_INPUTS,0.0,T_he_1x_min)  
        p_he_1x_min= self.state_he.p()
        p_he_1_min = p_he_1x_min
        T_he_1_min = T_he_1x_min-self.parameters['dT_he_cd_sc']
        self.state_he.update(CoolProp.PT_INPUTS,p_he_1_min,T_he_1_min)  
        h_he_1_min = self.state_he.hmass()
        s_he_1_min = self.state_he.smass()
        self.state_he.update(CoolProp.PSmass_INPUTS,p_he_2_min,s_he_1_min)  
        h_he_2s_min= self.state_he.hmass()
        
        T_he_4x_max= self.T_he_cs_ex+self.parameters['dT_he_cd_pp']
        self.state_he.update(CoolProp.QT_INPUTS,1.0,T_he_4x_max)  
        p_he_4x_max= self.state_he.p()
        p_he_1x_max= p_he_4x_max*(1-self.parameters['dp_he_cd'])
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_1x_max,0.0)  
        T_he_1x_max= self.state_he.T()
        p_he_1_max = p_he_1x_max
        T_he_1_max = T_he_1x_max-self.parameters['dT_he_cd_sc']
        self.state_he.update(CoolProp.PT_INPUTS,p_he_1_max,T_he_1_max)  
        h_he_1_max = self.state_he.hmass()
        s_he_1_max = self.state_he.smass()
        self.state_he.update(CoolProp.PSmass_INPUTS,p_he_2_min,s_he_1_max)  
        h_he_2s_max= self.state_he.hmass()
        # ... pump ............................................................
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            h_he_2_minp2_minT1 = h_he_1_min+(h_he_2s_min-h_he_1_min)\
                                                    /self.parameters['eta_pm']
            h_he_2_minp2_maxT1 = h_he_1_max+(h_he_2s_max-h_he_1_max)\
                                                    /self.parameters['eta_pm']
        # ... 'evaporator' ....................................................
        self.state_he.update(CoolProp.PT_INPUTS,p_he_3_min,T_he_3_min)  
        h_he_3_minp2_minT3 = self.state_he.hmass()
        self.state_he.update(CoolProp.PT_INPUTS,p_he_3_min,T_he_3_max)  
        h_he_3_minp2_maxT3 = self.state_he.hmass()
        den = (self.i_he_hs_su-self.i_he_hs_ex)
        r_he_minp2_minT1_minT3 = (h_he_3_minp2_minT3-h_he_2_minp2_minT1)/den
        r_he_minp2_minT1_maxT3 = (h_he_3_minp2_maxT3-h_he_2_minp2_minT1)/den
        r_he_minp2_maxT1_minT3 = (h_he_3_minp2_minT3-h_he_2_minp2_maxT1)/den
        r_he_minp2_maxT1_maxT3 = (h_he_3_minp2_maxT3-h_he_2_minp2_maxT1)/den
        # ... find the bounds and guess value for the optimisation problem ....
        r_he_min = min([r_he_minp2_minT1_minT3,r_he_minp2_minT1_maxT3,
                        r_he_minp2_maxT1_minT3,r_he_minp2_maxT1_maxT3])
        r_he_max = max([r_he_minp2_minT1_minT3,r_he_minp2_minT1_maxT3,
                        r_he_minp2_maxT1_minT3,r_he_minp2_maxT1_maxT3])
        r_he_gue = np.mean([r_he_min,r_he_max,r_he_max])#r_he_minp2_maxT1_maxT3
        # --- Maximum 'evaporator' pressure (indicative) ----------------------
        # ... skipped because we know that it is far from the actual value ...
        
        # --- Solve the optimisation problem ----------------------------------
        try:
            sol = scipy.optimize.minimize(self.resi_eta,x0=r_he_gue,
                                          bounds=[(max([0,   r_he_min]),
                                                   min([1e+5,r_he_max]))])
            res = sol.x[0]
            if res == r_he_gue:
                import warnings
                warnings.warn('In: TBORC.opt_eta, sol = guess. '+
                              'Full computation was done instead.')
                r_he_gue = np.arange(r_he_min,r_he_max,0.025)
                out      = np.zeros(len(r_he_gue))
                for i, r in enumerate(r_he_gue):
                    try:    out[i] = self.resi_eta([r])
                    except: out[i] = 1e+5
                imax = np.argmax(1/out)
                res = r_he_gue[imax]
                # refine the results
                r_he_gue = np.arange(0.96*res,1.04*res,0.001)
                out      = np.zeros(len(r_he_gue))
                for i, r in enumerate(r_he_gue):
                    try:    out[i] = self.resi_eta([r])
                    except: out[i] = 1e+5
                imax = np.argmax(1/out)
                res = r_he_gue[imax]
        except:
            import warnings
            warnings.warn('In: TBORC.opt_eta, minimize failed. '+
                          'Full computation was done instead.')
            r_he_gue = np.arange(r_he_min,r_he_max,0.025)
            out      = np.zeros(len(r_he_gue))
            for i, r in enumerate(r_he_gue):
                try:    out[i] = self.resi_eta([r])
                except: out[i] = 1e+5
            imax = np.argmax(1/out)
            res = r_he_gue[imax]
            # refine the results
            r_he_gue = np.arange(0.96*res,1.04*res,0.001)
            out      = np.zeros(len(r_he_gue))
            for i, r in enumerate(r_he_gue):
                try:    out[i] = self.resi_eta([r])
                except: out[i] = 1e+5
            imax = np.argmax(1/out)
            res = r_he_gue[imax]
        return res
    
    def resi_eta(self,iter_var):
        """
        Inverse of eta_cyclen (function to minimise).  
        Used to find the optimum m_rat in the TBORC.
        """
        # --- Guess values ----------------------------------------------------
        self.m_rat = iter_var[0]
        p_he_1x, p_he_3 = self.find_p()
        # --- Condenser -------------------------------------------------------
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_1x,0.0)
        T_he_1x = self.state_he.T()
        p_he_1  = p_he_1x
        T_he_1  = T_he_1x-self.parameters['dT_he_cd_sc']
        self.state_he.specify_phase(CoolProp.iphase_liquid)
        self.state_he.update(CoolProp.PT_INPUTS,p_he_1,T_he_1)
        h_he_1  = self.state_he.hmass()
        s_he_1  = self.state_he.smass()
        self.state_he.unspecify_phase()
        p_he_4x = p_he_1x/(1-self.parameters['dp_he_cd'])
        # --- Pump ------------------------------------------------------------
        p_he_2  = p_he_3/(1-self.parameters['dp_he_ev'])
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,p_he_2,s_he_1)
            h_he_2s = self.state_he.hmass()
            h_he_2  = h_he_1+(h_he_2s-h_he_1)/self.parameters['eta_pm']
        # --- 'Evaporator' ----------------------------------------------------
        h_he_3 = self.m_rat*(self.i_he_hs_su-self.i_he_hs_ex)+h_he_2
        self.state_he.specify_phase(CoolProp.iphase_supercritical)
        self.state_he.update(CoolProp.HmassP_INPUTS,h_he_3,p_he_3)
        T_he_3  = self.state_he.T()   
        s_he_3  = self.state_he.smass()
        self.state_he.unspecify_phase()
        # --- Expander --------------------------------------------------------
        p_he_4  = p_he_4x
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,p_he_4,s_he_3)
            h_he_4s = self.state_he.hmass()
            h_he_4  = h_he_3-(h_he_3-h_he_4s)*self.parameters['eta_max_ex']
        # --- Proceed to efficiency -------------------------------------------
        w_in_en     = h_he_2 -h_he_1 
        q_in_en     = h_he_3 -h_he_2
        w_out_en    = h_he_3 -h_he_4 
        w_net_en    = w_out_en-w_in_en
        eta_cyclen  = w_net_en/q_in_en
        out = 1/eta_cyclen
        if round(T_he_3,2)>round(self.T_he_hs_su
                                -self.parameters['dT_he_ev_pp'],2): out = 1e+5
        return out
    
    def check_consistency(self):
        """
        Check that the results are consistent.
        """
        if self.s_he_4x > self.s_he_4:                  self.error = True
        if self.s_he_1x < self.s_he_1:                  self.error = True
        if self.s_he_2  < self.s_he_1:                  self.error = True
        if self.s_he_2  > self.s_he_3:                  self.error = True
        if self.s_he_4  < self.s_he_3:                  self.error = True
        if self.s_he_4x < self.s_he_1x:                 self.error = True
        if self.T_he_2  > self.T_he_hs_ex-0.95*self.parameters['dT_he_ev_pp']: 
            self.error = True
        if self.T_he_3  > self.T_he_hs_su-0.95*self.parameters['dT_he_ev_pp']:
            self.error = True
        if self.T_he_4  < self.T_he_cs_ex+0.95*self.parameters['dT_he_cd_pp']: 
            self.error = True
        if self.T_he_1  < self.T_he_cs_su+0.95*self.parameters['dT_he_cd_pp']: 
            self.error = True
        if self.eta_he_cyclen < 0:                      self.error = True
        if self.T_hi-self.T_wi          < 0.95*self.parameters['dT_he_ev_pp']:    
            self.error = True
        if self.T_he_4x-self.T_ci       < 0.95*self.parameters['dT_he_cd_pp']:    
            self.error = True  
        
        if self.p_he_1x < 0.99*self.p_he_1x_min: self.error = True
        if self.p_he_3  < 1.00*self.p_he_3_min:  self.error = True
        
        if abs(self.resi_1) > 1: self.error = True
        if abs(self.resi_2) > 1: self.error = True
        
        if not bool(self.parameters['wet_ex']):
            p_vec = np.linspace(self.p_he_3,self.p_he_4)
            for i, p in enumerate(p_vec):
                self.state_he.update(CoolProp.PSmass_INPUTS,p,self.s_he_3) 
                h_is = self.state_he.hmass()
                h    = self.i_he_3-self.eta_is_ex*(self.i_he_3-h_is)
                self.state_he.update(CoolProp.HmassP_INPUTS,h,p)
                if 0 < self.state_he.Q() < 1:
                    # print(self.state_he.Q())
                    self.error = True
                if  self.error: break

    def export_states(self):
        """
        Export the thermodynamic states of the SBORC necessary to draw
        the T-s and p-h diagrams.
        """
        self.p_he = self.p_he_1,self.p_he_1x,self.p_he_2,\
                    self.p_he_3,self.p_he_4,self.p_he_4x
        self.T_he = self.T_he_1,self.T_he_1x,self.T_he_2,\
                    self.T_he_3,self.T_he_4,self.T_he_4x
        self.i_he = self.i_he_1,self.i_he_1x,self.i_he_2,\
                    self.i_he_3,self.i_he_4,self.i_he_4x
        self.s_he = self.s_he_1,self.s_he_1x,self.s_he_2,\
                    self.s_he_3,self.s_he_4,self.s_he_4x
        self.x_he = self.x_he_1,self.x_he_1x,self.x_he_2,\
                    self.x_he_3,self.x_he_4,self.x_he_4x
        if self.options['exergy']: 
            self.e_he = self.e_he_1,self.e_he_1x,self.e_he_2,\
                        self.e_he_3,self.e_he_4,self.e_he_4x
        else:
            self.e_he = 0,0,0,0,0,0
        out = self.p_he,self.T_he,self.i_he,self.s_he,self.x_he,self.e_he
        return out
    
    def compute_exergy(self):
        """
        Compute the exergy at each state of the TBORC.
        """
        self.compute_ref_states()
        
        # 1 - Exergy of the cycle ---------------------------------------------
        
        self.e_he_1  = (self.i_he_1 -self.i_he_0)\
          -self.T_he_0*(self.s_he_1 -self.s_he_0)
        self.e_he_1x = (self.i_he_1x-self.i_he_0)\
          -self.T_he_0*(self.s_he_1x-self.s_he_0)
        self.e_he_2  = (self.i_he_2 -self.i_he_0)\
          -self.T_he_0*(self.s_he_2 -self.s_he_0)
        self.e_he_3  = (self.i_he_3 -self.i_he_0)\
          -self.T_he_0*(self.s_he_3 -self.s_he_0)
        self.e_he_4  = (self.i_he_4 -self.i_he_0)\
          -self.T_he_0*(self.s_he_4 -self.s_he_0)
        self.e_he_4x = (self.i_he_4x-self.i_he_0)\
          -self.T_he_0*(self.s_he_4x-self.s_he_0)
        
        if round(self.T_he_1 ,2)==round(self.T_he_0,2): self.e_he_1  = 0
        if round(self.T_he_1x,2)==round(self.T_he_0,2): self.e_he_1x = 0
        if round(self.T_he_2 ,2)==round(self.T_he_0,2): self.e_he_2  = 0
        if round(self.T_he_3 ,2)==round(self.T_he_0,2): self.e_he_3  = 0
        if round(self.T_he_4 ,2)==round(self.T_he_0,2): self.e_he_4  = 0
        if round(self.T_he_4x,2)==round(self.T_he_0,2): self.e_he_4x = 0
        
        # 2 - Exergy of the boundaries ----------------------------------------
        
        self.compute_exergy_boundaries()

#%% MODEL: TRORC

class TRORC(TBORC):
    """
    Class for the simulation of Transcritical Recuperated Organic Rankine 
    Cycles.

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (heat engine boundaries and operational 
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
    # ==== HEAT ENGINE MODEL ================================================ #
    # ======================================================================= #
    
    def evaluate_cycle(self):
        """
        Evalutate the thermodynamic cycle and operations.
        
        Controlled to match (only for operational versions):
            * 'source': available heat flow rate at the hot source
            * 'sink':   available heat flow rate at the cold sink
            * 'power':  electrical power produced by the heat engine
        
        """        

        if self.T_he_hs_su-self.parameters['dT_he_ev_pp'] < \
           self.state_he.T_critical():
            raise ValueError('The heat source is below critical temperature!\
                             In: '+os.path.join(os.path.abspath(__file__),
                             'TRORC','evaluate_cycle'))
                
        # --- Cycle computation -----------------------------------------------
        if self.parameters['m_rat'] == 0: self.m_rat = self.opt_eta()
        else:                             self.m_rat = self.parameters['m_rat']

        # ..... saturated (liquid) ............................................
        self.p_he_1x, self.p_he_3 = self.find_p()
        self.x_he_1x = 0.0
        self.state_he.update(CoolProp.PQ_INPUTS,self.p_he_1x,self.x_he_1x)
        self.T_he_1x = self.state_he.T()
        self.i_he_1x = self.state_he.hmass()
        self.s_he_1x = self.state_he.smass()
        # ..... liquid ........................................................
        self.p_he_1  = self.p_he_1x
        self.T_he_1  = self.T_he_1x-self.parameters['dT_he_cd_sc']
        self.state_he.specify_phase(CoolProp.iphase_liquid)
        self.state_he.update(CoolProp.PT_INPUTS,self.p_he_1,self.T_he_1)
        self.i_he_1  = self.state_he.hmass()    
        self.x_he_1  = self.state_he.Q()
        self.s_he_1  = self.state_he.smass()
        self.state_he.unspecify_phase()
        # ..... saturated (gas) ...............................................
        self.p_he_4x = self.p_he_1x/(1-self.parameters['dp_he_cd'])
        self.x_he_4x = 1.0
        self.state_he.update(CoolProp.PQ_INPUTS,self.p_he_4x,self.x_he_4x)
        self.T_he_4x = self.state_he.T()
        self.i_he_4x = self.state_he.hmass()
        self.s_he_4x = self.state_he.smass()
        # ..... liquid (super-critical) .......................................
        self.p_he_2  = self.p_he_3/((1-self.parameters['dp_he_ev'])\
                                   *(1-self.parameters['dp_he_rg_lq']))
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,self.p_he_2,
                                                        self.s_he_1)
            self.i_he_2s = self.state_he.hmass()
            self.i_he_2  = self.i_he_1\
                         +(self.i_he_2s-self.i_he_1)/self.parameters['eta_pm']
            self.state_he.update(CoolProp.HmassP_INPUTS,
                                 self.i_he_2,self.p_he_2)
            self.T_he_2  = self.state_he.T()
            self.P_pm    = 0
            self.eta_is_pm = self.parameters['eta_pm']
        self.x_he_2  = self.state_he.Q()
        self.s_he_2  = self.state_he.smass()
        # ..... gas (super-critical) ..........................................
        self.p_he_2r = self.p_he_3/(1-self.parameters['dp_he_ev'])
        extra        = self.T_he_2,self.i_he_2,self.p_he_2r,\
                       self.p_he_3,self.p_he_4x
        self.i_he_2r = self.find_recup(extra)
        self.i_he_3  = self.m_rat*(self.i_he_hs_su-self.i_he_hs_ex)\
                     + self.i_he_2r  
        self.state_he.specify_phase(CoolProp.iphase_supercritical)
        self.state_he.update(CoolProp.HmassP_INPUTS,self.i_he_3,self.p_he_3)
        self.T_he_3  = self.state_he.T()   
        self.x_he_3  = self.state_he.Q()
        self.s_he_3  = self.state_he.smass()
        self.v_he_3  = self.state_he.rhomass()**(-1)
        self.state_he.unspecify_phase()
        # ..... gas ...........................................................
        self.p_he_4  = self.p_he_4x/(1-self.parameters['dp_he_rg_vp'])
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,
                                 self.p_he_4,self.s_he_3)
            self.i_he_4s = self.state_he.hmass()
            self.i_he_4  = self.i_he_3-(self.i_he_3-self.i_he_4s)\
                         * self.parameters['eta_max_ex']
            self.m_he    = 0
            self.P_ex    = 0
            self.eta_is_ex = self.parameters['eta_max_ex']
            self.ex_flag = 'na'
            self.state_he.update(CoolProp.HmassP_INPUTS,
                                 self.i_he_4,self.p_he_4)
            self.T_he_4  = self.state_he.T()
        self.s_he_4  = self.state_he.smass()
        self.x_he_4  = self.state_he.Q()
        self.v_he_4  = self.state_he.rhomass()**(-1)
        # ..... gas (potentially saturated) ...................................
        self.p_he_4r = self.p_he_4x
        self.T_he_4r = self.T_he_4-(self.T_he_4-self.T_he_2)\
                     * self.parameters['epsilon']
        if self.T_he_4r == self.T_he_4x:           
            self.state_he.update(CoolProp.QT_INPUTS,1.0,self.T_he_4r)
        else:
            if  self.p_he_4r < self.state_he.p_critical()\
            and self.T_he_4r < self.state_he.T_critical():
                self.state_he.specify_phase(CoolProp.iphase_gas)
            if  self.p_he_4r < self.state_he.p_critical()\
            and self.T_he_4r > self.state_he.T_critical():
                self.state_he.specify_phase(CoolProp.iphase_supercritical_gas)
            self.state_he.update(CoolProp.PT_INPUTS,self.p_he_4r,self.T_he_4r)
            self.state_he.unspecify_phase()
        self.i_he_4r = self.state_he.hmass() 
        self.s_he_4r = self.state_he.smass() 
        self.x_he_4r = self.state_he.Q()
        # ..... liquid (potentially saturated) ................................
        self.i_he_2r = self.i_he_2+(self.i_he_4-self.i_he_4r)
        self.state_he.update(CoolProp.HmassP_INPUTS,self.i_he_2r,self.p_he_2r)
        self.T_he_2r = self.state_he.T()
        self.s_he_2r = self.state_he.smass()
        self.x_he_2r = self.state_he.Q()
        
        # --- Evaluate the operational conditions -----------------------------
        
        if self.parameters['version'] == 'operational_light':    
            self.m_he    = self.find_mass_flow()
            self.P_ex    = self.m_he*(self.i_he_3-self.i_he_4)
            self.m_he_cs = self.find_flow_secondary_cs()
            self.m_he_hs = self.find_flow_secondary_hs()
        
        # --- Evaluate the Key Performance Indicators -------------------------
        
        self.retrieve_kpi()
    
    # ======================================================================= #
    # ==== SUB-FUNCTIONS FOR HEAT ENGINE MODEL ============================== #
    # ======================================================================= #  
    
    def resi_p(self,iter_var):
        """
        Residuals of pinch-point method in condenser and 'evaporator'.
        Used to find the condensing and 'evaporating' pressures in the TRORC.
        
        Cells discretization: 0.25K per cell (for secondary fluid)
        """
        # --- Guess values ----------------------------------------------------
        p_he_1x_guess, p_he_3_guess = iter_var
        # --- Condenser -------------------------------------------------------
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_1x_guess,0.0)
        T_he_1x = self.state_he.T()
        p_he_1  = p_he_1x_guess
        T_he_1  = T_he_1x-self.parameters['dT_he_cd_sc']
        self.state_he.specify_phase(CoolProp.iphase_liquid)
        self.state_he.update(CoolProp.PT_INPUTS,p_he_1,T_he_1)
        h_he_1  = self.state_he.hmass()
        s_he_1  = self.state_he.smass()
        self.state_he.unspecify_phase()
        p_he_4x = p_he_1x_guess/(1-self.parameters['dp_he_cd'])
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_4x,1.0)
        T_he_4x = self.state_he.T()
        h_he_4x = self.state_he.hmass()
        # --- Pump ------------------------------------------------------------
        p_he_2  = p_he_3_guess/((1-self.parameters['dp_he_ev'])\
                               *(1-self.parameters['dp_he_rg_lq']))
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,p_he_2,s_he_1)
            h_he_2s = self.state_he.hmass()
            h_he_2  = h_he_1+(h_he_2s-h_he_1)/self.parameters['eta_pm']
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_2,p_he_2)
            T_he_2  = self.state_he.T()
        # --- Recuperator -----------------------------------------------------
        p_he_2r = p_he_3_guess/(1-self.parameters['dp_he_ev'])
        extra   = T_he_2, h_he_2,p_he_2r,p_he_3_guess,p_he_4x
        h_he_2r = self.find_recup(extra)
        # --- 'Evaporator' ----------------------------------------------------
        cell_resolution   = 0.05 # K
        n_elem = int(1/cell_resolution*(self.T_he_hs_su-self.T_he_hs_ex)+1)
        T_he_hs_vec       = np.linspace(self.T_he_hs_su,self.T_he_hs_ex,n_elem)
        p_he_hs_vec       = np.linspace(self.p_he_hs_su,self.p_he_hs_ex,n_elem)
        h_he_hs_vec       = np.zeros(n_elem)
        for i, T in enumerate(h_he_hs_vec):
            self.state_hs.update(CoolProp.PT_INPUTS,p_he_hs_vec[i],
                                                    T_he_hs_vec[i])
            h_he_hs_vec[i]= self.state_hs.hmass()
        T_he_hs_vec       = np.flip(T_he_hs_vec)
        h_he_hs_vec       = np.flip(h_he_hs_vec)
        h_he_wf_vec       = np.zeros(n_elem)
        h_he_wf_vec[0]    = h_he_2r
        for i in range(len(h_he_hs_vec)-1):
            h_he_wf_vec[i+1] = self.m_rat*(h_he_hs_vec[i+1]-h_he_hs_vec[i])\
                             + h_he_wf_vec[i]
        p_he_wf_vec       = np.linspace(p_he_2r,p_he_3_guess,n_elem)
        T_he_wf_vec       = np.zeros(n_elem)
        for i, T in enumerate(T_he_wf_vec):
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_wf_vec[i],
                                                        p_he_wf_vec[i])
            T_he_wf_vec[i]= self.state_he.T()
        T_he_3  = T_he_wf_vec[-1]
        self.state_he.specify_phase(CoolProp.iphase_supercritical)
        self.state_he.update(CoolProp.PT_INPUTS,p_he_3_guess,T_he_3)
        h_he_3  = self.state_he.hmass()
        s_he_3  = self.state_he.smass()
        self.state_he.unspecify_phase()
        # --- Expander --------------------------------------------------------
        p_he_4  = p_he_4x/(1-self.parameters['dp_he_rg_vp'])
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,p_he_4,s_he_3)
            h_he_4s = self.state_he.hmass()
            h_he_4  = h_he_3-(h_he_3-h_he_4s)*self.parameters['eta_max_ex']
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_4,p_he_4)
            T_he_4  = self.state_he.T()
        # --- Regenerator -----------------------------------------------------
        p_he_4r = p_he_4x
        T_he_4r = T_he_4-(T_he_4-T_he_2)*self.parameters['epsilon']
        if T_he_4r == T_he_4x:           
            self.state_he.update(CoolProp.QT_INPUTS,1.0,T_he_4r)
        else:
            self.state_he.update(CoolProp.PT_INPUTS,p_he_4r,T_he_4r)
        h_he_4r = self.state_he.hmass()       
        h_he_2r = h_he_2+(h_he_4-h_he_4r)
        # --- Energy balance on evaporator ------------------------------------
        imin    = np.argmin(T_he_hs_vec-T_he_wf_vec)
        T_hi    = T_he_hs_vec[imin]
        T_wi    = T_he_wf_vec[imin]
        d_hi    = T_he_hs_vec[imin]-T_he_wf_vec[imin]\
                                   -self.parameters['dT_he_ev_pp']
        # --- Energy balance on condenser -------------------------------------
        frac_cs = (h_he_4r-h_he_1)/(self.i_he_cs_ex-self.i_he_cs_su)
        T_ci    = T_he_4x-self.parameters['dT_he_cd_pp']
        self.state_cs.update(CoolProp.PT_INPUTS,self.p_he_cs_ex,T_ci)
        h_ci    = self.state_cs.hmass()  
        d_ci    = frac_cs*(h_ci-self.i_he_cs_su)-(h_he_4x-h_he_1)
        if p_he_1x_guess < self.p_he_1x_min: # pinch at the sink inlet
            T_ci    = T_he_1-self.parameters['dT_he_cd_pp']
            self.state_cs.update(CoolProp.PT_INPUTS,self.p_he_cs_su,T_ci)
            h_ci    = self.state_cs.hmass()  
            d_ci    = h_ci-self.i_he_cs_su
        # --- Final balance ---------------------------------------------------
        out       = d_hi, d_ci
        self.T_ci = T_ci
        self.T_hi = T_hi
        self.T_wi = T_wi
        return out
    
    def find_recup(self,extra):
        """
        Find the condensing and evaporating pressures in the TBORC based
        on the pinch-point and subcooling, and for the selected m_rat.
        """
        # --- Find the pressure levels ----------------------------------------
        T_he_2,h_he_2,p_he_2r,p_he_3,p_he_4x = extra
        # --- Find the temperatures -------------------------------------------
        guess = h_he_2*1.1
        try:
            res   = scipy.optimize.fsolve       (self.resi_recup, guess, extra, 
                                          full_output=True,)
            h_sol = res[0]
        except:
            res   = scipy.optimize.least_squares(self.resi_recup, guess, extra, 
                                                bounds = (h_he_2,h_he_2*1e+3))
            h_sol = res.x
            import warnings
            warnings.warn('In: TRORC.find_recup, fsolve failed.'+
                          ' least_squares was tested instead. Message: '+
                          str(res.message))
        try:    h_he_2r = h_sol[0]
        except: h_he_2r = h_sol
        self.resi_rec = self.resi_recup(h_he_2r,*extra)
        return h_he_2r
    
    def resi_recup(self,iter_var,*extra):
        """
        Residuals of effectiveness method in recuperator.
        Used to find the exit temperatures in the TRORC.
        
        Cells discretization: 0.25K per cell (for secondary fluid)
        """
        # --- Guess values ----------------------------------------------------
        try:    h_he_2r_guess = iter_var[0]
        except: h_he_2r_guess = iter_var
        T_he_2,h_he_2,p_he_2r,p_he_3,p_he_4x = extra
        # --- Computation -----------------------------------------------------
        h_he_3 = self.m_rat*(self.i_he_hs_su-self.i_he_hs_ex)+h_he_2r_guess
        self.state_he.specify_phase(CoolProp.iphase_supercritical)
        self.state_he.update(CoolProp.HmassP_INPUTS,h_he_3,p_he_3)
        T_he_3  = self.state_he.T()   
        s_he_3  = self.state_he.smass()
        self.state_he.unspecify_phase()
        p_he_4  = p_he_4x/(1-self.parameters['dp_he_rg_vp'])
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,p_he_4,s_he_3)
            h_he_4s = self.state_he.hmass()
            h_he_4  = h_he_3-(h_he_3-h_he_4s)*self.parameters['eta_max_ex']
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_4,p_he_4)
            T_he_4  = self.state_he.T()
        T_he_4r = T_he_4-(T_he_4-T_he_2)*self.parameters['epsilon']
        self.state_he.update(CoolProp.PT_INPUTS,p_he_4x,T_he_4r)
        h_he_4r = self.state_he.hmass()
        # --- Final energy balance --------------------------------------------
        dh_h = (h_he_4       -h_he_4r)
        dh_c = (h_he_2r_guess-h_he_2 ) 
        out  = dh_h-dh_c
        return out
    
    def opt_eta(self):
        """
        Find the optimum m_rat in the TRORC based on the pinch-point and for 
        the selected subcooling.
        """
        T_he_3_min = self.state_he.T_critical()+0.1
        T_he_3_max = self.T_he_hs_su-self.parameters['dT_he_ev_pp']
        p_he_3_min = self.state_he.p_critical()+1
        p_he_2_min = p_he_3_min/((1-self.parameters['dp_he_ev'])\
                                *(1-self.parameters['dp_he_rg_lq']))
        
        # --- Minimum 'evaporator' pressure (firm) ----------------------------
        # ... condenser .......................................................
        T_he_1x_min= self.T_he_cs_su+self.parameters['dT_he_cd_pp']\
                                    +self.parameters['dT_he_cd_sc']
        self.state_he.update(CoolProp.QT_INPUTS,0.0,T_he_1x_min)  
        p_he_1x_min= self.state_he.p()
        p_he_1_min = p_he_1x_min
        T_he_1_min = T_he_1x_min-self.parameters['dT_he_cd_sc']
        self.state_he.update(CoolProp.PT_INPUTS,p_he_1_min,T_he_1_min)  
        h_he_1_min = self.state_he.hmass()
        s_he_1_min = self.state_he.smass()
        self.state_he.update(CoolProp.PSmass_INPUTS,p_he_2_min,s_he_1_min)  
        h_he_2s_min= self.state_he.hmass()
        
        T_he_4x_max= self.T_he_cs_ex+self.parameters['dT_he_cd_pp']
        self.state_he.update(CoolProp.QT_INPUTS,1.0,T_he_4x_max)  
        p_he_4x_max= self.state_he.p()
        p_he_1x_max= p_he_4x_max*(1-self.parameters['dp_he_cd'])
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_1x_max,0.0)  
        T_he_1x_max= self.state_he.T()
        p_he_1_max = p_he_1x_max
        T_he_1_max = T_he_1x_max-self.parameters['dT_he_cd_sc']
        self.state_he.update(CoolProp.PT_INPUTS,p_he_1_max,T_he_1_max)  
        h_he_1_max = self.state_he.hmass()
        s_he_1_max = self.state_he.smass()
        self.state_he.update(CoolProp.PSmass_INPUTS,p_he_2_min,s_he_1_max)  
        h_he_2s_max= self.state_he.hmass()
        # ... pump ............................................................
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            
            h_he_2_minp2_minT1 = h_he_1_min+(h_he_2s_min-h_he_1_min)\
                                                    /self.parameters['eta_pm']
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_2_minp2_minT1,
                                                        p_he_2_min)
            T_he_2_minp2_minT1 = self.state_he.T()
            
            h_he_2_minp2_maxT1 = h_he_1_max+(h_he_2s_max-h_he_1_max)\
                                                    /self.parameters['eta_pm']
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_2_minp2_maxT1,
                                                        p_he_2_min)
            T_he_2_minp2_maxT1 = self.state_he.T()
        
        h_he_2r_minp2_minT1 = h_he_2_minp2_minT1
        h_he_2r_minp2_maxT1 = h_he_2_minp2_maxT1
        
        # ... 'evaporator' ....................................................
        self.state_he.update(CoolProp.PT_INPUTS,p_he_3_min,T_he_3_min)  
        h_he_3_minp2_minT3 = self.state_he.hmass()
        self.state_he.update(CoolProp.PT_INPUTS,p_he_3_min,T_he_3_max)  
        h_he_3_minp2_maxT3 = self.state_he.hmass()
        den = (self.i_he_hs_su-self.i_he_hs_ex)
        # ... flow rate ratios ................................................ 
        r_he_minp2_minT1_minT3 = (h_he_3_minp2_minT3-h_he_2r_minp2_minT1)/den
        r_he_minp2_minT1_maxT3 = (h_he_3_minp2_maxT3-h_he_2r_minp2_minT1)/den
        r_he_minp2_maxT1_minT3 = (h_he_3_minp2_minT3-h_he_2r_minp2_maxT1)/den
        r_he_minp2_maxT1_maxT3 = (h_he_3_minp2_maxT3-h_he_2r_minp2_maxT1)/den
        
        # --- Maximum 'evaporator' pressure (indicative) ----------------------
        # ... skipped because we know that it is far from the actual value ...
        
        # ... find the bounds and guess value for the optimisation problem ....
        r_he_min = min([r_he_minp2_minT1_minT3,r_he_minp2_minT1_maxT3,
                        r_he_minp2_maxT1_maxT3])
        r_he_max = max([r_he_minp2_minT1_minT3,r_he_minp2_minT1_maxT3,
                        r_he_minp2_maxT1_maxT3])
        r_he_gue = np.mean([r_he_min,r_he_max])
        
        # --- Solve the optimisation problem ----------------------------------
        try:
            sol = scipy.optimize.minimize(self.resi_eta,x0=r_he_gue,
                                          bounds=[(max([0,   r_he_min]),
                                                   min([1e+5,r_he_max]))])
            res = sol.x[0]
            if res == r_he_gue:
                import warnings
                warnings.warn('In: TRORC.opt_eta, sol = guess. '+
                              'Full computation was done instead.')
                r_he_gue = np.arange(r_he_min,r_he_max,0.025)
                out      = np.zeros(len(r_he_gue))
                for i, r in enumerate(r_he_gue):
                    try:    out[i] = self.resi_eta([r])
                    except: out[i] = 1e+5
                imax = np.argmax(1/out)
                # refine the results
                res = r_he_gue[imax]
                r_he_gue = np.arange(0.96*res,1.04*res,0.001)
                out      = np.zeros(len(r_he_gue))
                for i, r in enumerate(r_he_gue):
                    try:    out[i] = self.resi_eta([r])
                    except: out[i] = 1e+5
                imax = np.argmax(1/out)
                res = r_he_gue[imax]
        except:
            import warnings
            warnings.warn('In: TRORC.opt_eta, minimize failed. '+
                          'Full computation was done instead.')
            r_he_gue = np.arange(r_he_min,r_he_max,0.025)
            out      = np.zeros(len(r_he_gue))
            for i, r in enumerate(r_he_gue):
                try:    out[i] = self.resi_eta([r])
                except: out[i] = 1e+5
            imax = np.argmax(1/out)
            # refine the results
            res = r_he_gue[imax]
            r_he_gue = np.arange(0.96*res,1.04*res,0.001)
            out      = np.zeros(len(r_he_gue))
            for i, r in enumerate(r_he_gue):
                try:    out[i] = self.resi_eta([r])
                except: out[i] = 1e+5
            imax = np.argmax(1/out)
            res = r_he_gue[imax]
        return res
    
    def resi_eta(self,iter_var):
        """
        Inverse of eta_cyclen (function to minimise).  
        Used to find the optimum m_rat in the TRORC.
        """
        # --- Guess values ----------------------------------------------------
        self.m_rat = iter_var[0]
        p_he_1x, p_he_3 = self.find_p()
        # --- Condenser -------------------------------------------------------
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_1x,0.0)
        T_he_1x = self.state_he.T()
        p_he_1  = p_he_1x
        T_he_1  = T_he_1x-self.parameters['dT_he_cd_sc']
        self.state_he.specify_phase(CoolProp.iphase_liquid)
        self.state_he.update(CoolProp.PT_INPUTS,p_he_1,T_he_1)
        h_he_1  = self.state_he.hmass()
        s_he_1  = self.state_he.smass()
        self.state_he.unspecify_phase()
        p_he_4x = p_he_1x/(1-self.parameters['dp_he_cd'])
        self.state_he.update(CoolProp.PQ_INPUTS,p_he_4x,1.0)
        T_he_4x = self.state_he.T()
        # --- Pump ------------------------------------------------------------
        p_he_2  = p_he_3/((1-self.parameters['dp_he_ev'])\
                         *(1-self.parameters['dp_he_rg_lq']))
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,p_he_2,s_he_1)
            h_he_2s = self.state_he.hmass()
            h_he_2  = h_he_1+(h_he_2s-h_he_1)/self.parameters['eta_pm']
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_2,p_he_2)
            T_he_2  = self.state_he.T()
        # --- Recuperator -----------------------------------------------------
        p_he_2r = p_he_3/(1-self.parameters['dp_he_ev'])
        extra   = T_he_2,h_he_2,p_he_2r,p_he_3,p_he_4x
        h_he_2r = self.find_recup(extra)
        # --- 'Evaporator' ----------------------------------------------------
        h_he_3 = self.m_rat*(self.i_he_hs_su-self.i_he_hs_ex)+h_he_2r
        self.state_he.specify_phase(CoolProp.iphase_supercritical)
        self.state_he.update(CoolProp.HmassP_INPUTS,h_he_3,p_he_3)
        T_he_3  = self.state_he.T()   
        s_he_3  = self.state_he.smass()
        self.state_he.unspecify_phase()
        # --- Expander --------------------------------------------------------
        p_he_4  = p_he_4x/(1-self.parameters['dp_he_rg_vp'])
        if self.parameters['version'] == 'thermodynamic_full'\
        or self.parameters['version'] == 'operational_light':
            self.state_he.update(CoolProp.PSmass_INPUTS,p_he_4,s_he_3)
            h_he_4s = self.state_he.hmass()
            h_he_4  = h_he_3-(h_he_3-h_he_4s)*self.parameters['eta_max_ex']
            self.state_he.update(CoolProp.HmassP_INPUTS,h_he_4,p_he_4)
            T_he_4  = self.state_he.T()   
        # --- Regenerator -----------------------------------------------------
        p_he_4r = p_he_4x
        T_he_4r = T_he_4-(T_he_4-T_he_2)*self.parameters['epsilon']
        if T_he_4r == T_he_4x:           
            self.state_he.update(CoolProp.QT_INPUTS,1.0,T_he_4r)
        else:
            self.state_he.update(CoolProp.PT_INPUTS,p_he_4r,T_he_4r)
        h_he_4r = self.state_he.hmass()       
        h_he_2r = h_he_2 + (h_he_4-h_he_4r)
        # --- Proceed to efficiency -------------------------------------------
        w_in_en     = h_he_2 -h_he_1 
        q_in_en     = h_he_3 -h_he_2r
        w_out_en    = h_he_3 -h_he_4 
        w_net_en    = w_out_en-w_in_en
        eta_cyclen  = w_net_en/q_in_en
        out = 1/eta_cyclen
        if round(T_he_3,2)>round(self.T_he_hs_su
                                -self.parameters['dT_he_ev_pp'],2): out = 1e+5
        return out

    def resi_mass_flow(self,iter_var):
        """
        Residuals of different parameters to control the heat engine.
        """
        # --- Evaluate the massflow -------------------------------------------
        cheat_m_he     = iter_var
        # --- Check if match the demanded power input -------------------------
        if  self.parameters['mode'] == 'power':
            w_in_he_en   = self.i_he_2   - self.i_he_1
            W_in_he_en   = cheat_m_he    * w_in_he_en
            w_out_he_en  = self.i_he_3   - self.i_he_4
            W_out_he_en  = cheat_m_he    * w_out_he_en
            cheat_P_he   = W_out_he_en   - W_in_he_en
            resi         = self.P_he     - cheat_P_he
        # --- Check if match the maximum heat output --------------------------
        if  self.parameters['mode'] == 'sink':
            q_available  = self.i_he_cs_ex  - self.i_he_cs_su
            Q_available  = self.m_he_cs     * q_available
            q_required   = self.i_he_4r     - self.i_he_1
            m_required   = Q_available/q_required
            resi         = cheat_m_he       - m_required
        # --- Check if match the maximum heat input ---------------------------
        if  self.parameters['mode'] == 'source':
            q_available  = self.i_he_hs_su  - self.i_he_hs_ex
            Q_available  = self.m_he_hs     * q_available
            q_required   = self.i_he_3      - self.i_he_2r
            m_required   = Q_available/q_required
            resi         = cheat_m_he       - m_required
        return resi
    
    def find_flow_secondary_hs(self):
        """
        Finds the mass flow rate of the hot secondary fluid to match the 
        specified temperature profile based on the refrigerant flow rate.
        """
        q_available  = self.i_he_3      - self.i_he_2r
        Q_available  = self.m_he        * q_available
        q_required   = self.i_he_hs_su  - self.i_he_hs_ex
        if q_required > 0:  m_required = Q_available/q_required
        else:               m_required = 0
        return m_required
    
    def find_flow_secondary_cs(self):
        """
        Finds the mass flow rate of the cold secondary fluid to match the 
        specified temperature profile based on the refrigerant flow rate.
        """
        q_available  = self.i_he_4r    - self.i_he_1
        Q_available  = self.m_he       * q_available
        q_required   = self.i_he_cs_ex - self.i_he_cs_su
        if q_required > 0:  m_required = Q_available/q_required
        else:               m_required = 0
        return m_required
    
    def check_consistency(self):
        """
        Check that the results are consistent.
        """
        if self.s_he_4x > self.s_he_4:                  self.error = True
        if self.s_he_4r > self.s_he_4:                  self.error = True
        if self.s_he_4x > self.s_he_4r:                 self.error = True
        if self.s_he_1x < self.s_he_1:                  self.error = True
        if self.s_he_2  < self.s_he_1:                  self.error = True
        if self.s_he_2  > self.s_he_3:                  self.error = True
        if self.s_he_4  < self.s_he_3:                  self.error = True
        if self.s_he_4x < self.s_he_1x:                 self.error = True
        if self.T_he_2r > self.T_he_hs_ex-0.95*self.parameters['dT_he_ev_pp']: 
            self.error = True
        if self.T_he_3  > self.T_he_hs_su-0.95*self.parameters['dT_he_ev_pp']:
            self.error = True
        if self.T_he_4r < self.T_he_cs_ex+0.95*self.parameters['dT_he_cd_pp']: 
            self.error = True
        if self.T_he_1  < self.T_he_cs_su+0.95*self.parameters['dT_he_cd_pp']: 
            self.error = True
        if self.i_he_2r < self.i_he_2:                  self.error = True
        if self.eta_he_cyclen < 0:                      self.error = True
        if self.T_he_4x-self.T_ci < 0.95*self.parameters['dT_he_cd_pp']:    
            self.error = True 
        
        if self.p_he_1x < 0.99*self.p_he_1x_min: self.error = True
        if self.p_he_3  < 1.00*self.p_he_3_min:  self.error = True
        
        if abs(self.resi_1)   > 1: self.error = True
        if abs(self.resi_2)   > 1: self.error = True
        if abs(self.resi_rec) > 1: self.error = True
        
        if not bool(self.parameters['wet_ex']):
            p_vec = np.linspace(self.p_he_3,self.p_he_4)
            for i, p in enumerate(p_vec):
                self.state_he.update(CoolProp.PSmass_INPUTS,p,self.s_he_3) 
                h_is = self.state_he.hmass()
                h    = self.i_he_3-self.eta_is_ex*(self.i_he_3-h_is)
                self.state_he.update(CoolProp.HmassP_INPUTS,h,p)
                if 0 < self.state_he.Q() < 1:
                    self.error = True
                if  self.error: break
    
    def retrieve_kpi(self):
        """
        Retrieve the KPI's of the TRORC
        """
        
        if self.options['exergy']: self.compute_exergy()
        
        # --- Basic quantities ------------------------------------------------
        self.w_in_he_en = self.i_he_2 -self.i_he_1 
        self.q_in_he_en = self.i_he_3 -self.i_he_2r
        self.w_out_he_en= self.i_he_3 -self.i_he_4 
        self.q_out_he_en= self.i_he_4r-self.i_he_1
        self.w_net_he_en= self.w_out_he_en-self.w_in_he_en
        self.f_he_hs    = self.q_in_he_en /(self.i_he_hs_su-self.i_he_hs_ex)
        self.f_he_cs    = self.q_out_he_en/(self.i_he_cs_ex-self.i_he_cs_su)
        
        if self.options['exergy']:
            self.w_in_he_ex = self.e_he_2 -self.e_he_1 
            self.q_in_he_ex = self.e_he_3 -self.e_he_2r
            self.w_out_he_ex= self.e_he_3 -self.e_he_4 
            self.q_out_he_ex= self.e_he_4r-self.e_he_1 
            self.w_net_he_ex= self.w_out_he_ex-self.w_in_he_ex
        
        # --- Efficiencies ----------------------------------------------------
        self.eta_he_cyclen = self.w_net_he_en/self.q_in_he_en
        
        if self.options['exergy']:
            self.eta_he_cyclex = self.w_net_he_en/self.q_in_he_ex
            self.eta_he_pumpex = self.w_in_he_ex /self.w_in_he_en
            self.eta_he_expex  = self.w_out_he_en/self.w_out_he_ex
            self.eta_he_regex  = (self.e_he_2r-self.e_he_2) \
                               / (self.e_he_4-self.e_he_4r)
            self.eta_he_evapen = self.f_he_hs**(-1) * self.q_in_he_en\
                               / (self.i_he_hs_su-self.i_he_hs_ex)
            self.eta_he_evapex = self.f_he_hs**(-1) * self.q_in_he_ex\
                               / (self.e_he_hs_su-self.e_he_hs_ex)
            self.eta_he_srcen  = (self.i_he_hs_su-self.i_he_hs_ex) \
                               / (self.i_he_hs_su-self.i_he_hs_ref)
            self.eta_he_srcex  = (self.e_he_hs_su-self.e_he_hs_ex) \
                               / (self.e_he_hs_su-self.e_he_hs_ref)
            self.eta_he_genven = self.eta_he_evapen*self.eta_he_srcen
            self.eta_he_genvex = self.eta_he_evapex*self.eta_he_srcex  
            self.eta_he_condex = self.f_he_cs*(self.e_he_cs_ex\
                               - self.e_he_cs_su) / self.q_out_he_ex
            self.eta_he_toten  = self.eta_he_cyclen*self.eta_he_genven
            self.eta_he_totex  = self.eta_he_cyclex*self.eta_he_genvex
        
        # --- Losses ----------------------------------------------------------
        if self.options['exergy']:
            self.loss_he_pumpex = self.w_in_he_en-self.w_in_he_ex
            self.loss_he_expex  = self.w_out_he_ex-self.w_out_he_en
            self.loss_he_regex  = (self.e_he_4-self.e_he_4r) \
                                - (self.e_he_2r-self.e_he_2)
            self.loss_he_evapex = self.f_he_hs*(self.e_he_hs_su\
                                - self.e_he_hs_ex) - self.q_in_he_ex
            self.loss_he_srcen  = self.f_he_hs*(self.i_he_hs_ex\
                                - self.i_he_hs_ref)
            self.loss_he_srcex  = self.f_he_hs*(self.e_he_hs_ex\
                                - self.e_he_hs_ref)
            self.loss_he_condex = self.q_out_he_ex-self.f_he_cs \
                                * (self.e_he_cs_ex-self.e_he_cs_su)
            self.loss_he_sinken = self.f_he_cs*(self.i_he_cs_ex\
                                                -self.i_he_cs_su)
            self.loss_he_sinkex = self.f_he_cs*(self.e_he_cs_ex\
                                                -self.e_he_cs_su)
            self.pow_he_supplen = self.f_he_hs*(self.i_he_hs_su\
                                                -self.i_he_hs_ref)
            self.pow_he_supplex = self.f_he_hs*(self.e_he_hs_su\
                                                -self.e_he_hs_ref)
        
        # --- Carnot and Lorentz efficiencies ---------------------------------
        T_cs_logAv =(self.T_he_cs_ex-self.T_he_cs_su)\
             /np.log(self.T_he_cs_ex/self.T_he_cs_su)
        T_hs_logAv =(self.T_he_hs_su-self.T_he_hs_ex)\
             /np.log(self.T_he_hs_su/self.T_he_hs_ex)
        ETA_lorenz =(     T_hs_logAv-     T_cs_logAv)/     T_hs_logAv
        ETA_carnot =(self.T_he_hs_su-self.T_he_cs_su)/self.T_he_hs_su
        self.eta_he_lorenz = ETA_lorenz
        self.eta_he_carnot = ETA_carnot
        self.psi_he_lorenz = self.eta_he_cyclen/ETA_lorenz
        self.psi_he_carnot = self.eta_he_cyclen/ETA_carnot
        
        # --- Preliminary design information ----------------------------------
        self.vol_coef       = self.v_he_3/self.w_out_he_en
        self.exp_ratio      = self.p_he_3/self.p_he_4
        self.vol_ratio      = self.v_he_4/self.v_he_3
        
        if self.parameters['version'] == 'operational_light':
            # --- Multiply the specific quantities by the mass flow -----------
            self.Q_in_he_en = self.m_he*self.q_in_he_en
            self.W_in_he_en = self.m_he*self.w_in_he_en
            self.Q_out_he_en= self.m_he*self.q_out_he_en
            self.W_out_he_en= self.m_he*self.w_out_he_en
            self.W_net_he_en= self.m_he*self.w_net_he_en 
    
    def export_states(self):
        """
        Export the thermodynamic states of the SRORC necessary to draw
        the T-s and p-h diagrams.
        """
        self.p_he = self.p_he_1, self.p_he_1x,self.p_he_2, self.p_he_2r,\
                    self.p_he_3 ,self.p_he_4, self.p_he_4x,self.p_he_4r
        self.T_he = self.T_he_1, self.T_he_1x,self.T_he_2, self.T_he_2r,\
                    self.T_he_3, self.T_he_4, self.T_he_4x,self.T_he_4r
        self.i_he = self.i_he_1, self.i_he_1x,self.i_he_2, self.i_he_2r,\
                    self.i_he_3, self.i_he_4, self.i_he_4x,self.i_he_4r
        self.s_he = self.s_he_1, self.s_he_1x,self.s_he_2, self.s_he_2r,\
                    self.s_he_3, self.s_he_4, self.s_he_4x,self.s_he_4r
        self.x_he = self.x_he_1, self.x_he_1x,self.x_he_2, self.x_he_2r,\
                    self.x_he_3, self.x_he_4, self.x_he_4x,self.x_he_4r
        if self.options['exergy']: 
            self.e_he = self.e_he_1, self.e_he_1x,self.e_he_2, self.e_he_2r,\
                        self.e_he_3, self.e_he_4, self.e_he_4x,self.e_he_4r
        else:
            self.e_he = 0,0,0,0,0,0,0,0
        out = self.p_he,self.T_he,self.i_he,self.s_he,self.x_he,self.e_he
        return out    
    
    def compute_exergy(self):
        """
        Compute the exergy at each state of the TRORC.
        """
        self.compute_ref_states()
        
        # 1 - Exergy of the cycle ---------------------------------------------
        
        self.e_he_1  = (self.i_he_1 -self.i_he_0)\
          -self.T_he_0*(self.s_he_1 -self.s_he_0)
        self.e_he_1x = (self.i_he_1x-self.i_he_0)\
          -self.T_he_0*(self.s_he_1x-self.s_he_0)
        self.e_he_2  = (self.i_he_2 -self.i_he_0)\
          -self.T_he_0*(self.s_he_2 -self.s_he_0)
        self.e_he_2r = (self.i_he_2r-self.i_he_0)\
          -self.T_he_0*(self.s_he_2r-self.s_he_0)
        self.e_he_3  = (self.i_he_3 -self.i_he_0)\
          -self.T_he_0*(self.s_he_3 -self.s_he_0)
        self.e_he_4  = (self.i_he_4 -self.i_he_0)\
          -self.T_he_0*(self.s_he_4 -self.s_he_0)
        self.e_he_4x = (self.i_he_4x-self.i_he_0)\
          -self.T_he_0*(self.s_he_4x-self.s_he_0)
        self.e_he_4r = (self.i_he_4r-self.i_he_0)\
          -self.T_he_0*(self.s_he_4r-self.s_he_0)
        
        if round(self.T_he_1 ,2)==round(self.T_he_0,2): self.e_he_1  = 0
        if round(self.T_he_1x,2)==round(self.T_he_0,2): self.e_he_1x = 0
        if round(self.T_he_2 ,2)==round(self.T_he_0,2): self.e_he_2  = 0
        if round(self.T_he_2r,2)==round(self.T_he_0,2): self.e_he_2r = 0
        if round(self.T_he_3 ,2)==round(self.T_he_0,2): self.e_he_3  = 0
        if round(self.T_he_4 ,2)==round(self.T_he_0,2): self.e_he_4  = 0
        if round(self.T_he_4x,2)==round(self.T_he_0,2): self.e_he_4x = 0
        if round(self.T_he_4r,2)==round(self.T_he_0,2): self.e_he_4r = 0
        
        # 2 - Exergy of the boundaries ----------------------------------------
        
        self.compute_exergy_boundaries()