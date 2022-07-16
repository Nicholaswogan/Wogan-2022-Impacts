import numpy as np
from pathos.multiprocessing import ProcessingPool as Pool

from ImpactAtmosphere import SteamAtm
from coupling import output2photochem
from photochem import Atmosphere, zahnle_earth

class Constants:
    yr = 365*24*60*60
cons = Constants()

def impact_evolve(init, settings_in, outfile, eddy, ztop, nz, zero_out, rainfall_rate, atol, rtol, t_eval):
    N_H2O_ocean = init['N_H2O_ocean']
    N_CO2 = init['N_CO2']
    N_N2  = init['N_N2']
    M_i = init['M_i']
    stm = SteamAtm('zahnle_earth_ct.yaml')
    sol_stm = stm.impact(N_H2O_ocean,N_CO2,N_N2,M_i)

    settings_out = outfile+"_settings.yaml"
    atmosphere_out = outfile+"_atmosphere.txt"

    output2photochem(stm, sol_stm, settings_in, settings_out, atmosphere_out, eddy, ztop, nz, zero_out, rainfall_rate)
    
    pc = Atmosphere(zahnle_earth,\
                    settings_out,\
                    "input/Sun_4.0Ga.txt",\
                    atmosphere_out)

    pc.var.atol = atol  
    pc.var.rtol = rtol      
    t_start = 0.0
    success = pc.evolve(outfile+'.dat',t_start, pc.wrk.usol, t_eval, overwrite=True)

def Ceres_nominal():
    params = {}

    init = {}
    init['N_H2O_ocean'] = 15.0e3
    init['N_CO2'] = 23.*0.1
    init['N_N2'] = 36.
    init['M_i'] = 2e24
    params['init'] = init

    params['settings_in'] = "input/settings_Hadean.yaml"
    params['outfile'] = "results/CO2=2.3e0_N2=3.6e1_M_i=2.0e24_eddy=1e6"
    params['eddy'] = 1e6
    params['ztop'] = 2200e5
    params['nz'] = 200
    params['zero_out'] = ['NH3']
    params['rainfall_rate'] = -1

    params['atol'] = 1e-27
    params['rtol'] = 1e-3
    params['t_eval'] = np.logspace(5,np.log10(cons.yr*10e6),500)

    return params

def Ceres_rain():
    params = {}

    init = {}
    init['N_H2O_ocean'] = 15.0e3
    init['N_CO2'] = 23.*0.1
    init['N_N2'] = 36.
    init['M_i'] = 2e24
    params['init'] = init

    params['settings_in'] = "input/settings_Hadean.yaml"
    params['outfile'] = "results/CO2=2.3e0_N2=3.6e1_M_i=2.0e24_eddy=1e6_rain"
    params['eddy'] = 1e6
    params['ztop'] = 2200e5
    params['nz'] = 200
    params['zero_out'] = []
    params['rainfall_rate'] = 1

    params['atol'] = 1e-27
    params['rtol'] = 1e-3
    params['t_eval'] = np.logspace(5,np.log10(cons.yr*10e6),500)

    return params

    
if __name__ == "__main__":
    impact_evolve(**Ceres_rain())

    # def wrap(fun):
    #     impact_evolve(**fun())
    # simulations = [Ceres_nominal, Ceres_rain]
    # p = Pool(2)
    # p.map(wrap, simulations)

    
    
    
    
    
