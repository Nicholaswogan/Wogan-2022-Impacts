import numpy as np
from scipy.optimize import root
from photochem.utils._format import FormatSettings_main, yaml, MyDumper, Loader

def pahlevan_poly(P):
    if P > 0.21:
        log10P = np.log10(P)
        T = 295.74510905 + 76.02965729*log10P + \
            (-3.60269118)*log10P**2.0 + 9.88900334*log10P**3.0
    else:
        T = 235.0
    return T

def sat_pressure_H2O(T):
    return np.exp(14.8 -5420.0/T -0.00795/T**2)*1e6


def gravity(radius, mass, z):
    G_grav = 6.67430e-11

    grav = np.empty(len(z))
    for i in range(len(z)):         
        grav[i] = G_grav * (mass/1.0e3) / ((radius + z[i])/1.0e2)**2.0
        grav[i] = grav[i]*1.0e2 # convert to cgs

    return grav
    
# @nb.njit
def obj_all(x, nz, z, mubar_dry, Ntot_dry, grav):
    P0 = 10.0**x[0]
    PH2O = 10.0**x[1:nz+1]
    
    k = 1.3807e-16
    Navo = 6.022e23
    mu_H2O = 18
    P_trop = 0.21 # bar
    
    dz = z[1] - z[0]
    
    P = np.empty((nz))
    TTT = pahlevan_poly(P0/1e6)
    mubar = np.empty((nz))
    mubar[0] = (1 - PH2O[0]/P0)*mubar_dry + (PH2O[0]/P0)*mu_H2O
    P[0] = P0*np.exp(- mubar[0]*grav[0]/(Navo*k*TTT)*(dz/2))
    for i in range(1,nz):
        TTT = pahlevan_poly(P[i-1]/1e6)
        mubar[i] = (1 - PH2O[i-1]/max(P[i-1],1.0e-100))*mubar_dry + (PH2O[i-1]/max(P[i-1],1.0e-100))*mu_H2O
        P[i] = P[i-1]*np.exp(- mubar[i]*grav[i]/(Navo*k*TTT)*(dz))
    
    T = np.array([pahlevan_poly(PP/1e6) for PP in P])
    
    trop_ind = np.argmin(np.abs(P-P_trop*1e6))
    PH2O_new = np.empty(nz)
    for i in range(trop_ind+1):
        PH2O_new[i] = sat_pressure_H2O(T[i])
    for i in range(trop_ind+1,nz):
        PH2O_new[i] = (PH2O_new[trop_ind]/P[trop_ind])*P[i]
    
    nH2O = PH2O/(k*T)
    NH2O = np.sum(nH2O*dz)/Navo
    
    n = P/(k*T)
    N = n*dz
    Ntot = np.sum(N)
    weight = N/Ntot
    
    P0_new = (NH2O + Ntot_dry)*np.sum(mubar*weight)*np.sum(grav*weight)
    
    resid = np.empty(nz+1)
    resid[0] = np.log10(P0)-np.log10(P0_new)
    resid[1:nz+1] = np.log10(PH2O)-np.log10(PH2O_new)
    
    return resid, P, T

# @nb.njit
def obj(x, nz, z, mubar_dry, Ntot_dry, grav):
    resid, P, T = obj_all(x, nz, z, mubar_dry, Ntot_dry, grav)
    return resid

def make_water_profile(stm, sol):

    mass = 5.972e27
    radius = 6.371e8
    nz = 200    
    ztop = 500e5

    Ntot = sol['Ntot'][-1]
    Ntot_dry = Ntot - Ntot*sol['H2O'][-1]
    usol_dry = np.array([sol[sp][-1]*(Ntot/Ntot_dry) for sp in stm.gas.species_names])
    H2O_ind = stm.gas.species_names.index('H2O')
    usol_dry[H2O_ind] = 0.0
    mubar_dry = np.sum(usol_dry*stm.gas.molecular_weights)

    dz = ztop/nz
    z = np.empty(nz)
    z[0] = dz/2
    for i in range(1,nz):
        z[i] = z[i-1]+dz

    grav = gravity(radius, mass, z)

    init_cond = np.append(5,np.ones(nz)*-5)
    
    s = root(obj,init_cond,args=(nz, z, mubar_dry, Ntot_dry, grav,))
    if not s.success:
        raise Exception("root solve failed")

    resid, P, T = obj_all(s.x,nz, z, mubar_dry, Ntot_dry, grav)
    if (P[-1]/1e6) > 0.21:
        raise Exception("!!!")

    P0 = 10.0**s.x[0]
    PH2O = 10.0**s.x[1:nz+1]
    mubar = (1 - PH2O/P)*mubar_dry + (PH2O/P)*18

    usol = np.empty((len(stm.gas.species_names),nz))
    for i in range(nz):
        usol_dry[H2O_ind] = PH2O[i]/P[i]
        usol[:,i] = usol_dry
        
    Tsurf = pahlevan_poly(P0/1e6)
    PH2O_surf = sat_pressure_H2O(Tsurf)
    
    return z, T, P, usol, P0, PH2O_surf
    
def make_atmosphere_txt(z, T, usol, eddy, species_names, filename, zero_out):
    fil = open(filename,'w')
    fil.write('{:25}'.format('alt'))
    fil.write('{:25}'.format('temp'))
    fil.write('{:25}'.format('eddy'))
    for i,sp in enumerate(species_names):
        fil.write('{:25}'.format(sp))
    fil.write('\n')
    for j in range(len(z)):
        fil.write('{:25}'.format('%.15e'%(z[j]/1e5)))
        fil.write('{:25}'.format('%.15e'%T[j]))
        fil.write('{:25}'.format('%.15e'%eddy))
        for i,sp in enumerate(species_names):
            if sp in zero_out:
                fil.write('{:25}'.format('%.15e'%0.0))
            else:
                fil.write('{:25}'.format('%.15e'%usol[i,j]))
        fil.write('\n')
    fil.close()
    
    
def make_settings(infile, outfile, ztop, nz, PH2O_surf, P0, z, P, rainfall_rate):
    fH2O_surf = PH2O_surf/P0

    fil = open(infile,'r')
    data = yaml.load(fil,Loader=Loader)
    fil.close()

    data['atmosphere-grid']['bottom'] = 0.0
    data['atmosphere-grid']['top'] = float(ztop)
    data['atmosphere-grid']['number-of-layers'] = int(nz)

    data['planet']['surface-pressure'] = float(P0/1e6)
    
    if rainfall_rate > 0:
        data['planet']['water']['gas-rainout'] = True
        data['planet']['water']['rainfall-rate'] = rainfall_rate
        ind = np.argmin(np.abs(P - 0.21e6))
        data['planet']['water']['tropopause-altitude'] = float(z[ind])
    else:
        data['planet']['water']['gas-rainout'] = False    
    
    nn = len(data['boundary-conditions'])
    found = False
    for i in range(nn):
        if data['boundary-conditions'][i]['name'] == "H2O":
            found = True
            break

    H2O_boundary = {'type':'mix','mix':float(fH2O_surf)}
    if found:
        data['boundary-conditions'][i]['lower-boundary'] = H2O_boundary
    else:
        raise Exception("!!!")
        
    data = FormatSettings_main(data)

    fil = open(outfile,'w')
    yaml.dump(data,fil,Dumper=MyDumper,sort_keys=False,width=70)
    fil.close()
    
def output2photochem(stm, sol, settings_in, settings_out, atmosphere_out, eddy, ztop, nz, zero_out, rainfall_rate):
    z, T, P, usol, P0, PH2O_surf = make_water_profile(stm, sol)
    make_atmosphere_txt(z, T, usol, eddy, stm.gas.species_names, atmosphere_out, zero_out)
    make_settings(settings_in, settings_out, ztop, nz, PH2O_surf, P0, z, P, rainfall_rate)
    
    
    