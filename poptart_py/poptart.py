import numpy as np
import collections

Atmosphere_layer = collections.namedtuple("Atmosphere_layer", 't_a p_a h_rel d')

# snell
def snell(theta1, v1, v2):
    return np.arcsin(v2*np.sin(theta1)/v1)

# speed of sound as a function of temperature
def sound_speed(t_a):
    return 331 + 0.6 * t_a

# calculate sound attenuation
def calc_sound_attenuation(f, atm_layer):
    t_a = atm_layer.t_a
    p_a = atm_layer.p_a
    h_rel = atm_layer.h_rel
    p_r = 101.325
    t_r = 293.15
    t_01 = 273.16
    v = 10.79586*(1-t_01/t_a)-5.02808*np.log10(t_a/t_01) \
        +1.50474*1e-4*(1-10**(-8.29692*(t_a/t_01-1))) \
        +0.42873e-3*1e-3*(-1+10**(4.76955*(1-t_01/t_a)))-2.2195983
    p_sat = p_r*10**v
    h = h_rel*(p_sat/p_r)*(p_a/p_r)**-1
    f_ro = p_a*(24+(4.04e4*h)*(0.02+h)/(0.391+h))
    f_rn = p_a*(t_a/t_r)**-0.5*(9+280*h*np.exp(-4.17*(t_a/t_r)**-1/3-1))
    a_f = 8.686*f**2*((1.84*1e-11*(p_a/p_r)**-1*(t_a/t_r)**0.5)
                      +(t_a/t_r)**-2.5*(0.01278*np.exp(-2239.1/t_a)
                                        *f_ro/(f_ro**2+f**2))
                      +0.1068*np.exp(-3352.0/t_a)*f_rn/(f_rn**2+f**2))
    return a_f
