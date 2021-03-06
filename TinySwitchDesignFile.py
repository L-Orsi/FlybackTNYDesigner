import numpy as np
from scipy.optimize import fsolve
import pprint

# Input Voltage
vac = 220  # input altern voltage [V]
vac_tol = 20  # input voltage tolerance [V]
fl = 50  # line frequency [Hz]

# Output Diode
vD = 0.5  # diode forward voltage [V]
vr = 40  # diode max reverse voltage [V]

# Output
vo = 5  # output voltage [V]
po = 15  # output power [W]

# Snubber Topology
snubber = ["RCD", "RC"][0]  # Snubber topology

# TNY Specs
i_limit = 0.512  # From TNY268 datasheet Min Value - [A]
ip2 = 0.588  # From TNY268 datasheet Max Value - [A]
fs = 132e3 * 0.94  # Min Switching Frequency [Hz]

# Transformer Properties
bp = 2000  # [Gauss | mT*10]
Ae = 0.6  # [cm2 | mm2/100]
k1 = 90  # [Adimensional]
k2 = -0.708  # [Adimensional]

# Optocoupler
vled = 1  # [V]


# Script


vac_min = vac*(1-vac_tol/100)
vac_max = vac*(1+vac_tol/100)
diode_loss = vD/vo*100
snubber_loss = 0.15 if snubber == "RCD" else 0.2
nu = (100 - (diode_loss + snubber_loss))/100
vd_max = np.sqrt(2)*vac_max
ci = po * 3e-6
tc = 3e-3
vd_min = ((2*vac_min**2)-(2*po*(1/(2*fl)-tc))/(nu*ci))**0.5
piv = 0.8*vr  # peak inverse voltage estimation [V]
vor = vd_max*(vo+vD)/(piv-vo)  # refleted output volage [V]
if vor > 150:
    print("Choose a Different Diode")
ip = i_limit*0.9
d_max = 2*po/(nu*vd_min*ip)
kdp = vor*(1-d_max)/(vd_min*d_max)
if kdp > 1:
    print(f"Mostly DCM : kdp={kdp}")
elif 0 < kdp < 1:
    ##print(f"Mostly CCM : kdp={kdp}")
    pass
else:
    print(f"Error - Not DCM nor CCM : kdp={kdp}")
if kdp > (1-d_max)/(0.67-d_max):
    print("Fully Discontinuos!")
    exit  # Not in the scope of this script
else:
    ##print("CCM Design Started!")
    pass
d_max = vor / (vor + vd_min)
krp = 2*(ip*d_max*nu*vd_min-po)/(ip*d_max*nu*vd_min)
if krp >= 0.6:
    #print(f"Good To Go! : krp={krp}")
    pass
else:
    print("Not in the scope of this script!")
    exit
Z = 1  # Loss Allocation Factor
Lp = (po*1e6)/(krp*(1-krp/2)*1/0.9*ip**2*fs)*(Z*(1-nu)+nu)/nu
Np = np.round(100 * ip2 * Lp/(bp*Ae))
Ns = np.round(Np*(vo+vD)/vor)


def func(x): return 40*np.pi*Ae*(Np**2/(1000*Lp)-1/(k1*x**k2))


Lg_guess = Np**2*1700e-9
print(Lg_guess)
Lg = fsolve(func, Lg_guess)[0]
Al = (k1*Lg**k2)

if not (0.1 < Lg < 2 and 60 < Al < 560):
    print("Gap formula may be innacurate")

iprms = ip2*(d_max*(krp**2/3-krp+1))**0.5
isrms = ip2*Np/Ns*((1-d_max)*(krp**2/3-krp+1))**0.5
io = po/vo  # Not Sure
iripple = (isrms**2-io**2)**0.5
vz = vo - vled  # Use TL431 for better accuracy

res = {
    "Max Magnetic Field B [mT]": bp/10,
    "Line Frequency [Hz]": 50,
    # "Output Ripple Current [A]": iripple,
    "Output Voltage [V]": vo,
    "Input VAC Voltage [Vrms]": vac,
    "Primary Inductance [uH]": Lp,
    "Primary Turns [N]": Np,
    "Secondary Turns [N]": Ns,
    "Gap Length [mm]": Lg,
    "Switching Frequency [kHz]": fs/1000,
    "RMS Primary Current [A]": iprms,
    "Rectifier Diodes Max Current [A]": iprms*1.5,
    "Max Duty Cycle [%]": d_max*100,
    "RMS Secondary Current [A]": isrms,
    "Control Zener Voltage [V]": vz,
    "Input Capacitor [uF]": ci,
    "Min Input Voltage [Vdc]": vd_min
}

pprint.pprint(res, indent=2)
