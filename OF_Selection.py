import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj, add_new_fuel

#USER SETTINGS
Pc = 350.0        # chamber pressure [psia]
OF_range = np.linspace(1.0, 10, 60)

#CREATE E98 FUEL
e98_card = """fuel C2H6 wt=0.98 fuel H2O  wt=0.02"""
add_new_fuel("E98", e98_card)

#CREATE CEA OBJECT
cea = CEA_Obj(oxName="N2O", fuelName="E98", fac_CR=None)  # E98 case

isp_list = []
tc_list = []
eps_list = []

for OF in OF_range:
    PcOvPe = Pc / 14.7

    eps = cea.get_eps_at_PcOvPe(Pc=Pc, MR=OF, PcOvPe=PcOvPe)
    isp = cea.get_Isp(Pc=Pc, MR=OF, eps=eps, frozen=0)
    tc  = cea.get_Tcomb(Pc=Pc, MR=OF)
    
    eps_list.append(eps)
    isp_list.append(isp)
    tc_list.append(tc)

#PLOTTING
plt.figure()
plt.plot(OF_range, isp_list)

plt.axvline(x=9, color="red")
plt.axvline(x=3, color="green")
plt.xlabel("O/F Ratio")
plt.ylabel("Vacuum Isp (s)")
plt.title("Isp vs O/F (N2O + E98)")
plt.grid()

plt.figure()
plt.plot(OF_range, tc_list)
plt.axvline(x=9, color="red")
plt.axvline(x=3, color="green")
plt.axhline(y=1743, color="blue")
plt.xlabel("O/F Ratio")
plt.ylabel("Chamber Temperature (K)")
plt.title("Tc vs O/F (N2O + E98)")
plt.grid()

plt.figure()
plt.plot(OF_range, eps_list)
plt.axvline(x=3, color="green")
plt.xlabel("O/F")
plt.ylabel("Required Expansion Ratio")
plt.grid()

plt.show()

