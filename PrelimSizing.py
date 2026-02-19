import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj

#USER SETTINGS
Pc = 350.0        # chamber pressure [bar]
OF_range = np.linspace(1.0, 10, 60)

#CREATE CEA OBJECT
ceaIPA = CEA_Obj(oxName="N2O", fuelName="C3H8", fac_CR=None)  # IPA case
ceaE98 = CEA_Obj(oxName="N2O", fuelName="C2H6", fac_CR=None)  # E98 case

IPA_isp_list = []
IPA_tc_list = []
IPA_eps_list = []

E98_isp_list = []
E98_tc_list = []
E98_eps_list = []

for OF in OF_range:
    PcOvPe = Pc / 14.7
    IPA_eps = ceaIPA.get_eps_at_PcOvPe(Pc=Pc, MR=OF, PcOvPe=PcOvPe)
    IPA_isp = ceaIPA.get_Isp(Pc=Pc, MR=OF, eps=IPA_eps, frozen=0)
    IPA_tc  = ceaIPA.get_Tcomb(Pc=Pc, MR=OF)

    E98_eps = ceaE98.get_eps_at_PcOvPe(Pc=Pc, MR=OF, PcOvPe=PcOvPe)
    E98_isp = ceaE98.get_Isp(Pc=Pc, MR=OF, eps=E98_eps, frozen=0)
    E98_tc  = ceaE98.get_Tcomb(Pc=Pc, MR=OF)
    
    IPA_isp_list.append(IPA_isp)
    IPA_tc_list.append(IPA_tc)
    IPA_eps_list.append(IPA_eps)
    
    E98_isp_list.append(E98_isp)
    E98_tc_list.append(E98_tc)
    E98_eps_list.append(E98_eps)

#PLOTTING
plt.figure()
plt.plot(OF_range, IPA_isp_list, color = "orange")
plt.plot(OF_range, E98_isp_list, color = "blue")
plt.axvline(x=9, color="red")
plt.xlabel("O/F Ratio")
plt.ylabel("Vacuum Isp (s)")
plt.title("Isp vs O/F (N2O + IPA)")
plt.grid()

plt.figure()
plt.plot(OF_range, IPA_tc_list, color = "orange")
plt.plot(OF_range, E98_tc_list, color = "blue")
plt.axvline(x=9, color="red")
#plt.axhline(y=1743, color="blue")
plt.xlabel("O/F Ratio")
plt.ylabel("Chamber Temperature (K)")
plt.title("Tc vs O/F (N2O + IPA)")
plt.grid()

plt.figure()
plt.plot(OF_range, IPA_eps_list)
plt.plot(OF_range, E98_eps_list)
plt.xlabel("O/F")
plt.ylabel("Required Expansion Ratio")
plt.grid()

plt.show()
