# -*- coding: utf-8 -*-
"""
Created on Sat Nov  3 12:11:13 2018

@author: kpett
""" 
import cantera as ct
from matplotlib import pyplot
from pandas import DataFrame

#Function returns reversible work and irriversibility of a process
def rev_irrev(hin, hout, sin, sout, Tin, Qin, Qout, mdotratio):
    To = 300
    irrev = mdotratio*(Qout+To*(sout-sin))
    rev = ((hin-hout)-To*(sin-sout))*mdotratio+Qin*(1-(To/Tin))*mdotratio
    return rev, irrev
  
#Function returns the actual enthalpy at the outlet of a pump     
def h_OutPump(n_pump, h_OutIs, h_In):
    h_OutAct = ((h_OutIs - h_In)/n_pump)+h_In
    return h_OutAct

#Function returns the actual enthalpy at the outlet of a compressor
def h_OutCompressor(n_compressor, h_OutIs, h_In):
    h_OutAct = ((h_OutIs - h_In)/n_compressor)+h_In
    return h_OutAct

#Function returns the actual enthalpy at the outlet of a turbine 
def h_OutTurbine(n_turb, h_OutIs, h_In):
    h_OutAct = -(n_turb)*(h_In-h_OutIs)+h_In
    return h_OutAct

#Define Lists
_pr = []    
_mDotRatio = []
_n_thermal = []
_Wnet = []

_W12=[]
_Q23=[]
_W34=[]
_Q41=[]
_W56=[]
_Q67=[]
_W78=[]
_Q89=[]

_nII12 = []
_nII34 = []
_nII56 = []
_nII78 = []

_i12 = []
_i23 = []
_i34 = []
_i41 = []
_i56 = []
_i67 = []
_i78 = []
_i89 = []

#Define Fluid States
air5 = ct.Solution('air.cti')
air6 = ct.Solution('air.cti')
air7 = ct.Solution('air.cti')
air8 = ct.Solution('air.cti')
air9 = ct.Solution('air.cti')
water1 = ct.Water()
water2 = ct.Water()
water3 = ct.Water()
water4 = ct.Water()

#Define Efficiencies
n_compressor = 0.8
n_turbineAir = 0.85
e_HRSG = 0.86
n_pump = 0.9
n_turb_w = 0.9

for pr in range(3,21):
    
    "State 5 - AIR - Inlet to Compressor"
    P5 = 101325
    T5 = 300             
    air5.TP = T5, P5    
    s5 = air5.s
    h5 = air5.h
    
    "State 6 - AIR - Outlet of Compressor/Inlet to Combustion Chamber"
    P6 = 101325*pr
    s6_is = s5
    air6.SP = s6_is, P6
    h6_is= air6.h
    h6 = h_OutCompressor(n_compressor, h6_is, h5)
    air6.HP = h6,P6
    s6 = air6.s
    
    "State 7 - AIR - Outlet of Combustion Chamber/Inlet to Turbine"    
    P7 = P6
    T7 = 1400
    air7.TP = T7,P7
    s7 = air7.s
    h7 = air7.h
     
    "State 8 - AIR - Outlet of Turbine/Inlet to HRSG"
    P8 = P5
    s8_is = s7
    air8.SP = s8_is, P8
    h8_is = air8.h
    h8 = h_OutTurbine(n_turbineAir, h8_is, h7)
    air8.HP = h8,P8
    T8 = air8.T
    
    "State 9 - AIR - Outlet of HRSG"
    T9 = 450
    P9 = P8
    air9.TP = T9, P9
    h9 = air9.h
        
    "State 1 - WATER - Outlet of Condenser/Inlet to pump"
    P1 = 5*10**3
    water1.PX = P1, 0
    h1 = water1.h
    s1 = water1.s
    
    "State 2 - WATER - Outlet of Pump/Inlet to HRSG"
    P2 = 7*10**6
    s2_is = s1
    water2.SP = s2_is, P2
    h2_is = water2.h
    h2 = h_OutPump(n_pump, h2_is, h1)
    water2.HP = h2,P2
    T2 = water2.T
    
    "State 3 - WATER - Outlet of HRSG/Inlet to Turbine"
    P3 = P2
    T3Perf = T8
    water3.TP = T3Perf,P3
    h3Perf = water3.h
    h3 = e_HRSG*(h3Perf-h2)+h2     
    water3.HP = h3, P3               #h3-h2/h3perf-h2 = effectiveness
    mDotRatio = (h8-h9)/(h3-h2)
    s3 = water3.s
    
    "State 4 -WATER - Outlet of Turbine/Inlet to Condenser"
    P4 = P1
    s4_is =  s3
    water4.SP = s4_is,P4
    h4_is = water4.h
    h4 = h_OutTurbine(n_turb_w, h4_is, h3)
    water4.HP = h4, P4
    
    #Define Remaining Variables
    T1 = water1.T
    T3 = water3.T
    T4 = water4.T
    T6 = air6.T
    s2 = water2.s
    s3 = water3.s
    s4 = water4.s
    s6 = air6.s
    s8 = air8.s
    s9 = air9.s
    
    #Works and Heats
    W12 = mDotRatio*(h1-h2)    # Work into Rankine pump per air mass flow rate (J/kg)                (-)
    Q23 = -mDotRatio*(h2-h3)   # Heat into water in HRSG per air mass flow rate (J/kg)               (+)
    W34 = mDotRatio*(h3-h4)    # Work out of Rankine turbine per air mass flow rate (J/kg)           (+)
    Q41 = -mDotRatio*(h4-h1)   # Heat out of water in condenser per air mass flow rate (J/kg)        (-)
    W56 = h5-h6                # Work into Brayton compressor per air mass flow rate (J/kg)          (-)   
    Q67 = -(h6-h7)             # Heat into combustion chamber Brayton per air mass flow rate (J/kg)  (+)
    W78 = h7-h8                # Work out of Brayton turbine per air mass flow rate (J/kg)           (+)
    Q89 = -(h8-h9)             # Heat out of air in HRSG per air mass flow rate (J/kg)               (-)
    
    Wnet = W56+W78+W12+W34    # Net work out of Cogen cycle
    n_thermal = Wnet/Q67      # Thermal efficiency of Cogen cycle
    
    #Reversible work and irreversibility
    Wrev12, i12 = rev_irrev(h1, h2, s1, s2, T1, 0, 0, mDotRatio)         # Pump - Rankine
    Wrev23, i23 = rev_irrev(h2, h3, s2, s3, T2, (h3-h2), 0, mDotRatio)   # HRSG - Rankine
    Wrev34, i34 = rev_irrev(h3, h4, s3, s4, T3, 0, 0, mDotRatio)         # Turbine - Rankine
    Wrev41, i41 = rev_irrev(h4, h1, s4, s1, T4, 0, -Q41, mDotRatio)      # Condenser - Rankine
    Wrev56, i56 = rev_irrev(h5, h6, s5, s6, T5, 0, 0, 1)         # Compressor - Brayton
    Wrev67, i67 = rev_irrev(h6, h7, s6, s7, T6, Q67, 0, 1)       # Combustion Chamber - Brayton
    Wrev78, i78 = rev_irrev(h7, h8, s7, s8, T7, 0, 0, 1)         # Turbine - Brayton
    Wrev89, i89 = rev_irrev(h8, h9, s8, s9, T8, 0, (h8-h9), 1)   # HRSG - Brayton
    
    #Second Law Efficiencies
    nII12 = Wrev12/W12    # Pump - Rankine
    nII34 = W34/Wrev34    # Turbine - Rankine
    nII56 = Wrev56/W56    # Compressor - Brayton
    nII78 = W78/Wrev78    # Turbine - Brayton
    
    #Lists for plots
    _pr.append(pr)    
    _mDotRatio.append(mDotRatio)
    _n_thermal.append(n_thermal)
    _Wnet.append(Wnet)
    
    _nII12.append(nII12)
    _nII34.append(nII34)
    _nII56.append(nII56)
    _nII78.append(nII78)
    
    _i12.append(i12)
    _i23.append(i23)
    _i34.append(i34)
    _i41.append(i41)
    _i56.append(i56)
    _i67.append(i67)
    _i78.append(i78)
    _i89.append(i89)
    
    _W12.append(W12)
    _Q23.append(Q23)
    _W34.append(W34)
    _Q41.append(Q41)
    _W56.append(W56)
    _Q67.append(Q67)
    _W78.append(W78)
    _Q89.append(Q89)
  
pyplot.figure('Thermal Efficiency v. Pressure Ratio')
pyplot.plot(_pr, _n_thermal)
pyplot.xlabel('Pressure Ratio')
pyplot.ylabel('Thermal Efficiency')
pyplot.title('Thermal Efficiency v. Pressure Ratio')
pyplot.tight_layout()
pyplot.savefig('ThermalEffvPR.jpg')

pyplot.figure('Mass Flow Ratio v. Pressure Ratio')
pyplot.plot(_pr, _mDotRatio)
pyplot.xlabel('Pressure Ratio')
pyplot.ylabel('Mass Flow Ratio')
pyplot.title('Mass Flow Ratio v. Pressure Ratio')
pyplot.tight_layout()
pyplot.savefig('MdotvPR.jpg')

pyplot.figure('Net Power per Mass Flow Rate of Air v. Pressure Ratio')
pyplot.plot(_pr, _Wnet)
pyplot.xlabel('Pressure Ratio')
pyplot.ylabel('Net Power per Mass Flow Rate of Air $(J/kg)$')
pyplot.title('Net Power per Mass Flow Rate of Air v. Pressure Ratio')
pyplot.tight_layout()
pyplot.savefig('NetWorkvPR.jpg')

pyplot.figure('II Law Efficiencies')
pyplot.plot(_pr, _nII78, label="Gas Turbine")
pyplot.plot(_pr, _nII56, label="Compressor")
pyplot.plot(_pr, _nII34, label="Steam Turbine")
pyplot.plot(_pr, _nII12, label="Pump")
pyplot.legend()
pyplot.xlabel('Pressure Ratio')
pyplot.ylabel('II Law Efficiencies')
pyplot.title('II Law Efficiencies v. Pressure Ratio')
pyplot.tight_layout()
pyplot.savefig('IILawEff.jpg')

fig, ax1 = pyplot.subplots()
ax1.plot(_pr, _Wnet, lw=2, color="blue")
ax1.set_ylabel("Net Power per Air Mass Flow Rate $(J/kg)$", fontsize=10, color="blue")
for label in ax1.get_yticklabels():
    label.set_color("blue")
ax2 = ax1.twinx()
ax2.plot(_pr, _n_thermal, lw=2, color="red")
ax2.set_ylabel("Thermal Efficiency", fontsize=10, color="red")
for label in ax2.get_yticklabels():
    label.set_color("red")
    
best_pr = _n_thermal.index(max(_n_thermal))+3
print('The max effeciency is: ', max(_n_thermal))
print('The most effecient Pressure Ratio is: ', best_pr)

_process = ["Compressor", "Combustion Chamber", "Gas Turbine", "Gas HRSG", "Pump", "Steam HRSG", "Steam Turbine", "Condensor"]
_i_best_pr = [_i12[best_pr],_i23[best_pr],_i34[best_pr],_i41[best_pr],_i56[best_pr],_i67[best_pr],_i78[best_pr],_i89[best_pr]]
_delta_h_best_pr = [-_W12[best_pr],_Q23[best_pr],-_W34[best_pr],_Q41[best_pr],-_W56[best_pr],_Q67[best_pr],-_W78[best_pr],_Q89[best_pr]]
_Q_best_pr = [0, _Q23[best_pr], 0, _Q41[best_pr], 0, _Q67[best_pr], 0, _Q89[best_pr]]
_W_best_pr = [_W12[best_pr], 0, _W34[best_pr], 0,_W56[best_pr], 0, _W78[best_pr], 0]

partITable = DataFrame({'Net Power Output per Air Mass Flow Rate (kJ/kg)': [x / 1000 for x in _Wnet],'Efficiency': _n_thermal, 'Mass Flow Ratio': _mDotRatio,'Pressure Ratio': _pr})
partITable.to_excel('Part1Table.xlsx', sheet_name = 'sheet1', index=False)

partIIabTable = DataFrame({ 'Irreversibility per Air Mass Flow Rate (kJ/kg)': [x / 1000 for x in _i_best_pr],'Change in Enthalpy per Air Mass Flow Rate (kJ/kg)': [x / 1000 for x in _delta_h_best_pr], 'Work per Air Mass Flow Rate (kJ/kg)': [x / 1000 for x in _W_best_pr], 'Heat Transfer per Air Mass Flow Rate (kJ/kg)': [x / 1000 for x in _Q_best_pr], 'Process': _process})
partIIabTable.to_excel('Part2Table.xlsx', sheet_name = 'sheet1', index=False)   
