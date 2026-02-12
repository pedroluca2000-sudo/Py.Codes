from CoolProp.CoolProp import PropsSI
FLUID="NitrousOxide"

def get_rho(Pressure, Temperature):
    rho = PropsSI("D", "P", Pressure, "T", Temperature, FLUID)

    return rho 

enthalpy = PropsSI("H", "P", 50e5, "Q", 0, FLUID)
entropy = PropsSI("S", "P", 1e5, "Q", 0, FLUID)
rho_exp = get_rho(50e5, 280)
print(f"rho={rho_exp} kg/m3 \nenthalpy={enthalpy/1e3} kJ/kg \nentropy={entropy/1e3} kJ/K")
