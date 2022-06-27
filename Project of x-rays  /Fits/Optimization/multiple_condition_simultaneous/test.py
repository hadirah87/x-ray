import numpy as np
import func

t=45 #grating thickness
code='sample.wls' #name of the mathematica code
d=130 #periodic length of grating
rho=2.06 #sld
w=0.8 #wavelength in nm
X=[d,rho,t,w]
func.simulation4(code,d,rho,t,w)
