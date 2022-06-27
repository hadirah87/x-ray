import os 
import numpy as np

def simulation(code_name,t,tc,d,f):
    os.system('wolframscript %s %s %s %s %s'%(code_name,t,tc,d,f))
    return

def simulation3(code_name,d,rho,t):
    os.system('wolframscript %s %s %s %s'%(code_name,d,rho,t))
    return

def simulation4(code_name,d,rho,t,w):
    os.system('wolframscript %s %s %s %s %s'%(code_name,d,rho,t,w))
    return

def simulation5(code_name,d,fd,fbd,ftd,t):
    os.system('wolframscript %s %s %s %s %s %s'%(code_name,d,fd,fbd,ftd,t))
    return

def simulation6(code_name,d,fd,fbd,ftd,t,tc):
    os.system('wolframscript %s %s %s %s %s %s %s'%(code_name,d,fd,fbd,ftd,t,tc))
    return

def sim(code_name,t):
    os.system('./%s  %s'%(code_name,t))
    return

def fitness(simulation,experiment):
    sim=np.loadtxt(simulation)
    exp=np.loadtxt(experiment)
    simy=sim[:,1]
    simx=sim[:,0]
    expy=exp[:,1]
    expx=exp[:,0]
    expinterp=np.interp(simx,expx,expy)
    fit=np.trapz((simy-expinterp)*(simy-expinterp),simx)/np.trapz(simy*simy,simx)
    return fit

def seq(Z):
    sarr = [str(a) for a in Z]
    return (", " . join(sarr))

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
