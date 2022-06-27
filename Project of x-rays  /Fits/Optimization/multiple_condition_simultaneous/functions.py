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

def simulation7(code_name,d,fd,fbd,ftd,t,tc,rho):
    os.system('wolframscript %s %s %s %s %s %s %s %s'%(code_name,d,fd,fbd,ftd,t,tc,rho))
    return


def simulation8(code_name,d,fd,tl,tc,tcb,w,trough,troughb):
    os.system('wolframscript %s %s %s %s %s %s %s %s %s'%(code_name,d,fd,tl,tc,tcb,w,trough,troughb))
    return

def simulation10(code_name,tl,fd,tcb,rhod,ft,w,wb,trough,troughp,troughb):
    os.system('wolframscript %s %s %s %s %s %s %s %s %s %s %s'%(code_name,tl,fd,tcb,rhod,ft,w,wb,trough,troughp,troughb))
    return

def simulation11(code_name,d,fd,tl,tcb,rhod,ft,w,wb,trough,troughp,troughb):
    os.system('wolframscript %s %s %s %s %s %s %s %s %s %s %s %s'%(code_name,d,fd,tl,tcb,rhod,ft,w,wb,trough,troughp,troughb))
    return

def simulation12(code_name,d,fd,tl,tcb,rhod,ft,w,wb,trough,troughp,troughb,dthani):
    os.system('wolframscript %s %s %s %s %s %s %s %s %s %s %s %s %s'%(code_name,d,fd,tl,tcb,rhod,ft,w,wb,trough,troughp,troughb,dthani))
    return

def simulation13(code_name,d,fd,tl,tcb,fbp,w,wb,trough,troughp,troughb,bkg,rhod,yshift):
    os.system('wolframscript %s %s %s %s %s %s %s %s %s %s %s %s %s %s'%(code_name,d,fd,tl,tcb,fbp,w,wb,trough,troughp,troughb,bkg,rhod,yshift))
    return



def simulation14(code_name,d,fd,tl,tc,tcb,fbp,w,wb,trough,troughp,troughb,bkg,rhod,xshift):
    os.system('wolframscript %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s'%(code_name,d,fd,tl,tc,tcb,fbp,w,wb,trough,troughp,troughb,bkg,rhod,xshift))
    return

def simulation15(code_name,d,fd,tl,tc,tcb,fbp,w,wb,trough,troughp,troughb,bkg,rhod,xshift,yshift):
    os.system('wolframscript %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s'%(code_name,d,fd,tl,tc,tcb,fbp,w,wb,trough,troughp,troughb,bkg,rhod,xshift,yshift))
    return


def simulation16(code_name,d,fd,tl,tc,tcb,fbp,w,wb,trough,troughp,troughb,bkg,rhod,xshift,yshift,troughabo):
    os.system('wolframscript %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s'%(code_name,d,fd,tl,tc,tcb,fbp,w,wb,trough,troughp,troughb,bkg,rhod,xshift,yshift,troughabo))
    return

def simulation18(code_name,d,fd,tl,tc,tcb,fbp,w,wb,trough,troughp,troughb,bkg,rhod,xshift,yshift,troughabo,to,trougho):
    os.system('wolframscript %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s'%(code_name,d,fd,tl,tc,tcb,fbp,w,wb,trough,troughp,troughb,bkg,rhod,xshift,yshift,troughabo,to,trougho))
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


def fity(a,b):
    sim=np.loadtxt(a)
    sim_of_sorted=sim[np.argsort(sim[:,0])]
    exp=np.loadtxt(b)
    exp_of_sorted=exp[np.argsort(exp[:,0])]
    simy=sim_of_sorted[:,1]
    simx=sim_of_sorted[:,0]
    expy=exp_of_sorted[:,1]
    expx=exp_of_sorted[:,0]
    expinterp=np.interp(simx,expx,expy)
    siminterp=np.interp(expx,simx,simy)
    fit_of_sorted=np.trapz((expy-siminterp)*(expy-siminterp),expx)/np.trapz(expy*expy,expx)
    return fit_of_sorted




def findmine(expn,simn,delta):
    mat=[]
    for i in range(0,50):
        r=delta*i
        sim=np.loadtxt(simn)
        sim_of_sorted=sim[np.argsort(sim[:,0])]
        exp=np.loadtxt(expn)
        exp_of_sorted=exp[np.argsort(exp[:,0])]
        simy=sim_of_sorted[:,1]+r
        simx=sim_of_sorted[:,0]
        expy=exp_of_sorted[:,1]
        expx=exp_of_sorted[:,0]
        expinterp=np.interp(simx,expx,expy)
        siminterp=np.interp(expx,simx,simy)
        fit_of_sorted=np.trapz((expy-siminterp)*(expy-siminterp),expx)/np.trapz(expy*expy,expx)
        mat.append(fit_of_sorted)
    mindelt=delta*np.argmin(mat)+delta
    return(mindelt)

def fitylog(a,b):
    sim=np.loadtxt(a)
    sim_of_sorted=sim[np.argsort(sim[:,0])]
    exp=np.loadtxt(b)
    exp_of_sorted=exp[np.argsort(exp[:,0])]
    simy=sim_of_sorted[:,1]
    simylog=np.log10(simy)
    simx=sim_of_sorted[:,0]
    expy=exp_of_sorted[:,1]
    expylog=np.log10(expy)
    expx=exp_of_sorted[:,0]
    expinterp=np.interp(simx,expx,expylog)
    fit_of_sorted=np.trapz((simylog-expinterp)*(simylog-expinterp),simx)/np.trapz(simylog*simylog,simx)
    return fit_of_sorted

def seq(Z):
    sarr = [str(a) for a in Z]
    return (", " . join(sarr))

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
