import numpy as np 
import pandas as pd 
g = 981          #gravity constant, 
rgaz = 82.0575   #individual gas constant
rc = 0.0584      # critical radius
dif = 1.5e-05    #diffusion
hg = 7.9e05      #Henry law constant
pa = 1.



def dvisc_f(s,t,p):
    def mu0(t):
        mu0=0.001002*10.**((22.6872-t*(1.0978+0.001827*t))/(t+89.93))
    return mu0
    
    c=(rstp(s,t,p,a1)*s)/1806.55
    dvisc=(1.+(1.0675e-4+t*5.185e-5)*np.sqrt(c)+(2.591e-3+t*3.3e-5)*c)*mu0(t)
    pa=p+10.13
    dp=pa*(-1.8266e-8+t*(1.3817e-9-2.6362e-11*t)+pa*(9.8965e-13+t*(-6.325e-14+1.2115e-15*t)))
    dvisc=dvisc+dp
    print (dvisc)
    return dvisc

def rstp(s,t,p,aa):
   if p == 0.:
       rst0 =(-0.157406+t*(6.793952e-2+t*(-9.095290e-3+                          
       t*(1.001685e-4+t*(-1.120083e-6+t*6.536332e-9)))) +
       s*(8.24493e-1+t*(-4.0899e-3+t*(7.6438e-5+                    
       t*(-8.2467e-7+t*5.3875e-9)))+(-5.72466e-3+                             
       t*(1.0227e-4-1.6546e-6*t))*np.sqrt(abs(s))+4.8314e-4*s)+1000.)
   return rst0  

# Selection and calculation of bubble velocity
def vb(r,tt): # bubble
    #real:: r
    #real:: m1,m2,j1,j2,vbmin 
    if (r<4.e03): 
        m1 = -0.849
        m2 = -0.815
        j1 = 0.733
        j2 = 4.792e-4
        vbmin = 22.16
    elif (r >= 4.e03 and r<4.e04 ): 
        m1 = 0.0
        m2 = 0.0
        j1 = 11.05
        j2 = 0.0
        vbmin = 19.15
    else :
        print  ('Incorrect radius')
        #break
    
    vb=(vbmin+j1*(r-rc)**m1)*np.exp(j2*tt*(r-rc)**m2) # cm/sec
    return vb

def fr(zr,r):    
    # Interpolation t,s to zr
    for i in range(1,nz-1):
        if(zr>=z[i] and zr<z[i+1]) :
            aa=(zr-z[i])/(z[i+1]-z[i])
            tt=t[i]+(t[i+1]-t[i])*aa
            ss=s[i]+(s[i+1]-s[i])*aa
            dd=d[i]+(d[i+1]-d[i])*aa
            break 
        
    for i in range(1,nm-1):
        if(zr>=zm[i] and zr < zm[i+1]): 
            cc=cm[i]+(cm[i+1]-cm[i])*(zr-zm[i])/(zm[i+1]-zm[i])
            break # exit

    # Velocity of bubble rising
    vbb=vb(r,tt)
    dvisc = dvisc_f(ss,tt,0.)
    print (dvisc)
    v= dvisc/((dd*1000.)*1.e04 )# cm^2/sec
    rcm=r*1.e-04
    re=2.*rcm*vbb/v # Re, d/less
    kb=sqrt(2/pi*(1-2.89*re**-0.5)*dif*vbb/rcm) # cm/sec
    dm=0. # mean density
    
    for i in range( 1,nz) : #do i=1,nz
        if (zr>z(i)) :
            dm=dm+(d(i+1)+d(i))*0.5*(z(i+1)-z(i))
        else :
            break #exit

    dm=dm/zr
    sur=surften(ss,tt)/(r*1.e-06)
    pb=pa+(0.1*g*dm*zr+2.*sur)*9.87167e-6
    a1=3.*(pa*101325.+dm*1000.*g/100.*zr/100.)+4.*sur
    a2=r*1.e-06*dm*1000.*g/100.*vbb/100.
    bb=kb*(cc-pb/hg)
    a3=3*rgaz*(tt+273.15)*bb/100.*101325.
    fr=(a2+a3)/a1*1.e6
    return (fr, aa,tt,ss,dd)

#
def mt(r):
    #use methan
    #implicit none
    #real:: r,rcm
    rcm=r*1.e-04
    mt=bb*4.*3.14256*rcm*rcm
    return 
    #end function mt 


df1 = pd.read_csv(r'C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re\woa_sel_13_decav_00.dat',
            delimiter = '\t',names = ['depth','temp','sal'])

z = df1.depth * 100 # in centimeter 
nz = len(z)
t = df1.temp  # Celsium
s = df1.sal    # 
aa = 0 
d = rstp(s,t,0.,aa)/1000. # in gr/cm^3

# Counter for concentration of methane 
df2 = pd.read_csv(
    r'C:\Users\elp\OneDrive - NIVA\Documents\Projects\PERMAFLUX.TRK\Bubbles\re\met_conc_woa_sel_13_decav_00.dat',
    header = None,delimiter = ' ',names = ['depth','conc'])
zm = df2.depth*100  #centimeter 
#print (df2)
nm = len(zm)
cm = df2.conc/(1.e9) # mole 

# Initial parameters 

rp=7750. # r bubble in [micrometer]
zr = bot =0.79e04 # depth in [cm]
dt=0.1  # time step [second]
mp = 7.367261933181644754e-04 #7750 volume of methane ( where r=0)

#def time_cycle():
### Time cycle 
n=0.
time=0.
while True: 
    n=n+1
    time=time+dt
    if (n<4) :    
        # Radius calculation 
        # Adams-Bashforth/Mouldton method 4 order 
        frr[n],aa,tt,ss,dd = fr(zr,rp)
        vb = vb(rp) # get bubble velocity 
        zr=zr-dt*vb    
        rn=rp+dt*frr(n)
        mn=mp+4.*np.pi*rp*rp*bb*1.e-8*dt  
        rp=rn
        mm[n]=mp-mn
        mp=mn
    
    else: 
        zr=zr-dt*vb(rp)
        frr[n]=fr(zr,rp)
        rnn=rp+dt*(55.*frr(n)-59.*frr(n-1)+37.*frr(n-2)-9.*frr(n-3))/24.
        zr=zr-dt*vb(rnn)
        fnn=fr(zr,rnn)
        rn=rp+dt*(9.*fnn+19*frr(n)-5.*frr(n-1)+frr(n-2))/24.
        mn=mp+4.*pi*rp*rp*bb*1.e-8*dt
        rp=rn
        mm[n]=mp-mn
        mp=mn 
    
    if(zr<100.) :  
        break
  
  
    
#write(22,'(7e13.5)') time,zr/100.,rn/1000.,mn,mm(n), kb, vb(rp)
    print ( '(7e13.5)', time,zr/100.,rn/1000.,mn,mm(n), kb, vb(rp) )
#print (df1.head())
#print (df2.head())
