#integer:: nz,nm
#real:: kb,bb,tt # k mass transfer


#def methane(): #program methane
    #real:: fr,mt,vb
    #real:: rstp,time,aa
    #real,dimension(500000):: frr,mtt,ft,mm
    #real:: rn,rp,dt,zr,bot,mn,mp,rnn,mnn
    #real:: k1,k2,k3,k4,fnn
    #integer:: n,i
    ## Reading z,t,s
    #open(20,file='woa_sel_13_decav_00.dat')
    #open(21,file='met_conc_woa_sel_13_decav_00.dat')
    #open(22,file='output_13_7750.dat')


#nz=0
#do while(.not.eof(20))
#for n in rnge(0,20)
# nz=nz+1
# Reading depth(horizont)(meters), temperature (Celsium degrees)
#read(20,*) z(nz),t(nz),s(nz)
#d(nz)=rstp(s(nz),t(nz),0.,aa)/1000. # in gr/cm^3




#Initial volume of methane ( where r=0)
'''
#mp = 2.472982421933350780e-08 #250
#mp = 1.978385937546680624e-07 #500
#mp = 6.677052539220046543e-07 #750
#mp = 1.582708750037344499e-06 #1000
#mp = 3.091228027416687975e-06 #1250
#mp = 5.341642031376037234e-06 #1500
#mp = 8.482329707231391994e-06 #1750
#mp = 1.266167000029875599e-05 #2000
#mp = 1.802804185589412789e-05 #2250
#mp = 2.472982421933350380e-05 #2500
#mp = 3.291539603593290187e-05 #2750
#mp = 4.273313625100829788e-05 #3000
#mp = 5.433142380987571505e-05 #3250
#mp = 6.785863765785113595e-05 #3500
#mp = 8.346315674025059058e-05 #3750
#mp = 1.012933600023900479e-04 #4000
#mp = 1.214976263895855245e-04 #4250
#mp = 1.442243348471530231e-04 #4500
#mp = 1.696218643204085196e-04 #4750
#mp = 1.978385937546680304e-04 #5000
#mp = 2.290229020952475991e-04 #5250
#mp = 2.633231682874632149e-04 #5500
#mp = 3.008877712766307860e-04 #5750
#mp = 3.418650900080663830e-04 #6000
#mp = 3.864035034270860494e-04 #6250
#mp = 4.346513904790057204e-04 #6500
#mp = 4.346513904790057204e-04 #6750
#mp = 5.428691012628090876e-04 #7000
#mp = 6.031356828853248709e-04 #7250
#mp = 6.677052539220047246e-04 #7500'''


'''
# Time cycle
'''
n=0; time=0.
do while(.true.)
     n=n+1
     time=time+dt
    if (n<4) :    
        # Radius calculation 
        #### With Adams–Bashforth/Moulton method 4 order
        #frr(n)=fr(zr,rp)
        zr=zr-dt*vb(rp)   
        rn=rp+dt*frr(n)
        mn=mp+4.*np.pi*rp*rp*bb*1.e-8*dt  
        rp=rn
        mm(n)=mp-mn
        mp=mn
    
     else: 
        zr=zr-dt*vb(rp)
        frr(n)=fr(zr,rp)
        rnn=rp+dt*(55.*frr(n)-59.*frr(n-1)+37.*frr(n-2)-9.*frr(n-3))/24.
        zr=zr-dt*vb(rnn)
        fnn=fr(zr,rnn)
        rn=rp+dt*(9.*fnn+19*frr(n)-5.*frr(n-1)+frr(n-2))/24.
        mn=mp+4.*pi*rp*rp*bb*1.e-8*dt
        rp=rn
        mm(n)=mp-mn
        mp=mn 
    
    if(zr<100.) :  
        break
    #write(22,'(7e13.5)') time,zr/100.,rn/1000.,mn,mm(n), kb, vb(rp)
    print ('(7e13.5)'). time,zr/100.,rn/1000.,mn,mm(n), kb, vb(rp))

'''
        



   # Selection and calculation of bubble velocity
   
