	program dorvmat
c-----
c	Output random vibration theory prediction of peak motion
c-----
c	common blocks
	parameter(NDMAX=100)
c-----
c       integer numdd
c       real ddt(NDMAX)
 	common/propq/q0,eta,velq
		real q0,eta,velq
	common/propg/ndist,rr(NDMAX), rpow(NDMAX)
		integer ndist
		real rr, rpow
	common/propd/ntime, rt(NDMAX),  tt(NDMAX)
		integer ntime	
		real rt, tt
	common/props/nsite, fs(NDMAX),  ss(NDMAX), kappa
		integer nsite
		real fs, ss, kappa
	common/srcfm/fmax, jsrc, sigma, betas, dens
		real fmax, sigma, betas, dens
		integer jsrc
	common/proc/comment
		character comment*80
	character mdfile*80
	integer jout
	logical dofilt
	logical dofourier
	logical doftbar
	integer instyp
	real xmommg
	logical isvert
	real dist
	real fn
	real damp
	real rref, dt
	integer niter
	character str*4
        real amaxrf,amaxr,amaxr05,amaxr95
c-----
c	get command line arguments
c-----
c	call gcmdln('rvcal',mdfile,jout,dofilt,dofourier,
c	call comvalues('rvcal',mdfile,jout,dofilt,dofourier,
c    1		doftbar,instyp,xmommg,isvert,dist,fn,damp,
c    1		rref,dt,niter,kseed,str)

c-----
c	get the propogational model
c-----
c      call getinfo(mdfile)
c-----
c	perform the simulation
c-----
        call initp()
c       write(*,*)'Qo'
        read(*,*)q0
c       write(*,*)'Eta'
        read(*,*)eta
c       write(*,*)'velq'
        read(*,*)velq
c       write(*,*)'ndist'
        read(*,*)ndist
        do 41 i=1,ndist
        read(*,*)rr(i),rpow(i) 
   41   continue
c       write(*,*)'ntime'
        read(*,*)ntime
        do 42 i=1,ntime
        read(*,*)rt(i),tt(i) 
   42   continue
c       write(*,*)'nsite'
        read(*,*)nsite
        do 43 i=1,nsite
        read(*,*)fs(i),ss(i) 
   43   continue
c       write(*,*)'kappa'
        read(*,*)kappa
c       write(*,*)'fmax'
        read(*,*)fmax
c       write(*,*)'jsrc'
        read(*,*)jsrc
c       write(*,*)'sigma'
        read(*,*)sigma
c       write(*,*)'betas'
        read(*,*)betas
c       write(*,*)'dens'
        read(*,*)dens
c       fn=1.1
c       dist=420.0
c       rref=200.0
c       xmommg=4.0
c       jout=1
c       str='VMAX'
c       damp=0.05
c       dofilt=.true.
c	common/propq/q0,eta,velq
ccommon/propg/ndist,rr(NDMAX), rpow(NDMAX)
ccommon/propd/ntime, rt(NDMAX),  tt(NDMAX)
ccommon/props/nsite, fs(NDMAX),  ss(NDMAX), kappa
ccommon/srcfm/fmax, jsrc, sigma, betas, dens
ccommon/proc/comment
	call rvmat(amaxrf,amaxr,amaxr05,amaxr95)
        write(*,*) amaxrf,amaxr,amaxr05,amaxr95
c	write(*,*)fn,dist,rref,amaxrf,amaxr,amaxr05,amaxr95,
c    1          xmommg,jout,str,damp,dofilt 
c	call output(fn,dist,rref,amaxrf,amaxr,amaxr05,amaxr95,
c    1          xmommg,jout,str,damp,dofilt)
	end
