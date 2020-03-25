	subroutine tssim(niter,kseed,  ddt, 
     1		dofilt, dofourier,fn,zeta,fmax,instyp,
     1		doftbar,avgmx,xmommg,sigma,vssrc,rhsrc,
     1		jsrc,dist,q0,eta,velq,kappa,dobin,str)
c-----
c	niter	I	- number of iterations
c	iseed	I	- random  number seed 
c	dt	R	- sampling interval
c	acc	R	- acceleration array
c	vel	R	- velocity array
c	dis	R	- displacment array
c	dofilt	
c	dofourier
c	fn
c	damp
c	fmax
c	instyp
c	doftbar
c	avgmx
c	dist	R	- distance
c-----
	integer niter, kseed
	real  ddt
	logical dofilt, dofourier, doftbar
	real fn, damp, fmax
	integer instyp
	real kappa
	logical dobin
	character str*(*)


	
	integer N2
	parameter (N2=131072)
	real acc(N2), vel(N2), dis(N2)
	real out(N2)

	logical rmv_trend
	complex  s, zamp
	real cfac
c-----
c	common to boore's routines
c		indxwind = 1  exponential
c		iran_type = 0 Gaussian
c		fcut = fmax
c-----
	include 'smsim.fi'

	iran_type = 0
c-----
c	common/magdist
c-----
	amag = xmommg
	r = dist
c-----
c	source_scale_params
c	common /source_scale_params/ stressc, dlsdm, amagc,  
c    :                       stress, fa, fb, am0, am0b_m0, fbdfa
c-----
	stressc = sigma
	dlsdm = 0.0
	fbdfa = 4.0
	amagc = 7.0
c-----
c      common /const_params/ rho, beta, prtitn, rtp, fs, pi, twopi, 
c     :                      const, iaorins, idva, freq_indep_factor
c-----
	beta = vssrc
	rho = rhsrc
	prtitn = 0.707
	rtp = 0.55
	fs = 2.0
	r_ref = 1.0
	iaorins = 1
	idva = 2
	pi = 3.1415927
	twopi = 6.2831853
c-----
c      common /site_dimin_params/ fm, akappa, dkappadmag, amagkref
c-----
	akappa = kappa
	dkappadmag = 0.0
	amagkref = 0.0
	fm = fmax
c-----
c      common /source_shape_params/ numsource, pf, pd
c-----
	numsource = jsrc
	pf = 2.0
	pd = 1.0
c-----
c	tdsim_params
c-----
	dt = ddt
	tshift = 5.0
	nsims = niter
	seed = kseed
	iseed = -abs(kseed)
	dur_fctr = 1.0
c-----
c	lowcut_filter__params
c-----
	fcut =  0.05
	norder = 2
c-----
c	common rs
c-----
	damp = zeta
c-----
c      common /q_params/ fr1, qr1, s1, ft1, ft2, fr2, qr2, s2, c_q
c-----
	fr1 = 1.0
	qr1 = q0
	s1 = eta
	ft1 = 1.0
	ft2 = 1.0
	fr2 = 1.0
	qr2 = q0
	s2 = eta
	c_q = velq
c-----
c      common /wind_params/indxwind,taper,twdtmotion,eps_wind,eta_wind
c-----
	indxwind = 1
	taper = 0.05
	twdtmotion = 1.0
	eps_wind = 0.2
	eta_wind = 0.05

	new_mr = .true.



c-----
c	If filter is used, set up the filter parameters
c-----
	if(dofilt)then
		nbut = 8
		fl = fn*1.414
		fh = fn/1.414
c-----
c		get normalization at the center frequency fn
c-----		
		fac = 6.2831853*fn*dt
c-----
c		     i omega dt
c		z = e
c-----
		s = cmplx(cos(fac),sin(fac))
		zamp = cmplx(1.0,0.0)
		call zpassf(nbut,2,.false.,fl,fh,s,zamp,dt,1)
		cfac = cabs(zamp)
	endif

	rmv_trend = .false.

	NSIMS=1
	LS=lgstr(str)

	do 1000 isum = 1, nsims
c-----
c		create the acceleration
c-----
		call getacc(acc)  
		v0 = 0.0
		d0 = 0.0
c-----
c		get ground velocity, displacement
c-----
		 call acc2vd(acc, npts, dt, rmv_trend, v0, d0, vel, dis)  
c-----
c		get pga, pgv, pgd
c-----
c-----
c		compute specific motions
c-----
		if(instyp.eq.1)then
			if ( dofilt)then
c-----
c			get filtered acceleration
c-----
			call zpass(acc,nbut,2,.false.,fl,fh,dt,npts)
			call doout(acc,npts,dt,str,dofilt,dobin)
			else
			call doout(acc,npts,dt,str,dofilt,dobin)
			endif
		else if(instyp.eq.2)then
			if( dofilt)then
c-----
c			get filtered velocity
c-----
			call zpass(vel,nbut,2,.false.,fl,fh,dt,npts)
			call doout(vel,npts,dt,str,dofilt,dobin)
			else
			call doout(vel,npts,dt,str,dofilt,dobin)
			endif
		else if(instyp.eq.3)then
			if(dofilt)then
c-----
c			get filtered displacement
c-----
			call zpass(dis,nbut,2,.false.,fl,fh,dt,npts)
			call doout(dis,npts,dt,str,dofilt,dobin)
			else
			call doout(dis,npts,dt,str,dofilt,dobin)
			endif
		else if(instyp.ge.9 .and.instyp.le.13)then
c-----
c			filter the acceleration and select 
c			sd, sv, sa,  psv or psa
c-----
			call psrv(fn,damp,psv,acc,out,sd,sv,sa,dt,npts) 
			call doout(out,npts,dt,str,dofilt,dobin)
C			if(instyp.eq.9)then
C				vmax = sd
C			else if(instyp.eq.10)then
C				vmax = sv
C			else if(instyp.eq.11)then
C				vmax = sa
C			else if(instyp.eq.12)then
C				vmax = psv
C			else if(instyp.eq.13)then
C				vmax = 6.2831853*fn*psv
C			endif
		endif
 1000	continue
c-----
c	compute the averages
c-----
	
	return
	end


	subroutine domxmn(x,npts,depmax,depmin,depmen)
c-----
c	get extremal values of the time series
c-----
	real*4 x(*)
	real*4 depmax,depmin,depmen
	integer*4 npts
	depmax = -1.0e+38
	depmin =  1.0e+38
	sum = 0.0
	do 1000 i=1, npts
		if( x(i) .gt. depmax) depmax = x(i)
		if( x(i) .lt. depmin) depmin = x(i)
		sum = sum + x(i)
 1000	continue
	if(npts.gt.0)then
		depmen = sum / npts
	else
		depmen = -12345.
	endif
	return
	end

	subroutine doout(out,npts,dt,str,dofilt,dobin)
	real out(npts), dt
	integer npts
	character str*(*)
	logical dofilt, dobin

	character stnm*8
	character fname*16
	real depmax, depmin, depmen


	call domxmn(out,NPTS,depmax,depmin,depmen)
	call newhdr()
	call setnhv('NPTS', NPTS, nerr)
	call setfhv('DELTA', dt, nerr)
	call setfhv('DEPMIN  ',depmin,nerr)
	call setfhv('DEPMAX  ',depmax,nerr)
	call setfhv('DEPMEN  ',depmen,nerr)
	call setfhv('B', 0.0, nerr)
	call setfhv('E', ( NPTS-1)*dt, nerr)
	call setlhv('LEVEN   ',.true.,nerr)
	call setlhv('LSPOL   ',.true.,nerr)
	call setlhv('LOVROK  ',.true.,nerr)
	call setlhv('LCALDA  ',.false.,nerr)
	call setihv('IFTYPE   ','ITIME   ',nerr)
	call setihv('IZTYPE   ','IB      ',nerr)
	stnm =' '
	stnm(1:4)=str(1:4)
	call setkhv('KSTNM', stnm, nerr)
	ls = lgstr(str)
	fname=str(1:ls)//'.sac'
	if(dobin)then
		call bwsac(1,NPTS,fname,out)
	else
		call awsac(1,NPTS,fname,out)
	endif
	return
	end
