	program tscal
c-----
c	Output time domain simulation series in SAC format
c-----
c	common blocks
c-----
	parameter(NDMAX=100)
	common/propq/q0,eta,velq
		real q0,eta,velq
	common/propg/ndist, rr(NDMAX),rpow(NDMAX)
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
c-----
c	compute mean peak motion using time tomain simulation
c-----
	character mdfile*80
	integer jout
	logical dofilt
	logical dofourier
	logical doftbar
	logical dobin
	integer instyp
	real xmommg
	logical isvert
	real dist
	real fn
	real damp
	real rref, dt
	integer niter
	character str*4
c-----
c	get command line arguments
c-----
	call gcmdln('tdcal',mdfile,jout,dofilt,dofourier,
     1		doftbar,instyp,xmommg,isvert,dist,fn,damp,
     1		rref,dt,niter,kseed,str,dobin)

c-----
c	get the propogational model
c-----
	call getinfo(mdfile)
c-----
c	perform the simulation
c-----
	call tssim(niter,kseed,  dt, 
     1		dofilt, dofourier,fn,damp,fmax,instyp,
     1		doftbar,avgmx,xmommg,sigma,betas,dens,
     1		jsrc,dist,q0,eta,velq,kappa,dobin,str)
	end


	subroutine gcmdln(pname,mdfile,jout,dofilt,dofourier,
     1		doftbar,instyp,xmommg,isvert,dist,fn,damp,
     1		rref,dt,niter,kseed,str,dobin)
c-----
c	parse command line arguments to get program options
c-----
c	pname	Ch*(*)	- name fo this program
c	mdfile	Ch*(*)	- name of file with propagational model
c	jout	I	- pointer to output format
c			-OG	1 D(r) format forced to be zero at reference distacne
c					DIST D(r) DX DY
c			-OS	2 Excitation at the reference distance (log of amplitude)
c					FREQ VALUE DX DY
c					m/sec for filtered velocity
c					m/sec for Fourier acceleration spectra
c					m     for Fourier velocity     spectra
c
c			-OB	3 PERIOD PEAK_VALUE
c			-OR	4 MOMENT_MAGNITUDE PEAK_MOTION (meters)
c			-OU	5 MOMENT_MAGNITUDE DISTANCE PEAK_MOTION (meters)
c			-OT	6 TSTAR FREQ_ZERO_CROSS ETA
c			-0A	7 R MW PEAK PARAMETER
c	dofilt	L	- .true. Do Bandpass (default=.false.) 
c	dofourier L	- .true. output only the spectra  not the time domain peak 
c	doftbar	L	- .true. get average within the filtered window           
c	instyp	I	- instrument filter
c				 1 Peak ground acceleration
c				 2 Peak ground velocity
c				 3 Peak ground displacement
c				 4 Wood Anderson Standard Instrument peak displacement
c				 5 WWSSN short period
c				 9 SD peak displacment of SDOF oscillator
c				10 SV peak velocity of SDOF oscillator
c				11 SA peak acceleration of SDOF oscillator
c				12 PSV == omega sub n SD
c				13 PSA == omega * PSV
c	xmommg	R	- moment magnitude of the event
c	isvert	L	- .true. make vertical else horizontal
c	dist	R	- obtain motion at this distance
c	fn	R	- oscillator frequency for SDOF oscillator, else
c			  bandpass filter center frequency
c	damp	R	- oscillator damping (default = 0.05)
c	rref	R	- reference distance for normalization, else true unnormalized
c				peak motion is output. The is useful for -OD and -OS
c	dt	R	- sample interval for time domain synthetics - used only in
c				time domain simulation and not random vibration theory
c				default 0.01 sec
c	niter	I	- number of iterations for time domain stack (not used by rpcal)
c				(default 100)
c	kseed	I	- random number seed (default 12345)
c	str	C*4	- identification string for the output
c	dobin	L	- SAC output format (default .true.)
c				.true. = binary
c				.false. = SAC alpha
c-----
	character pname*(*)
	character mdfile*(*)
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
	integer niter, kseed
	character str*(*)
	logical dobin
c-----
c	internal variables
c-----
	integer i, nmarg
	character name*80

c-----
c	set up defaults
c-----
	mdfile = ' '
	jout = -1
	dofilt = .false.
	dofourier = .false.
	doftbar = .false.
	xmommg = 5.0
	isvert = .false.
	dist = 1.0
	fn = -1.0
	damp = 0.05
	rref = -1.0
	dt = 0.01
	niter = 100
	instyp = 1
	xmommg = 4.0
	kseed = 12345
	str = ' '
	dobin = .true.
	
	nmarg = mnmarg()
	if(nmarg .eq. 0)call usage(pname)
	i = 0
   11	continue
	i = i + 1
	if(i.gt.nmarg)go to 13
		call mgtarg(i,name)
		if(name(1:2).eq.'-A' .or. name(1:2).eq.'-a')then
			if(name(1:4).eq.'-ASC' .or. name(1:4).eq.'-asc')then
				dobin = .false.
			else
				instyp = 1
				str='AMAX'
			endif
		else if(name(1:2).eq.'-V' .or. name(1:2).eq.'-v')then
			instyp = 2
			str='VMAX'
		else if(name(1:2).eq.'-D' .or. name(1:2).eq.'-d')then
			if(name(1:3).eq.'-DT' .or. name(1:3).eq.'-dt')then
				i = i + 1
				call mgtarg(i,name)
				read(name,'(bn,f10.0)')dt
			else if(name(1:5).eq.'-DAMP'.or.name(1:5).eq.'-damp')then
				i = i + 1
				call mgtarg(i,name)
				read(name,'(bn,f10.0)')damp
			else
				instyp = 3
				str='DMAX'
			endif
C		else if(name(1:3).eq.'-WA' .or. name(1:3).eq.'-wa')then
C			instyp = 4
C		else if(name(1:6).eq.'-WWSSN' .or. name(1:6).eq.'-wwssn')then
C			instyp = 5
C		else if(name(1:5).eq.'-USGS' .or. name(1:5).eq.'-usgs')then
C			instyp = 6
C		else if(name(1:6).eq.'-ECTNO' .or. name(1:6).eq.'-ectno')then
C			instyp = 7
C		else if(name(1:6).eq.'-ECTNN' .or. name(1:6).eq.'-ectnn')then
C			instyp = 8
C		else if(name(1:5).eq.'-LRSM' .or. name(1:5).eq.'-lrsm')then
C			instyp = 14
		else if(name(1:2).eq.'-S'.or.name(1:2).eq.'-s' )then
			if(name(1:3).eq.'-SD'.or.name(1:3).eq.'-sd')then
				instyp = 9
				str='SD  '
			else if(name(1:3).eq.'-SV'.or.name(1:3).eq.'-sv')then
				instyp = 10
				str='SV  '
			else if(name(1:3).eq.'-SA'.or.name(1:3).eq.'-sa')then
				instyp = 11
				str='SA  '
			else
				i = i + 1
				call mgtarg(i,name)
				read(name,'(i10)')kseed
			endif
		else if(name(1:4).eq.'-PSV'.or.name(1:4).eq.'-psv')then
			instyp = 12
				str='PSV '
		else if(name(1:4).eq.'-PSA'.or.name(1:4).eq.'-psa')then
			instyp = 13
				str='PSA '
		else if(name(1:2).eq.'-Z'.or.name(1:2).eq.'-z')then
			isvert = .true.
		else if(name(1:3).eq.'-FN'.or.name(1:3).eq.'-fn')then
			i = i + 1
			call mgtarg(i,name)
			read(name,'(bn,f10.0)')fn
		else if(name(1:3).eq.'-MW'.or.name(1:3).eq.'-mw')then
			i = i + 1
			call mgtarg(i,name)
			read(name,'(bn,f10.0)')xmommg
		else if(name(1:3).eq.'-MO'.or.name(1:3).eq.'-mo')then
			i = i + 1
			call mgtarg(i,mdfile)
C		else if(name(1:6).eq.'-NITER'.or.name(1:2).eq.'-niter')then
C			i = i + 1
C			call mgtarg(i,name)
C			read(name,'(bn,i10)')niter
		else if(name(1:2).eq.'-R'.or.name(1:2).eq.'-r' )then
C			if(name(1:5).eq.'-RREF' .or. name(1:5).eq.'-rref')then
C				i = i + 1
C				call mgtarg(i,name)
C				read(name,'(bn,f10.0)')rref
C			else
				i = i + 1
				call mgtarg(i,name)
				read(name,'(bn,f10.0)')dist
C			endif
		else if(name(1:2).eq.'-B'.or.name(1:2).eq.'-b')then
			dofilt = .true.
		else if(name(1:2).eq.'-?'.or.name(1:2).eq.'-h')then
			call usage(pname)
C		else if(name(1:3).eq.'-OG')then
C			jout = 1
C		else if(name(1:3).eq.'-OS')then
C			jout = 2
C		else if(name(1:3).eq.'-OB')then
C			jout = 3
C		else if(name(1:3).eq.'-OR')then
C			jout = 4
C		else if(name(1:3).eq.'-OU')then
C			jout = 5
C		else if(name(1:3).eq.'-OT')then
C			jout = 6
C		else if(name(1:3).eq.'-OA')then
C			jout = 7
C		else if(name(1:4).eq.'-FAS')then
C			dofourier = .true.
C		else if(name(1:5).eq.'-FBAR')then
C			doftbar = .true.
		else if(name(1:2).eq.'-h')then
			call usage(pname)
		else if(name(1:2).eq.'-?')then
			call usage(pname)
		endif
	go to 11
   13	continue
c-----
c	error checking
c-----
	if(damp.le.0.0)damp = 0.01
	if(mdfile.eq.' ')call usage(pname)
	if(fn.lt.0.0)then
		if(instyp.ge.9 .and.instyp.le.13)call usage(pname)
	endif
	return
	end

	subroutine usage(pname)
c-----
c	output the usage
c-----
c	pname	Ch*(*)	- name of this program
c-----
	character pname*(*)

	integer LER
	parameter (LER=0)
	
	ls = lgstr(pname)
	write(LER,*)'Usage: ',pname(1:ls),' [options]'
	write(LER,*)'  Ground motion type '
	write(LER,*)
     1	'  -A     (default true ) Generate Amax'
	write(LER,*)
     1	'  -V                     Generate Vmax'
	write(LER,*)
     1	'  -D                     Generate Dmax'
	write(LER,*)
     1	'  -SD                    Generate SD  '
	write(LER,*)
     1	'  -SV                    Generate SV  '
	write(LER,*)
     1	'  -SA                    Generate SA  '
	write(LER,*)
     1	'  -PSV                   Generate PSV '
	write(LER,*)
     1	'  -PSA                   Generate PSA '
	write(LER,*)' '
	write(LER,*)
     1	'  -MW    (default 4.0)                '
	write(LER,*)
     1	'  -MODEL mfile (required) Propagation model'
	write(LER,*)
     1	'  -B     (default false) Bandpass Amax,Vmax,Dmax'
	write(LER,*)
     1	'  -H     (default 0.05)  SDOF oscillator damping'
	write(LER,*)
     1	'  -FN    (Required for -B or SDOF) filter freq'
	write(LER,*)
     1	'  -R     (default none)  Desired distance'
	write(LER,*)
     1	'  -S     (default 12345) Seed for time domain'
	write(LER,*)' '
	write(LER,*)
     1	'  -?    This help message'
	write(LER,*)
     1	'  -h    This help message'
	write(LER,*)
     1	'  -PROP Help message for the propagation model'

	stop
	end

	subroutine musage(pname)
c-----
c	output the description of the model file
c-----
c	pname	Ch*(*)	- name of this program
c-----
	character pname*(*)

	integer LER
	parameter (LER=0)
	
	ls = lgstr(pname)
	write(LER,*)'Propagation model description: ',pname(1:ls),' [options]'
	write(LER,*)' description is given at right in [] s '
	write(LER,*)' '
	write(LER,*)' The KEYWORD defines the next parameter(s). If KEYWORD'
	write(LER,*)' is not used, then the default values is taken'
	write(LER,*)
     1	'RVTDCAL1.0                   [File type description]'
	write(LER,*)
     1	'COMMENT	                     [Comment definition   ]'
	write(LER,*)
     1	'  Frankel 1996               [this is the comment  ]'
	write(LER,*)
     1	'KAPPA'
	write(LER,*)
     1	'        0.010                [default     0.0      ]'
	write(LER,*)
     1	'QETA'
	write(LER,*)
     1	'        680     0.36         [default 10000 0.0    ]'
	write(LER,*)
     1	'QVELOCITY                    [ exp(-pi f r / q vel ]'
	write(LER,*)
     1	'        3.6                  [default 3.5          ]'
	write(LER,*)
     1	'SHEAR                        [source depth Vs      ]'
	write(LER,*)
     1	'        3.6                  [default 3.5          ]'
	write(LER,*)
     1	'DENSITY                      [source density       ]'
	write(LER,*)
     1	'        2.8                  [default 2.7          ]'
	write(LER,*)
     1	'DURATION                     [distance dependent duration]'
	write(LER,*)
     1	'2                            [number of pairs      ]'
	write(LER,*)
     1	'     0.000     0.000         [distance     duration]'
	write(LER,*)
     1	'  1000.000    50.000'
	write(LER,*)
     1	'DISTANCE                     [geometrical spreading]'
	write(LER,*)
     1	'        3                    [number of segments   ]'
	write(LER,*)
     1	'     1.000    -1.000         [>=1 km, 1/R  ,then   ]'
	write(LER,*)
     1	'    70.000     0.000         [>=70 km, 1/R^0, then '
	write(LER,*)
     1	'   130.000    -0.500         [>= 130 km, 1/R^0.5   ]'
	write(LER,*)
     1	'SITE                         [Site response        ]'
	write(LER,*)
     1	'  3                          [number freq, resp pair]'
	write(LER,*)
     1	'     3.010     1.000         [freq   response]'
	write(LER,*)
     1	'     6.341     2.282'
	write(LER,*)
     1	'    82.000     2.463'
	write(LER,*)
     1	'FMAX'
	write(LER,*)
     1	'        100.0               [default 100 Hz       ]'
	write(LER,*)
     1	'SIGMA                       [constant stress drop ]'
	write(LER,*)
     1	'        150.0               [stress drop in bars  ]         '
	write(LER,*)
     1	'---------------------------------------------------'
	write(LER,*)
     1	'Other source types are AT93 for ENA, AB98CA for CA,'
	write(LER,*)
     1	'AS2000 for California'

	stop
	end
