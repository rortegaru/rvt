	subroutine brsac (IRU,LN,name,data,nerr)
c-----
c	IRU	I*4	logical unit for IO
c	LN	I*4	length of data array
c	name	C*	Name of file to be opened
c	rhdr	R*4	Real header
c	ihdr	I*4	Integer Header
c	chdr	C*	Character Header
c	data	R*4	Array of trace values
c	nerr	I*4	-1 file does not exist
c			-2 data points in file exceed dimension
c
c	NOTE IF THE BINARY FILE HAS MORE THAN LN POINTS, THEN ONLY
c	LN POINTS ARE USED
c-----
c
c  This routine reads waveform data written in SAC binary format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
c-----
	real*4 data(LN)

	logical ext

	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*(*) name
c-----
c  Read real and integer header blocks to find actual number
c  of waveform data points which is stored in ihdr(10).
c
c-----
	inquire(file=name,exist=ext)
	if(.not. ext)then
		ihdr(10) = 0
		nerr = -1
		return
	endif
		nerr = 0
		open (IRU,file=name,form='unformatted',
     &			access='direct',recl=440,status='old')
			read (IRU,rec=1) (rhdr(i), i=1,70),
     &					 (ihdr(i), i=1,40)
		close (IRU)
c-----
c
c  Read header and waveform data blocks using recored length of 158*4=632.
c
c-----
		if(ihdr(10).gt.LN)then
			maxpts = LN
			ihdr(10) = LN
			nerr = -2
		else 
			maxpts = ihdr(10)
			nerr = 0
		endif
		nrec=632+4*maxpts
		nread = 0
c-----
C-----
C	The true simple code is commented out. The code to work
C	with an early buggy SUNOS 4.1 FORTRAN follows
C-----
C		nrec=632+4*ihdr(10)
C		open (IRU,file=name,form='unformatted',
C     &			access='direct',recl=nrec)
C			read (IRU,rec=1) (rhdr(i), i=1,70),
C     &					 (ihdr(i), i=1,40),
C     &					 (chdr(i), i=1,24),
C     &					 (data(i), i=1,ihdr(10))
C		close (IRU)
C-----
		nrec=632+4*maxpts
		nread = 0
c-----
c	because of SUNOS Fortran problems with IO transfers 
c	more than 2048 bytes, read these  chunks in 
c----- 
		ndat = maxpts
		if(nrec.gt.2048)then
			open (IRU,file=name,form='unformatted',
     &				access='direct',recl=2048)
			ndat1 = (2048 - 632) / 4
			irec = 1
			read (IRU,rec=irec,err=1001) (rhdr(i), i=1,70),
     &					 (ihdr(i), i=1,40),
     &					 (chdr(i), i=1,24),
     &					 (data(i), i=1,ndat1)
			nread = nread + ndat1
 1000			continue
			nl = nread + 1
			nh = nl + 512 - 1
			if(nh.gt.ndat)then
				nh = ndat
			endif
			if(nl.gt.ndat)go to 1001
			irec = irec + 1
			read (IRU,rec=irec,err=1001) (data(i), i=nl,nh)
			nread = nread + (nh-nl+1)

			go to 1000
 1001			continue
		close (IRU)
		else
			open (IRU,file=name,form='unformatted',
     &				access='direct',recl=nrec)
			read (IRU,rec=1) (rhdr(i), i=1,70),
     &					 (ihdr(i), i=1,40),
     &					 (chdr(i), i=1,24),
     &					 (data(i), i=1,ndat)
		close (IRU)
		endif
		if(ihdr(10).gt.LN)then
			maxpts = LN
			ihdr(10) = LN
		else 
			maxpts = ihdr(10)
		endif
	return
	end

	subroutine arsac (IRU,LN,name,data,nerr)
c-----
c	IRU	I*4	logical unit for IO
c	LN	I*4	length of data array
c	name	C*	Name of file to be opened
c	rhdr	R*4	Real header
c	ihdr	I*4	Integer Header
c	chdr	C*	Character Header
c	data	R*4	Array of trace values
c	nerr	I*4	-1 file does not exist
c			-2 data points in file exceed dimension
c
c	NOTE IF THE BINARY FILE HAS MORE THAN LN POINTS, THEN ONLY
c	LN POINTS ARE USED
c-----
c
c  This routine reads waveform data written in SAC binary format.
c
c  This routine reads files written in SAC alpha (ASCII) format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
c-----
	logical ext
	real*4 data(LN)
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*(*) name
c-----
	inquire(file=name,exist=ext)
	if(.not. ext)then
		ihdr(10) = 0
		nerr = -1
		return
	endif
	nerr = 0
		open (IRU,file=name,status='old',access='sequential')
c----		rewind IRU
c  Read real header block.
c-----
		j1=1
		j2=5
		do 1110  i=1,14
			read (IRU,'(5g15.0)') (rhdr(j), j=j1,j2)
			j1=j1+5
			j2=j2+5
 1110		continue
c-----
c  Read integer header block.
c-----
		j1=1
		j2=5
		do 1120 i=1,8
			read (IRU,'(5i10)') (ihdr(j), j=j1,j2)
			j1=j1+5
			j2=j2+5
 1120		continue
c-----
c  Read character header block.
c-----
		j1=1
		j2=3
		do 1130 i=1,8
			read (IRU,'(3a8)') (chdr(j), j=j1,j2)
			j1=j1+3
			j2=j2+3
 1130		continue
		if(ihdr(10).gt.LN)then
			maxpts = LN
			ihdr(10) = LN
			nerr = -2
		else 
			maxpts = ihdr(10)
			nerr = 0
		endif
c-----
c  Read waveform data organized in a five columns block.
c-----
		nrow=(maxpts/5)+1
		do 1140 i=1,nrow
			l=i*5
			k=l-4
			read (IRU,'(5g15.0)') (data(j), j=k,l)
 1140		continue
		close (IRU)
		if(ihdr(10).gt.LN)then
			maxpts = LN
			ihdr(10) = LN
		else 
			maxpts = ihdr(10)
		endif
	return
	end

	subroutine getfhv(strcmd,fval,nerr)
c-----
c	Get float header value
c
c	strcmd	C*8	String to key on
c	val	R*4	Real value returned
c	nerr	I*4	Error condition
c				0 no error
c				1336 Value not defined
c				1337 Header variable does not exist
c-----
	character strcmd*(*)
	logical streql
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output real header
c-----
	fval = -12345.
	nerr = -1
	do 1000 i=1,70
		if(streql(strcmd,rstr(i))) then
			nerr = 0
			fval = rhdr(i)
		endif
 1000	continue
	return
	end

	subroutine getnhv(strcmd,ival,nerr)
c-----
c	Get integer header value
c
c	str	C*8	String to key on
c	ival	R*4	integer value returned
c	nerr	I*4	Error condition
c				0 no error
c				1336 Value not defined
c				1337 Header variable does not exist
c-----
	character strcmd*(*)
	logical streql
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output integer header
c-----
	ival = -12345
	nerr = -1
	do 2000 i=1,40
		if(streql(strcmd,istr(i))) then
			nerr = 0
			ival = ihdr(i)
		endif
 2000	continue
	return
	end

	subroutine getkhv(strcmd,cval,nerr)
c-----
c	Get character header value
c
c	strcmd	C*8	String to key on
c	cval	C*8	character value returned
c	nerr	I*4	Error condition
c				0  no error
c				1336 Value not defined
c				1337 Header variable does not exist
c-----
	character strcmd*(*)
	logical streql

	character cval*8
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output character header
c-----
	nerr = -1
	do 3000 i=1,24
		if(streql(strcmd,cstr(i))) then
			nerr = 0
			cval = chdr(i)
		endif
 3000	continue
	return
	end

	subroutine getlhv(strcmd,lval,nerr)
	character strcmd*(*)
	logical lval
	
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	input logical header
c-----
	call getnhv(strcmd,ival,nerr)
	if(ival.eq.0)then
		lval = .false.
	else
		lval = .true.
	endif
	return
	end

c---------------------------------------------------------
	subroutine bwsac (IWU,LN,name,data)
c---------------------------------------------------------

c
c  This routine writes out a waveform data in SAC binary format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character name*(*)
	real*4 data(LN)
c
c  The actual number of waveform data points is stored in integer
c  header 10. The file recored length is 158*4=632.
c
		nrec=632+4*ihdr(10)
		open (IWU,file=name,form='unformatted',
     &			access='direct',recl=nrec,status='unknown')
		rewind IWU
		write (IWU,rec=1) (rhdr(i),i=1,70),
     &				  (ihdr(k),k=1,40),
     &				  (chdr(j),j=1,24),
     &				  (data(l), l=1,ihdr(10))
		close (IWU)
	return
	end
c---------------------------------------------------------
	subroutine awsac (IWU,LN,name,data)
c---------------------------------------------------------
c
c  This routine writes out files in SAC alpha (ASCII) format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
	real*4 data(LN)
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character name*(*)
c
		open (IWU,file=name,status='unknown',
     &				 access='sequential')
		rewind IWU
c
c  Write real header block.
c
		j1=1
		j2=5
		do 1100 i=1,14
			write (IWU,'(5g15.7)') (rhdr(j), j=j1,j2)
			j1=j1+5
			j2=j2+5
 1100		continue
c
c  Write integer header block.
c
		j1=1
		j2=5
		do 1110 i=1,8
			write (IWU,'(5i10)') (ihdr(j), j=j1,j2)
			j1=j1+5
			j2=j2+5
 1110		continue
c
c  Write character header block.
c
		j1=1
		j2=3
		do 1120 i=1,8
			write (IWU,'(3a8)') (chdr(j), j=j1,j2)
			j1=j1+3
			j2=j2+3
 1120		continue
c
c  Ensure the last row is padded with zeros, if actual number of
c  waveform points is less than the product of number of rows
c  and number of columns which constitutes the data block.
c
		nrow=(ihdr(10)/5)+1
		nrc5=nrow*5
		if (nrc5 .gt. ihdr(10)) then
			nrcx=ihdr(10)+1
			do 1140 i=nrcx,nrc5
				data(i)=0.0
 1140			continue
		end if
c
c  Write waveform data in five columns format.
c
		do 1150 i=1,nrow
			k=i*5
			j=k-4
			write (IWU,'(5g15.7)') (data(l), l=j,k)
 1150		continue
		close (IWU)
	return
	end

	subroutine setfhv(strcmd,fval,nerr)
c-----
c	Set float header value
c
c	strcmd	C*8	String to key on
c	fval	C*8	real value set
c	nerr	I*4	Error condition
c				0  no error
c				1337 Header variable does not exist
c-----
	character strcmd*(*)
	logical streql
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output real header
c-----
	do 1000 i=1,70
		if(streql(strcmd,rstr(i))) rhdr(i) = fval
 1000	continue
	return
	end

	subroutine setnhv(strcmd,ival,nerr)
c-----
c	Set integer header value
c
c	strcmd	C*8	String to key on
c	ival	C*8	integer value set
c	nerr	I*4	Error condition
c				0  no error
c				1337 Header variable does not exist
c-----
	character strcmd*(*)
	integer ival, nerr

	logical streql
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output integer header
c-----
	do 2000 i=1,40
		if(streql(strcmd,istr(i))) ihdr(i) = ival
 2000	continue
	return
	end

	subroutine setkhv(strcmd,cval,nerr)
c-----
c	Set character header value
c
c	strcmd	C*8	String to key on
c	cval	C*8	character value set
c	nerr	I*4	Error condition
c				0  no error
c				1337 Header variable does not exist
c-----
	character strcmd*(*)
	character cval*8
	logical streql
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output character header
c-----
	do 3000 i=1,24
		if(streql(strcmd,cstr(i))) chdr(i) = cval
 3000	continue
	return
	end


	subroutine setlhv(strcmd,lval,nerr)
	character strcmd*(*)
	logical lval
	
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	character*8 rstr(70), istr(40), cstr(24)
	data (rstr(i),i=1,45)/
     1	'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1	'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1	'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1	'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1	'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1	'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1	'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1	'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'FHDR40  ', 
     1	'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
	data (rstr(i),i=46,70)/
     1	'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1	'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1	'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1	'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'FHDR65  ', 
     1	'FHDR66  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
	data istr/
     1	'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1	'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1	'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1	'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1	'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1	'IDHR11  ', 'IDHR12  ', 'IDHR13  ', 'IDHR14  ', 'IDHR15  ', 
     1	'IDHR16  ', 'IDHR17  ', 'IDHR18  ', 'IDHR19  ', 'IDHR20  ', 
     1	'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
     	data cstr/
     1	'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1	'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1	'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1	'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1	'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1	'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1	/
c-----
c	output logical header
c-----
	if(lval)then
		call setnhv(strcmd,1,nerr)
	else
		call setnhv(strcmd,0,nerr)
	endif
	return
	end


	subroutine newhdr()
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	call inihdr()
	return
	end

	subroutine inihdr()
	common/sachdr/rhdr,ihdr,chdr
	real*4 rhdr(70)
	integer*4 ihdr(40)
	character*8 chdr(24)
	do 1100 i=1,70
		rhdr(i)= -12345.0
 1100	continue
	do 1110 i=1,35
		ihdr(i)= -12345
 1110	continue
	ihdr(7)=6
	ihdr(8)=0
	ihdr(9)=0
	do 1120 i=36,40
		ihdr(i)=0
 1120	continue
	do 1130 i=1,24
		chdr(i)='-12345  '
 1130	continue
	return
	end


	logical function streql(str1,str2)
	character str1*(*), str2*(*)
	character nstr1*8, nstr2*8
c-----
c	determine if two strings are equal
c-----
	nstr1 = ' '
	l1 = lgstr(str1)
	nstr1(1:l1) = str1(1:l1)
	nstr2 = ' '
	l2 = lgstr(str2)
	nstr2(1:l2) = str2(1:l2)
c-----
	if(nstr1 .eq. nstr2)then
		streql = .true.
	else
		streql = .false.
	endif
	return 
	end

	subroutine getihv(strcmd,strval,nerr)
c-----
c	Get enumerated header value
c
c	strcmd	C*8	String to key on
c	strval	C*8	real value set
c	nerr	I*4	Error condition
c			0  no error
c			1336 Header variable undefined
c			1337 Header variable does not exist
c-----
	character strcmd*(*), strval*8
	parameter (NEVAL=50)
	character*8 eval(NEVAL)
c-----
c	header integer equivalents of enumerated values
c	e.g., IDISP == 2
c-----
      data eval/'ITIME   ','IRLIM   ','IAMPH   ','IXY     ','IUNKN   ', 
     1	'IDISP   ', 'IVEL    ', 'IACC    ', 'IB      ', 'IDAY    ', 
     2	'IO      ', 'IA      ', 'IT0     ', 'IT1     ', 'IT2     ', 
     3	'IT3     ', 'IT4     ', 'IT5     ', 'IT6     ', 'IT7     ', 
     4	'IT8     ', 'IT9     ', 'IRADNV  ', 'ITANNV  ', 'IRADEV  ', 
     5	'ITANEV  ', 'INORTH  ', 'IEAST   ', 'IHORZA  ', 'IDOWN   ', 
     6	'IUP     ', 'ILLLBB  ', 'IWWSN1  ', 'IWWSN2  ', 'IHGLP   ', 
     7	'ISRO    ', 'INUCL   ', 'IPREN   ', 'IPOSTN  ', 'IQUAKE  ', 
     8	'IPREQ   ', 'IPOSTQ  ', 'ICHEM   ', 'IOTHER  ', 'IGOOD   ', 
     9	'IGLCH   ', 'IDROP   ', 'ILOWSN  ', 'IRLDTA  ', 'IVOLTS  '/
c-----
		call getnhv(strcmd,nval,nerr)
		strval = '        '
		if(nerr.eq.0)then
			if(nval.ge.1 .and. nval.le.NEVAL)then
				strval = eval(nval)
			endif
		endif
	return
	end

	subroutine setihv(strcmd,strval,nerr)
c-----
c	Set enumerated header value
c
c	strcmd	C*8	String to key on
c	strval	C*8	real value set
c	nerr	I*4	Error condition
c			0  no error
c			1336 Header variable undefined
c			1337 Header variable does not exist
c-----
	character strcmd*(*), strval*8
	parameter (NEVAL=50)
	character*8 eval(NEVAL)
	character*8 ival(NEVAL)
	logical streql
c-----
c	header integer equivalents of enumerated values
c	e.g., IDISP == 2
c-----
      data eval/'ITIME   ','IRLIM   ','IAMPH   ','IXY     ','IUNKN   ', 
     1	'IDISP   ', 'IVEL    ', 'IACC    ', 'IB      ', 'IDAY    ', 
     2	'IO      ', 'IA      ', 'IT0     ', 'IT1     ', 'IT2     ', 
     3	'IT3     ', 'IT4     ', 'IT5     ', 'IT6     ', 'IT7     ', 
     4	'IT8     ', 'IT9     ', 'IRADNV  ', 'ITANNV  ', 'IRADEV  ', 
     5	'ITANEV  ', 'INORTH  ', 'IEAST   ', 'IHORZA  ', 'IDOWN   ', 
     6	'IUP     ', 'ILLLBB  ', 'IWWSN1  ', 'IWWSN2  ', 'IHGLP   ', 
     7	'ISRO    ', 'INUCL   ', 'IPREN   ', 'IPOSTN  ', 'IQUAKE  ', 
     8	'IPREQ   ', 'IPOSTQ  ', 'ICHEM   ', 'IOTHER  ', 'IGOOD   ', 
     9	'IGLCH   ', 'IDROP   ', 'ILOWSN  ', 'IRLDTA  ', 'IVOLTS  '/
c-----
c	equivalence of header field position and enumerated value
c	e.g., IDEP can be IUNKN, IDISP, IVEL, IVOLTS or IACC
c-----
      data ival/'IFTYPE  ','IFTYPE  ','IFTYPE  ','IFTYPE  ','IDEP    ', 
     1	'IDEP    ', 'IDEP    ', 'IDEP    ', 'IZTYPE  ', 'IZTYPE  ', 
     2	'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 
     3	'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 
     4	'IZTYPE  ', 'IZTYPE  ', 'IRADNV  ', 'ITANNV  ', 'IRADEV  ', 
     5	'ITANEV  ', 'INORTH  ', 'IEAST   ', 'IHORZA  ', 'IDOWN   ', 
     6	'IUP     ', 'ILLLBB  ', 'IWWSN1  ', 'IWWSN2  ', 'IHGLP   ', 
     7	'ISRO    ', 'IEVTYP  ', 'IEVTYP  ', 'IEVTYP  ', 'IEVTYP  ', 
     8	'IEVTYP  ', 'IEVTYP  ', 'IEVTYP  ', 'IQUAL   ', 'IQUAL   ', 
     9	'IQUAL   ', 'IQUAL   ', 'IQUAL   ', 'ISYNTH  ', 'IDEP    '/

c-----
c	now do the work, parse the table for the match
c		strcmd = ival and strval = eval, then 
c		using the table index, I, 
c			do a call setnhv(strcmd,I,nerr)
c
c	However, the IUNKN is used in both IDEP and IZTYPE 
c	and IOTHER is used in both IEVTYP and IQUAL
c	
c-----
	nerr = 0
	if(streql(strcmd,'IDEP    ')
     1			.and. streql(strval,'IUNKN   '))then
		call setnhv('IDEP    ',5,nerr)
	else if(streql(strcmd,'IZTYPE  ')
     1			.and. streql(strval,'IUNKN   '))then
		call setnhv('IZTYPE  ',5,nerr)
	else if(streql(strcmd,'IEVTYP  ')
     1			.and. streql(strval,'IOTHER  '))then
		call setnhv('IEVTYP  ',44,nerr)
	else if(streql(strcmd,'IQUAL   ')
     1			.and. streql(strval,'IOTHER  '))then
		call setnhv('IQUAL   ',44,nerr)
	else
		nerr = 1336
c-----
c	IFTYPE
c-----
		do 100 i=1,NEVAL
			if(
     1				streql(strval,eval(i)))then
				call setnhv(strcmd,i,nerr)
			endif
  100	continue
	endif
	return
	end
		
		
	
