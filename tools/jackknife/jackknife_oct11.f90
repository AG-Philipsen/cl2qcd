!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    jackknife_nov07.f90
!
!    Program that implements jackknife statistical analysis
!
!    Usage: See below
!
!    Copyright (C) 2011 Lars Zeidlewicz
!                  2011 Christopher Pinke
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    e-mail: zeidlewicz@gmx.net
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! lz, last modified: nov 07, 2007
!
!  input:
!
!              --    files    --
!          'data':        contains a list of measurement
!
!              -- standard in --
!          statlength:    length of array for statist. analysis
!          datalength:    length of data
!                         (if datalength > statlength, the last 
!                          statlength entries of 'data' are used)
!          binnumber:     number of bins to be used   
!          decide:        if decide==0: only analysis for given
!                                       binnumber
!                         if decide==1: output 'jackbinerrs.dat'
!
!  output:
!
!             --     files    --
!          'jackbinerrs.dat':   (only if decide==1)
!                         contains a list of error, error of 
!                         susceptibility and integr. autocorrelation
!                         time as a function of binsize
!
!             -- standard out --
!          mean value, error estimate
!          susceptibility, error estimate
!          integrated autocorrelation time
!
!
!
!  remarks:
!
!          the program uses jackknife binning as described in
!          Montvay/Muenster: Quantum Fields on a Lattice (CUP)
!          Berg: Markov Chain Monte Carlo Simulations and their
!                statistical analysis (world scientific)
!
!          integrated autocorrelation time is estimated as the
!          ration of error**2 and (naive error)**2 = 
!          susceptibility/statlength  (multiplied by 0.5, this 
!          is consistent with definition in Montvay/Muenster)
!         
!          the error of the integrated autocorrelation time
!          is calculated using the chisquare distribution
!          (see Berg, ch. 2.4)
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! cp, last modified: oct 18, 2011
! 
! changes:	I included the module "getopt()" by Mark Gates (see below)
!		in order to use command line arguments for the jackknife
!		parameters instead of by-hand-input.
!
!		Optionally, you can also give the name of file containing
!		the data.
!
!		The command-line arguments can be called by -xARG or -x arg
!		
!		and one has the arguments:
!		-d	datalength
!		-s	statlength
!		-o 	optimal_binsize
!		-m	decice *optional
!		-f	filename (up to 14 characters) *optional
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CP: This is the getopt-module by Mark Gates, which I found at
! http://ews.uiuc.edu/~mrgates2/research/getopt.f90
!
! ------------------------------------------------------------
! Copyright 2008 by Mark Gates
!
! This program is free software; you can redistribute or modify it under
! the terms of the GNU general public license (GPL), version 2 or later.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! merchantability or fitness for a particular purpose.
!
! If you wish to incorporate this into non-GPL software, please contact
! me regarding licensing terms.
!
! ------------------------------------------------------------
! Fortran 95 getopt() and getopt_long(), similar to those in standard C library.
!
! ch = getopt( optstring, [longopts] )
! Returns next option character from command line arguments.
! If an option is not recognized, it returns '?'.
! If no options are left, it returns a null character, char(0).
!
! optstring contains characters that are recognized as options.
! If a character is followed by a colon, then it takes a required argument.
! For example, "x" recognizes "-x", while "x:" recognizes "-x arg" or "-xarg".
!
! optopt is set to the option character, even if it isn't recognized.
! optarg is set to the option's argument.
! optind has the index of the next argument to process. Initially optind=1.
! Errors are printed by default. Set opterr=.false. to suppress them.
!
! Grouped options are allowed, so "-abc" is the same as "-a -b -c".
!
! If longopts is present, it is an array of type(option_s), where each entry
! describes one long option.
!
!    type option_s
!        character(len=80) :: name
!        logical           :: has_arg
!        character         :: val
!    end type
!
! The name field is the option name, without the leading -- double dash.
! Set the has_arg field to true if it requires an argument, false if not.
! The val field is returned. Typically this is set to the corresponding short
! option, so short and long options can be processed together. (But there
! is no requirement that every long option has a short option, or vice-versa.)
!
! -----
! EXAMPLE
! program test
!     use getopt_m
!     implicit none
!     character:: ch
!     type(option_s):: opts(2)
!     opts(1) = option_s( "alpha", .false., 'a' )
!     opts(2) = option_s( "beta",  .true.,  'b' )
!     do
!         select case( getopt( "ab:c", opts ))
!             case( char(0))
!                 exit
!             case( 'a' )
!                 print *, 'option alpha/a'
!             case( 'b' )
!                 print *, 'option beta/b=', optarg
!             case( '?' )
!                 print *, 'unknown option ', optopt
!                 stop
!             case default
!                 print *, 'unhandled option ', optopt, ' (this is a bug)'
!         end select
!     end do
! end program test
!
! Differences from C version:
! - when options are finished, C version returns -1 instead of char(0),
!   and thus stupidly requires an int instead of a char.
! - does not support optreset
! - does not support "--" as last argument
! - if no argument, optarg is blank, not NULL
! - argc and argv are implicit
!
! Differences for long options:
! - optional argument to getopt(), rather than separate function getopt_long()
! - has_arg is logical, and does not support optional_argument
! - does not support flag field (and thus always returns val)
! - does not support longindex
! - does not support "--opt=value" syntax, only "--opt value"
! - knows the length of longopts, so does not need an empty last record

module getopt_m
	implicit none
	character(len=80):: optarg
	character:: optopt
	integer:: optind=1
	logical:: opterr=.true.
	
	type option_s
		character(len=80) :: name
		logical           :: has_arg
		character         :: val
	end type
	
	! grpind is index of next option within group; always >= 2
	integer, private:: grpind=2

contains

! ----------------------------------------
! Return str(i:j) if 1 <= i <= j <= len(str),
! else return empty string.
! This is needed because Fortran standard allows but doesn't *require* short-circuited
! logical AND and OR operators. So this sometimes fails:
!     if ( i < len(str) .and. str(i+1:i+1) == ':' ) then
! but this works:
!     if ( substr(str, i+1, i+1) == ':' ) then

character function substr( str, i, j )
	! arguments
	character(len=*), intent(in):: str
	integer, intent(in):: i, j
	
	if ( 1 <= i .and. i <= j .and. j <= len(str)) then
		substr = str(i:j)
	else
		substr = ''
	endif
end function substr


! ----------------------------------------
character function getopt( optstring, longopts )
	! arguments
	character(len=*), intent(in):: optstring
	type(option_s),   intent(in), optional:: longopts(:)
	
	! local variables
	character(len=80):: arg
	
	optarg = ''
	if ( optind > iargc()) then
		getopt = char(0)
	endif
	
	call getarg( optind, arg )
	if ( present( longopts ) .and. arg(1:2) == '--' ) then
		getopt = process_long( longopts, arg )
	elseif ( arg(1:1) == '-' ) then
		getopt = process_short( optstring, arg )
	else
		getopt = char(0)
	endif
end function getopt


! ----------------------------------------
character function process_long( longopts, arg )
	! arguments
	type(option_s),   intent(in):: longopts(:)
	character(len=*), intent(in):: arg
	
	! local variables
	integer:: i
	
	! search for matching long option
	optind = optind + 1
	do i = 1, size(longopts)
		if ( arg(3:) == longopts(i)%name ) then
			optopt = longopts(i)%val
			process_long = optopt
			if ( longopts(i)%has_arg ) then
				if ( optind <= iargc()) then
					call getarg( optind, optarg )
					optind = optind + 1
				elseif ( opterr ) then
					 print '(a,a,a)', "Error: option '", trim(arg), "' requires an argument"
				endif
			endif
			return
		endif
	end do
	! else not found
	process_long = '?'
	if ( opterr ) then
		print '(a,a,a)', "Error: unrecognized option '", trim(arg), "'"
	endif
end function process_long


! ----------------------------------------
character function process_short( optstring, arg )
	! arguments
	character(len=*), intent(in):: optstring, arg
	
	! local variables
	integer:: i, arglen
	
	arglen = len( trim( arg ))
	optopt = arg(grpind:grpind)
	process_short = optopt
	
	i = index( optstring, optopt )
	if ( i == 0 ) then
		! unrecognized option
		process_short = '?'
		if ( opterr ) then
			print '(a,a,a)', "Error: unrecognized option '-", optopt, "'"
		endif
	endif
	if ( i > 0 .and. substr( optstring, i+1, i+1 ) == ':' ) then
		! required argument
		optind = optind + 1
		if ( arglen > grpind ) then
			! -xarg, return remainder of arg
			optarg = arg(grpind+1:arglen)
		elseif ( optind <= iargc()) then
			! -x arg, return next arg
			call getarg( optind, optarg )
			optind = optind + 1
		elseif ( opterr ) then
			print '(a,a,a)', "Error: option '-", optopt, "' requires an argument"
		endif
		grpind = 2
	elseif ( arglen > grpind ) then
		! no argument (or unrecognized), go to next option in argument (-xyz)
		grpind = grpind + 1
	else
		! no argument (or unrecognized), go to next argument
		grpind = 2
		optind = optind + 1
	endif
end function process_short

end module getopt_m
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module jackglbls
implicit none
integer          :: statlength
integer          :: nbins, binsize
double precision :: average
end module

program main
use jackglbls
use getopt_m
implicit none
integer                                     :: j, optimal_binsize, decide
integer                                     :: datalength
double precision, allocatable, dimension(:) :: array, binnedarray, susbinnedarray
double precision                            :: error, suserror, saveerr, savesuserr
double precision                            :: susceptibility, tauint, savetauint, tauinterr
integer 				    :: check
integer					    :: filenamecheck
character(LEN=14)			    :: filename
external                                    :: binit, geterror, fortauerr

! CP: these are the "old" stdin-argument readins
!read(*,*) statlength
!read(*,*) datalength
!read(*,*) optimal_binsize
!read(*,*) decide

! CP: these are the "new" ones

character:: ch
type(option_s):: opts(5)

opts(1) = option_s( "statlength", .true., 's' )
opts(2) = option_s( "datalength",  .true.,  'd' )
opts(3) = option_s( "optimal_binsize",  .true.,  'o' )
opts(4) = option_s( "decide",  .true.,  'm' )
opts(5) = option_s( "file", .true., 'f' )

! CP: this is a cheap check if all options are given since I dont want to go too much into fortran
check = 0
10      format (1a)

do
	select case( getopt( "s:d:o:m:f:", opts ))
	        case( char(0))
			exit
		case( 'f' )
			read(optarg,*) filename
			print *, '# option file/f = ', filename
			check = check + 1
		case( 's' )
			read(optarg,*) statlength
			print *, '# option statlength/s = ', statlength
			check = check + 1
		case( 'd' )
			read(optarg,*) datalength
			print *, '# option datalength/d = ', datalength
			filenamecheck = 1
		case( 'o' )
			read(optarg,*) optimal_binsize
			print *, '# option optimal_binsize/o = ', optimal_binsize
			check = check + 1
		case( 'm' )
			read(optarg,*) decide
			print *, '# option decide/m = ', decide
		case( '?' )
			print *, 'unknown option ', optopt
 		stop
		case default
			print *, 'unhandled option ', optopt, ' (this is a bug)'
	end select
end do

if(check .NE. 3) then
	print *, 'not enough or too many options set. Aborting...'
	call EXIT(-1)
endif

optimal_binsize = int(statlength/optimal_binsize)

if(datalength.lt.statlength) then
  print*, 'datalength smaller than statlength'
  stop
endif

if(filenamecheck .EQ. 1) then
	open(14, file = filename)
else
	open(14, file = 'data')
endif

if(datalength.gt.statlength) then
  do j = 1, datalength-statlength
    read(14,*)
  enddo
endif

allocate(array(statlength))

do j = 1, statlength
  read(14,*) array(j)
enddo

close(14)

average        = sum(array)/dble(statlength)
susceptibility = sum(array*array)/dble(statlength) - average**2.D0

if(decide == 1) then
  open(15, file = 'jackbinerrs.dat')
  do binsize = 1, int(statlength/10)
    nbins = int(statlength/binsize)
    allocate(binnedarray(nbins))
    allocate(susbinnedarray(nbins))
    call binit(array,binnedarray,susbinnedarray)
    call geterror(binnedarray,average,error)
    call geterror(susbinnedarray,susceptibility,suserror)
    tauint = .5D0*error**2/(susceptibility/dble(statlength))
    if(binsize == optimal_binsize) then
      saveerr    = error
      savesuserr = suserror
      savetauint = tauint
    endif
    write(15,*) error, suserror, tauint
    deallocate(binnedarray)
    deallocate(susbinnedarray)
  enddo
  close(15)
else
  binsize = optimal_binsize
  nbins = int(statlength/binsize)
  allocate(binnedarray(nbins))
  allocate(susbinnedarray(nbins))
  call binit(array,binnedarray,susbinnedarray)
  call geterror(binnedarray,average,saveerr)
  call geterror(susbinnedarray,susceptibility,savesuserr)
  savetauint = .5D0*saveerr**2/(susceptibility/dble(statlength))
  deallocate(binnedarray)
  deallocate(susbinnedarray)
endif

deallocate(array)

print*, average, saveerr
print*, susceptibility, savesuserr
call fortauerr(nbins,tauinterr)
tauinterr = abs(savetauint-savetauint/tauinterr)
print*, savetauint, tauinterr

end

subroutine binit(array,binnedarray,susbinnedarray)
use jackglbls
implicit none
integer                                 :: j,k
double precision, dimension(statlength) :: array
double precision, dimension(nbins)      :: binnedarray, susbinnedarray
double precision                        :: partsum, wholesum
double precision                        :: squpartsum, squwholesum
wholesum    = sum(array)
squwholesum = sum((array - average)**2.)
do j = 1, nbins
  partsum    = 0.D0
  squpartsum = 0.D0
  do k = 1, binsize
    partsum    = partsum + array((j-1)*binsize + k)
    squpartsum = squpartsum + (array((j-1)*binsize + k)-average)**2.
  enddo
  binnedarray(j)    = (wholesum - partsum)/dble(statlength-binsize)
  susbinnedarray(j) = (squwholesum - squpartsum)/dble(statlength-binsize)
enddo
end subroutine

subroutine geterror(array, val, dval)
use jackglbls
implicit none
integer                                 :: k
double precision, dimension(statlength) :: array
double precision                        :: val, dval
dval = 0.D0
do k = 1, nbins
  dval = dval + (array(k) - val)**2.D0
enddo
dval = dval*dble(nbins-1)/dble(nbins)
dval = sqrt(dval)
end subroutine

subroutine fortauerr(N,output)
implicit none
double precision            :: output, chipdf
integer                     :: N
!double precision, parameter :: q = 0.025D0 !for p=95% confidence interval q=(1-p)/2
double precision, parameter :: q = 0.15D0 !for p=70% confidence interval q=(1-p)/2
double precision, external  :: gammaln
double precision            :: riemannsum, dx, normalization
integer                     :: j1,j2
double precision            :: xi,chi
N=N-1
normalization=exp(-gammaln(dble(N)/2.D0))
do j1=1,20000
 chi = dble(j1)/500.D0
dx = dble(N)*chi/(2.D0*10000.D0)
riemannsum=0.D0
do j2=1,10000
xi = dble(j2)*dx
riemannsum = riemannsum + dx*exp(-xi)*xi**(dble(N)/2.D0-1.D0)*normalization
enddo
if ( riemannsum.gt.dble(q) ) exit
enddo
output = chi
end subroutine

function gammaln(x)
! from NR in fortran (CUP),  ch 6.1
implicit none
double precision               :: gammaln,x
double precision               :: ser, stp, tmp, y
double precision, dimension(6) :: cof
integer                        :: j
data cof,stp/76.18009172947146D0,-86.50532032941677D0,24.01409824083091D0, &
             -1.231739572450155D0,0.1208650973866179D-2,-0.5395239384953D-5, &
	     2.5066282746310005D0/
y=x
tmp=x+5.5D0
tmp=(x+0.5D0)*log(tmp)-tmp
ser=1.000000000190015D0
do j=1,6
y=y+1.D0
ser = ser + cof(j)/y
enddo
gammaln=tmp+log(stp*ser/x)
end function
