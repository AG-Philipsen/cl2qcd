!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    jackknife_nov07.f90
!
!    Program that implements jackknife statistical analysis
!
!    Usage: See below
!
!    Copyright (C) 2011 Lars Zeidlewicz
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

module jackglbls
implicit none
integer          :: statlength
integer          :: nbins, binsize
double precision :: average
end module

program main
use jackglbls
implicit none
integer                                     :: j, optimal_binsize, decide
integer                                     :: datalength
double precision, allocatable, dimension(:) :: array, binnedarray, susbinnedarray
double precision                            :: error, suserror, saveerr, savesuserr
double precision                            :: susceptibility, tauint, savetauint, tauinterr
external                                    :: binit, geterror, fortauerr

read(*,*) statlength
read(*,*) datalength
read(*,*) optimal_binsize
read(*,*) decide

optimal_binsize = int(statlength/optimal_binsize)

if(datalength.lt.statlength) then
  print*, 'datalength smaller than statlength'
  stop
endif

open(14, file = 'data')

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
