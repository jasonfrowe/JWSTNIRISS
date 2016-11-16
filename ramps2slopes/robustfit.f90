subroutine robustfit(naxes,Imagecube,bpix,sat,zpt,slope,Image)
! naxes(1,2,3) x,y,nframes
use precision
implicit none
! import vars
integer, dimension(3) :: naxes
real(double) :: bpix,sat
real(double), dimension(:,:) :: zpt,slope,Image
real(double), dimension(:,:,:) :: Imagecube
! local vars
integer :: i,j,k,npt
real(double) :: a,b,dnaxes3,abdev
real(double), allocatable, dimension(:) :: x,y

allocate(x(naxes(1)),y(naxes(1)))
dnaxes3=dble(naxes(1)) !pre-compute double 

write(6,*) size(Imagecube,1),size(Imagecube,2),size(Imagecube,3)

!loop over 2D image to measure ramps
do i=1,naxes(2)
   do j=1,naxes(3)
      npt=0 !initialize number of points to fit to zero.
      do k=1,naxes(1) !collect all values for a single pixel
         if((Imagecube(k,i,j).lt.sat).and.(Imagecube(k,i,j).gt.bpix))then
            npt=npt+1
            x(npt)=dble(npt)
            y(npt)=Imagecube(k,i,j)
         endif
      enddo

!     robust fitting
      if(npt.ge.2)then
         call medfit(x,y,npt,a,b,abdev)
!        create Image
         zpt(i,j)=a
         slope(i,j)=b
         Image(i,j)=b*dnaxes3 !construct image
      else
         zpt(i,j)=Imagecube(1,i,j)
         slope(i,j)=0.0d0
         Image(i,j)=Imagecube(naxes(3),i,j) !construct image
      endif

   enddo
enddo


return
end