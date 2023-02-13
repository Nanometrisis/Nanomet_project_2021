module globalpar
    use lmsvars,only:dp
    integer,parameter:: nmaxfun = 130
    integer,parameter:: igd = 1
    real(dp),allocatable:: output_data(:,:) 
    character(len=12) :: output_labels(nmaxfun)
    integer:: n_output_data,n_output_total
    character(len=100) :: output_description(nmaxfun)
    complex(dp),parameter:: czero=(0.d0,0.d0)
    complex(dp),parameter:: cone=(1.d0,0.d0)
    complex(dp),parameter:: ctwo=(2.d0,0.d0)
    complex(dp),parameter:: ci=(0.d0,1.d0)
    real(dp),parameter:: pi = 4.d0*atan(1.d0)

  contains 
    subroutine makefunlabels
      
      character(len=4) s1,s2
      character(len=10) lab
      integer :: iout

      !
      iout = 1
      output_labels(iout) = 'OMEGA'  
      output_description(iout) = 'Frequency'

      iout = iout + 1
      output_labels(iout) = 'QX'
      output_description(iout) = 'Incoming wavevector (2pi/alaty scaled),  x'
      
      iout = iout + 1
      output_labels(iout) = 'QY'
      output_description(iout) = 'Incoming wavevector (2pi/alaty scaled),  y'

      iout = iout + 1
      output_labels(iout) = 'TRANS_S'
      output_description(iout) = 'Total Transmission S-polarization'

      iout = iout + 1
      output_labels(iout) = 'TRANS_P'
      output_description(iout) = 'Total Transmission P-polarization'

      iout = iout + 1
      output_labels(iout) = 'REFLE_S'
      output_description(iout) = 'Total Reflection S-polarization'
      
      iout = iout + 1
      output_labels(iout) = 'REFLE_P'
      output_description(iout) = 'Total Reflection P-polarization'
      
      iout = iout + 1
      output_labels(iout) = 'ABSOR_S'
      output_description(iout) = 'Total Absorption S-polarization'

      iout = iout + 1
      output_labels(iout) = 'ABSOR_P'
      output_description(iout) = 'Total Absorption P-polarization'
      
     
      n_output_data= iout    
      n_output_total= iout
      write(6,*) 'Total number of output columns :', iout
      if (iout > nmaxfun) STOP ' Increase nmaxfun'
      write(6,*) ' AnisoLay 1.0, Jan 2023 '
      write(6,*) ' Output produced in file anisolayers_out.txt '
      
      
      write(21,"(100A12)") output_labels(1:n_output_data)
      
100     format(A10,'|',A80) 
      
      
    end subroutine makefunlabels
    
    integer function get_funid(label)
      implicit none
      integer:: i,nmaxfun1
      character(*),intent(in):: label

      nmaxfun1 = size(output_labels)
      get_funid = 0 
      do i=1,nmaxfun1
         if (trim(output_labels(i))==trim(label)) then
            get_funid = i
            exit
         end if
      end do
    end function get_funid
    
  end module globalpar
