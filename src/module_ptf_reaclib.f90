module module_ptf_reaclib
  implicit none
  
  ! data for Reaclib
  integer :: nct_reaclib
  integer :: nz_reaclib
  integer :: na_reaclib
  character(5),allocatable :: name_reaclib(:),ref_reaclib(:,:)
  logical,allocatable :: bhf_reaclib(:),bex_reaclib(:)
  integer,allocatable :: npt_reaclib(:),nnt_reaclib(:),naw_reaclib(:)
  real(8),allocatable :: ams_reaclib(:),spn_reaclib(:),exc_reaclib(:),ptf_reaclib(:,:)
  integer,allocatable :: jnuc_reaclib(:,:)
  integer,allocatable :: iadz_reaclib(:),naz_reaclib(:)
  
  !real(8) :: t9_reaclib(24)
  real(8) :: t9_reaclib(24) = (/0.10d0,0.15d0,0.20d0,0.30d0,0.40d0,0.50d0,0.60d0,0.70d0,0.80d0,0.90d0,1.00d0,1.50d0,2.00d0,2.50d0,3.00d0,3.50d0,4.00d0,4.50d0,5.00d0,6.00d0,7.00d0,8.00d0,9.00d0,10.0d0/)

contains
  subroutine init_ptf_reaclib(fn)!, nct_in, nz_in, na_in)
    character(*),intent(in) :: fn
    ! integer,intent(in) :: nct_in, nz_in, na_in

    character(5) :: str1,str2
    integer :: k,i

    character(256) :: str
    
    open(10,file=fn,status="old",action="read")
    read(10,*)
    read(10,*)
    i=0
    str1=""
    loop_count:do
       str2=str1
       read(10,'(a5)') str1
       ! write(6,*) str1
       if(str1==str2) exit loop_count
       i = i + 1
    enddo loop_count
    nct_reaclib = i

    write(6,*) "# of Reaclib species:", nct_reaclib
    
    ! nct_reaclib = nct_in
    
    allocate ( name_reaclib(nct_reaclib),ref_reaclib(2,nct_reaclib))
    allocate ( bhf_reaclib(nct_reaclib),bex_reaclib(nct_reaclib))
    allocate ( npt_reaclib(nct_reaclib),nnt_reaclib(nct_reaclib),naw_reaclib(nct_reaclib))
    allocate ( ams_reaclib(nct_reaclib),spn_reaclib(nct_reaclib),exc_reaclib(nct_reaclib),ptf_reaclib(nct_reaclib,24))
    
    do k=1,nct_reaclib
       ! read(10,'(a5,f12.3,i4,i4,f6.1,f10.3,1x,a5)') name_reaclib(k),ams_reaclib(k),npt_reaclib(k),nnt_reaclib(k),spn_reaclib(k),exc_reaclib(k),ref_reaclib(1,k)
       read(10,'(a5,f12.3,i4,i4,f6.1,f10.3)') name_reaclib(k),ams_reaclib(k),npt_reaclib(k),nnt_reaclib(k),spn_reaclib(k),exc_reaclib(k)!,ref_reaclib(1,k)
       !write(6,*) k,name_reaclib(k),ams_reaclib(k),npt_reaclib(k),nnt_reaclib(k),spn_reaclib(k),exc_reaclib(k)!,ref_reaclib(1,k)
       !ref_reaclib(2,k) = ref_reaclib(1,k)
       read(10,*) (ptf_reaclib(k,i),i=1 ,8 )
       read(10,*) (ptf_reaclib(k,i),i=9 ,16)
       read(10,*) (ptf_reaclib(k,i),i=17,24)
       
       ! if(npt_reaclib(k)==26 .and. nnt_reaclib(k) ==30)then
       !    write(6,*) str1
       !    write(6,*) ref_reaclib(1,k)
       !    write(6,*) exc_reaclib(k)
       
       !    do i=1,24
       !       write(6,*) t9_reaclib(i),ptf_reaclib(k,i)
       !    enddo
       !    stop
       ! endif
         
       naw_reaclib(k) = nint(ams_reaclib(k))
    enddo
    close(10)
    
    ! nz_reaclib = nz_in
    ! na_reaclib = na_in
    nz_reaclib = maxval(npt_reaclib(:))
    na_reaclib = maxval(naw_reaclib(:))
    write(6,*) "Z, A max :", nz_reaclib, na_reaclib
    
    allocate ( jnuc_reaclib(1:na_reaclib,0:nz_reaclib))
    allocate ( iadz_reaclib(0:nz_reaclib),naz_reaclib(0:nz_reaclib))

    jnuc_reaclib(:,:)=0
    do k=1,nct_reaclib
       jnuc_reaclib(naw_reaclib(k),npt_reaclib(k)) = k
    enddo

  ! do k=1,nct_reaclib
  !    if(ams_reaclib(k)<60d0)then
  !       write(99,*)
  !       do i=1,24
  !          write(99,'(99es12.4)') ams_reaclib(k),dble(npt_reaclib(k)),spn_reaclib(k),t9_reaclib(i),ptf_reaclib(k,i)
  !       enddo
  !    endif
  ! enddo
  ! stop

  end subroutine init_ptf_reaclib

  subroutine get_ptf_reaclib(t9, k, pf)
    real(8),intent(in) :: t9
    integer,intent(in) :: k
    real(8),intent(out) :: pf

    real(8) :: t4(4),pf4(4),dpf
    real(8) :: g0, logpf
    integer :: nt

    g0 = (2d0*spn_reaclib(k)+1d0)
    
    if(t9_reaclib(2) <= t9 .and. t9 <= t9_reaclib(22))then
       call locate(t9_reaclib, 24,t9,nt)
    elseif(t9 < t9_reaclib(2))then
       nt = 2
    elseif(t9 > t9_reaclib(22))then
       nt = 22
    endif

    t4(:) = t9_reaclib(nt-1:nt+2)

    pf4(:)=log10(ptf_reaclib(k,nt-1:nt+2))
    
    call polint(t4,pf4,4,t9,pf,dpf)
    
    if( 3<=nt .and. nt<=21 .and. (pf > max(pf4(2),pf4(3)) .or. pf < min(pf4(2),pf4(3))) ) then
       logpf=(pf4(3)-pf4(2))/(t4(3)-t4(2))*(t9-t4(2))+pf4(2)
    endif

    pf = g0 * 10d0**logpf
    
  end subroutine get_ptf_reaclib

  subroutine calc_ptf_reaclib(t9,g)
    real(8),intent(in) :: t9
    real(8),intent(out) :: g(nct_reaclib)
    
    real(8) :: t4(4),pf4(4),pf,dpf
    integer :: nt

    integer :: k,l

    do k=1,nct_reaclib
       g(k) = (2d0*spn_reaclib(k)+1d0)
    enddo
    ! return

    if(t9>t9_reaclib(24))then
       do k=1,nct_reaclib
          g(k) = g(k) * ptf_reaclib(k,24)
       enddo
       return
    endif
    
    if(t9_reaclib(2) <= t9 .and. t9 <= t9_reaclib(22))then
       call locate(t9_reaclib, 24,t9,nt)
    elseif(t9 < t9_reaclib(2))then
       nt = 2
    elseif(t9 > t9_reaclib(22))then
       nt = 22
    endif

    t4(:) = t9_reaclib(nt-1:nt+2)
    
    do k=1,nct_reaclib
       ! do l=1,4
       !    pf4(l)=ptf_reaclib(k,nt-2+l)
       ! enddo
       pf4(:)=log10(ptf_reaclib(k,nt-1:nt+2))

       call polint(t4,pf4,4,t9,pf,dpf)

       if( 3<=nt .and. nt<=21 .and. (pf > max(pf4(2),pf4(3)) .or. pf < min(pf4(2),pf4(3))) ) then
          pf=(pf4(3)-pf4(2))/(t4(3)-t4(2))*(t9-t4(2))+pf4(2)
       endif

       g(k) = g(k) * 10d0**pf
       ! if(pf4(1) == 1.d0 .or. (nptf(k).ne.0.and.pf4(3).eq.1.d5))then
       !    pf=(pf4(3)-pf4(2))/(t4(3)-t4(2))*(t9-t4(2))+pf4(2)
       ! else
       !    call polint(t4,pf4,4,t9,pf,dpf)
       ! endif
    enddo
    
  end subroutine calc_ptf_reaclib
  
  subroutine output_reduced_table
    integer :: nn_lim, nz_lim
    character(100) :: fn
    integer :: k,i
    integer :: nct_reduced
   
    nn_lim = 90
    nz_lim = 55
    
    fn = "reduced"
    open(11,file=fn,status="replace",action="write")
    write(11,*)
    write(11,*) "010015020030040050060070080090100150200250300350400450500600700800900100"

    nct_reduced = 0
    do k=1,nct_reaclib
       if(npt_reaclib(k)<=nz_lim .and. nnt_reaclib(k)<=nn_lim)then
          write(11,'(a5)') name_reaclib(k)
          nct_reduced = nct_reduced + 1
       endif
    enddo
    write(11,*)
    write(6,*) nct_reduced

    do k=1,nct_reaclib
       if(npt_reaclib(k)<=nz_lim .and. nnt_reaclib(k)<=nn_lim)then
          write(11,'(a5,f12.3,i4,i4,f6.1,f10.3,1x,a5)') name_reaclib(k),ams_reaclib(k),npt_reaclib(k),nnt_reaclib(k),spn_reaclib(k),exc_reaclib(k),ref_reaclib(1,k)
          ref_reaclib(2,k) = ref_reaclib(1,k)
          write(11,'(8f9.2)') (ptf_reaclib(k,i),i=1 ,8 )
          write(11,'(8f9.2)') (ptf_reaclib(k,i),i=9 ,16)
          write(11,'(8f9.2)') (ptf_reaclib(k,i),i=17,24)
       endif
    enddo

    close(11)

  end subroutine output_reduced_table

end module module_ptf_reaclib
