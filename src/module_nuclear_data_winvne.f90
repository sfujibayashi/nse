module module_nuclear_data_winvne
  implicit none
  
  ! data of
  integer :: nct_winvne
  integer :: nz_winvne
  integer :: na_winvne
  character(5),allocatable :: name_winvne(:),ref_winvne(:,:)
  logical,allocatable :: bhf_winvne(:),bex_winvne(:)
  integer,allocatable :: npt_winvne(:),nnt_winvne(:),naw_winvne(:)
  real(8),allocatable :: ams_winvne(:),spn_winvne(:),exc_winvne(:),ptf_winvne(:,:)
  integer,allocatable :: jnuc_winvne(:,:)
  integer,allocatable :: iadz_winvne(:),naz_winvne(:)
  
  !real(8) :: t9_winvne(24)
  real(8) :: t9_winvne(24) = (/0.10d0,0.15d0,0.20d0,0.30d0,0.40d0,0.50d0,0.60d0,0.70d0,0.80d0,0.90d0,1.00d0,1.50d0,2.00d0,2.50d0,3.00d0,3.50d0,4.00d0,4.50d0,5.00d0,6.00d0,7.00d0,8.00d0,9.00d0,10.0d0/)

contains
  subroutine init_winvne(fn)!, nct_in, nz_in, na_in)
    character(*),intent(in) :: fn
    ! integer,intent(in) :: nct_in, nz_in, na_in

    character(5) :: str1,str2
    integer :: k,i

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
    nct_winvne = i

    write(6,*) "# of Winvne species:", nct_winvne
    
    ! nct_winvne = nct_in
    
    allocate ( name_winvne(nct_winvne),ref_winvne(2,nct_winvne))
    allocate ( bhf_winvne(nct_winvne),bex_winvne(nct_winvne))
    allocate ( npt_winvne(nct_winvne),nnt_winvne(nct_winvne),naw_winvne(nct_winvne))
    allocate ( ams_winvne(nct_winvne),spn_winvne(nct_winvne),exc_winvne(nct_winvne),ptf_winvne(24,nct_winvne))
    
    do k=1,nct_winvne
       ! read(10,'(a5,f12.3,i4,i4,f6.1,f10.3,1x,a5)') name_winvne(k),ams_winvne(k),npt_winvne(k),nnt_winvne(k),spn_winvne(k),exc_winvne(k),ref_winvne(1,k)
       read(10,'(a5,f12.3,i4,i4,f6.1,f10.3)') name_winvne(k),ams_winvne(k),npt_winvne(k),nnt_winvne(k),spn_winvne(k),exc_winvne(k)!,ref_winvne(1,k)
       !write(6,*) k,name_winvne(k),ams_winvne(k),npt_winvne(k),nnt_winvne(k),spn_winvne(k),exc_winvne(k)!,ref_winvne(1,k)
       !ref_winvne(2,k) = ref_winvne(1,k)
       read(10,*) (ptf_winvne(i,k),i=1 ,8 )
       read(10,*) (ptf_winvne(i,k),i=9 ,16)
       read(10,*) (ptf_winvne(i,k),i=17,24)
       
       ! if(npt_winvne(k)==26 .and. nnt_winvne(k) ==30)then
       !    write(6,*) str1
       !    write(6,*) ref_winvne(1,k)
       !    write(6,*) exc_winvne(k)
       
       !    do i=1,24
       !       write(6,*) t9_winvne(i),ptf_winvne(k,i)
       !    enddo
       !    stop
       ! endif
         
       naw_winvne(k) = nint(ams_winvne(k))
    enddo
    close(10)
    
    ! nz_winvne = nz_in
    ! na_winvne = na_in
    nz_winvne = maxval(npt_winvne(:))
    na_winvne = maxval(naw_winvne(:))
    write(6,*) "Z, A max :", nz_winvne, na_winvne
    
    allocate ( jnuc_winvne(1:na_winvne,0:nz_winvne))
    allocate ( iadz_winvne(0:nz_winvne),naz_winvne(0:nz_winvne))

    jnuc_winvne(:,:)=0
    do k=1,nct_winvne
       jnuc_winvne(naw_winvne(k),npt_winvne(k)) = k
    enddo

  ! do k=1,nct_winvne
  !    if(ams_winvne(k)<60d0)then
  !       write(99,*)
  !       do i=1,24
  !          write(99,'(99es12.4)') ams_winvne(k),dble(npt_winvne(k)),spn_winvne(k),t9_winvne(i),ptf_winvne(k,i)
  !       enddo
  !    endif
  ! enddo
  ! stop

  end subroutine init_winvne

  subroutine get_ptf_winvne(t9, k, pf)
    real(8),intent(in) :: t9
    integer,intent(in) :: k
    real(8),intent(out) :: pf

    real(8) :: t4(4),pf4(4),dpf
    integer :: nt

    if(t9_winvne(2) <= t9 .and. t9 <= t9_winvne(22))then
       call locate(t9_winvne, 24,t9,nt)
    elseif(t9 < t9_winvne(2))then
       nt = 2
    elseif(t9 > t9_winvne(22))then
       nt = 22
    endif

    t4(:) = t9_winvne(nt-1:nt+2)

    pf4(:)=log10(ptf_winvne(nt-1:nt+2,k))
    
    call polint(t4,pf4,4,t9,pf,dpf)
    
    if( 3<=nt .and. nt<=21 .and. (pf > max(pf4(2),pf4(3)) .or. pf < min(pf4(2),pf4(3))) ) then
       pf=(pf4(3)-pf4(2))/(t4(3)-t4(2))*(t9-t4(2))+pf4(2)
    endif

    pf = 10d0**pf
    
  end subroutine get_ptf_winvne

  subroutine get_stat_weight_winvne(t9, k, g)
    
    real(8), intent(in)  :: t9
    integer, intent(in)  :: k
    real(8), intent(out) :: g
    
    real(8) :: pf
    
    call get_ptf_winvne(t9, k, pf)
    
    g = (2d0*spn_winvne(k) + 1d0)*pf
    
  end subroutine get_stat_weight_winvne
  
  subroutine calc_ptf_winvne(t9,g)
    real(8),intent(in) :: t9
    real(8),intent(out) :: g(nct_winvne)
    
    real(8) :: t4(4),pf4(4),pf,dpf
    integer :: nt

    integer :: k

    do k=1,nct_winvne
       g(k) = (2d0*spn_winvne(k)+1d0)
    enddo
    ! return

    if(t9>t9_winvne(24))then
       do k=1,nct_winvne
          g(k) = g(k) * ptf_winvne(24,k)
       enddo
       return
    endif
    
    if(t9_winvne(2) <= t9 .and. t9 <= t9_winvne(22))then
       call locate(t9_winvne, 24,t9,nt)
    elseif(t9 < t9_winvne(2))then
       nt = 2
    elseif(t9 > t9_winvne(22))then
       nt = 22
    endif

    t4(:) = t9_winvne(nt-1:nt+2)
    
    do k=1,nct_winvne
       ! do l=1,4
       !    pf4(l)=ptf_winvne(k,nt-2+l)
       ! enddo
       pf4(:)=log10(ptf_winvne(nt-1:nt+2,k))

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
    
  end subroutine calc_ptf_winvne
  
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
    do k=1,nct_winvne
       if(npt_winvne(k)<=nz_lim .and. nnt_winvne(k)<=nn_lim)then
          write(11,'(a5)') name_winvne(k)
          nct_reduced = nct_reduced + 1
       endif
    enddo
    write(11,*)
    write(6,*) nct_reduced

    do k=1,nct_winvne
       if(npt_winvne(k)<=nz_lim .and. nnt_winvne(k)<=nn_lim)then
          write(11,'(a5,f12.3,i4,i4,f6.1,f10.3,1x,a5)') name_winvne(k),ams_winvne(k),npt_winvne(k),nnt_winvne(k),spn_winvne(k),exc_winvne(k),ref_winvne(1,k)
          ref_winvne(2,k) = ref_winvne(1,k)
          write(11,'(8f9.2)') (ptf_winvne(k,i),i=1 ,8 )
          write(11,'(8f9.2)') (ptf_winvne(k,i),i=9 ,16)
          write(11,'(8f9.2)') (ptf_winvne(k,i),i=17,24)
       endif
    enddo

    close(11)

  end subroutine output_reduced_table

end module module_nuclear_data_winvne
