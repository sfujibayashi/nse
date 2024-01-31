subroutine index_rank(n_rank, n_arr, array_in, index_out, fac)
  implicit none
  integer,intent(in) :: n_rank, n_arr
  real(8),intent(in) :: array_in(n_arr)
  integer,intent(out) :: index_out(n_rank)
  real(8),intent(in) :: fac

  real(8) :: array(n_arr)
  
  integer :: i, irank, irank_prim, ind
  integer :: index_buf(n_rank)
  real(8) :: rank_val(n_rank)
  
  array(:) = array_in(:)*fac

  rank_val(:) = minval(array(:))
  index_buf(:) = 0

  do i=1,n_arr

     loop_rank:do irank=1,n_rank
        if( array(i) >= rank_val(irank) )then
           do irank_prim=n_rank,irank+1,-1
              rank_val(irank_prim) = rank_val(irank_prim-1)
              index_buf(irank_prim) = index_buf(irank_prim-1) 
           enddo
           rank_val(irank) = array(i)
           index_buf(irank) = i
           exit loop_rank
        endif
     enddo loop_rank
     
     ! write(6,'(i8,es12.4)') i, array_in(i)
     ! write(6,'(99i12)')    index_buf(:)
     ! write(6,'(99es12.4)') rank_val(:)
     ! write(6,'(99es12.4)') (array_in(index_buf(irank)),irank=1,min(n_rank,i))
     
  enddo

  index_out(:) = index_buf(:)

end subroutine index_rank
