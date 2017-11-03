!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM q2r
  !----------------------------------------------------------------------------
  !
  !  q2r.x:
  !     reads force constant matrices C(q) produced by the phonon code
  !     for a grid of q-points, calculates the corresponding set of
  !     interatomic force constants (IFC), C(R)
  !
  !  Input data: Namelist "input"
  !     fildyn     :  input file name (character, must be specified)
  !                   "fildyn"0 contains information on the q-point grid
  !                   "fildyn"1-N contain force constants C_n = C(q_n)
  !                   for n=1,...N, where N is the number of q-points
  !                   in the irreducible brillouin zone
  !                   Normally this should be the same as specified
  !                   on input to the phonon code
  !                   In the non collinear/spin-orbit case the files 
  !                   produced by ph.x are in .xml format. In this case 
  !                   fildyn is the same as in the phonon code + the .xml 
  !                   extension.
  !     flfrc      :  output file containing the IFC in real space
  !                   (character, must be specified)
  !     zasr       :  Indicates type of Acoustic Sum Rules used for the Born
  !                   effective charges (character):
  !                   - 'no': no Acoustic Sum Rules imposed (default)
  !                   - 'simple':  previous implementation of the asr used
  !                     (3 translational asr imposed by correction of
  !                     the diagonal elements of the force-constants matrix)
  !                   - 'crystal': 3 translational asr imposed by optimized
  !                      correction of the IFC (projection).
  !                   - 'one-dim': 3 translational asr + 1 rotational asr
  !                     imposed by optimized correction of the IFC (the
  !                     rotation axis is the direction of periodicity; it
  !                     will work only if this axis considered is one of
  !                     the cartesian axis).
  !                   - 'zero-dim': 3 translational asr + 3 rotational asr
  !                     imposed by optimized correction of the IFC.
  !                   Note that in certain cases, not all the rotational asr
  !                   can be applied (e.g. if there are only 2 atoms in a
  !                   molecule or if all the atoms are aligned, etc.).
  !                   In these cases the supplementary asr are cancelled
  !                   during the orthonormalization procedure (see below).
  !     loto_2d    :  set to .true. to activate two-dimensional treatment of LO-TO splitting. 
  !
  !  If a file "fildyn"0 is not found, the code will ignore variable "fildyn"
  !  and will try to read from the following cards the missing information
  !  on the q-point grid and file names:
  !     nr1,nr2,nr3:  dimensions of the FFT grid formed by the q-point grid
  !     nfile      :  number of files containing C(q_n), n=1,nfile
  !  followed by nfile cards:
  !     filin      :  name of file containing C(q_n)
  !  The name and order of files is not important as long as q=0 is the first
  !
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm
  USE mp_global,  ONLY : mp_startup, mp_global_end
  USE io_global, ONLY : ionode_id, ionode, stdout
  USE environment, ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  CHARACTER(len=256) :: fildyn, filin, filj, filf, flfrc
  LOGICAL           :: la2F, loto_2d
  !
  NAMELIST / input / fildyn, flfrc, prefix, zasr, la2F, loto_2d
  !
  CALL mp_startup()
  CALL environment_start('Q2R')
  !
  IF (ionode) CALL input_from_file ( )
     !
  fildyn = ' '
  flfrc = ' '
  prefix = ' '
  loto_2d=.false.
  zasr = 'no'
     !
  la2F=.false.
     !
     !
  IF (ionode)  READ ( 5, input, IOSTAT =ios )

  CALL mp_bcast(ios, ionode_id, world_comm)
  CALL errore('q2r','error reading input namelist', abs(ios))

  CALL mp_bcast(fildyn, ionode_id, world_comm)
  CALL mp_bcast(flfrc, ionode_id, world_comm)
  CALL mp_bcast(prefix, ionode_id, world_comm)
  CALL mp_bcast(zasr, ionode_id, world_comm)
  CALL mp_bcast(loto_2d, ionode_id, world_comm)
  CALL mp_bcast(la2F, ionode_id, world_comm)
  !
  CALL do_q2r(fildyn, flfrc, prefix, zasr, la2F, loto_2d)
  !
  CALL environment_end('Q2R')

  CALL mp_global_end()
  !
END PROGRAM q2r
!
!----------------------------------------------------------------------------
SUBROUTINE gammaq2r( nqtot, nat, nr1, nr2, nr3, at )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE fft_scalar, ONLY : cfft3d
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp_images, ONLY : intra_image_comm
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nqtot, nat, nr1, nr2, nr3
  REAL(DP), INTENT(IN) :: at(3,3)
  !
  INTEGER, ALLOCATABLE :: nc(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: gaminp(:,:,:,:,:), gamout(:,:,:,:,:)
  !
  REAL(DP), PARAMETER :: eps=1.D-5, eps12=1.d-12
  INTEGER  :: nsig = 10, isig, filea2F, nstar, count_q, nq, nq_log, iq, &
       icar, ipol, m1,m2,m3, m(3), nr(3), j1,j2, na1, na2, nn
  LOGICAL :: lq
  REAL(DP) :: deg, ef, dosscf
  REAL(DP) :: q(3,48), xq, resi
  character(len=256) :: name
  CHARACTER(LEN=256) :: elph_dir
  CHARACTER(LEN=6) :: int_to_char
  LOGICAL :: exst
  INTEGER :: ios

  !
  ALLOCATE (gaminp(3,3,nat,nat,48), gamout(nr1*nr2*nr3,3,3,nat,nat) )
  ALLOCATE ( nc (nr1,nr2,nr3) )
  write (stdout,*)
  write (stdout,*) '  Preparing gamma for a2F '
  write (stdout,*)
  !
  nr(1) = nr1
  nr(2) = nr2
  nr(3) = nr3
  elph_dir='elph_dir/'
!  IF (ionode) INQUIRE(FILE=TRIM(elph_dir), EXIST=exst)
!  CALL mp_bcast(exst, ionode_id, intra_image_comm)
!  IF (.NOT. exst) CALL errore('gammaq2r','elph_dir directory not exists',1)
  !
  DO isig=1, nsig
     filea2F = 50 + isig
     nc = 0
     DO count_q=1,nqtot
        name= TRIM(elph_dir) // 'a2Fq2r.' // TRIM(int_to_char(filea2F)) &
                          // '.' // TRIM(int_to_char(count_q))
        IF (ionode) open(unit=filea2F, file=name, STATUS = 'old', &
                                 FORM = 'formatted', IOSTAT=ios)
        CALL mp_bcast(ios, ionode_id, intra_image_comm)
        IF (ios /= 0) CALL errore('gammaq2r','problem opening file' &
                                   //TRIM(name), 1)
     !
     ! to pass to matdyn, for each isig, we read: degauss, Fermi energy and DOS
     !
        !
        IF (ionode) THEN
           READ(filea2F,*) deg, ef, dosscf
           READ(filea2F,*) nstar
        ENDIF
        CALL mp_bcast(deg, ionode_id, world_comm)
        CALL mp_bcast(ef, ionode_id, world_comm)
        CALL mp_bcast(dosscf, ionode_id, world_comm)
        CALL mp_bcast(nstar, ionode_id, world_comm)
        !
        CALL read_gamma ( nstar, nat, filea2F, q, gaminp )
        !
        do nq = 1,nstar
           lq = .true.
           do ipol=1,3
              xq = 0.0d0
              do icar=1,3
                 xq = xq + at(icar,ipol) * q(icar,nq) * nr(ipol)
              end do
              lq = lq .AND. (ABS(NINT(xq) - xq) < eps)
              iq = NINT(xq)
              !
              m(ipol)= mod(iq,nr(ipol)) + 1
              if (m(ipol) < 1) m(ipol) = m(ipol) + nr(ipol)
           end do !ipol
           IF (.NOT.lq) CALL errore('gammaq2r','q not allowed',1)
           !
           if(nc(m(1),m(2),m(3)) == 0) then
              nc(m(1),m(2),m(3)) = 1
              CALL TRASL( gamout, gaminp, nq, nr1, nr2, nr3, nat, m(1), m(2), m(3) )
           else
              call errore('gammaq2r',' nc already filled: wrong q grid or wrong nr',1)
           end if
        enddo ! stars for given q-point
     ENDDO ! q-points
     !
     nq_log = SUM (nc)
     if (nq_log == nr1*nr2*nr3) then
        write (stdout,*)
        write (stdout,'(" Broadening = ",F10.3)') deg
        write (stdout,'(5x,a,i4)') ' q-space grid ok, #points = ',nq_log
     else
        call errore('gammaq2r',' missing q-point(s)!',1)
     end if
     do j1=1,3
        do j2=1,3
           do na1=1,nat
              do na2=1,nat
                 call cfft3d ( gamout(:,j1,j2,na1,na2), &
                      nr1,nr2,nr3, nr1,nr2,nr3, 1, 1 )
              end do
           end do
        end do
     end do
     gamout = gamout / DBLE (nr1*nr2*nr3)
     !
     IF (ionode) close(filea2F)
     !
     filea2F = 60 + isig
     name = TRIM(elph_dir) // 'a2Fmatdyn.'// TRIM(int_to_char(filea2F)) 
     IF (ionode) THEN
     open(unit=filea2F, file=name, STATUS = 'unknown')
     !
     WRITE(filea2F,*) deg, ef, dosscf
     write(filea2F,'(3i4)') nr1, nr2, nr3

     do j1=1,3
        do j2=1,3
           do na1=1,nat
              do na2=1,nat
                 write(filea2F,'(4i4)') j1,j2,na1,na2
                 nn=0
                 DO m3=1,nr3
                    DO m2=1,nr2
                       DO m1=1,nr1
                          nn=nn+1
                          write(filea2F,'(3i4,2x,1pe18.11)')   &
                               m1,m2,m3, DBLE(gamout(nn,j1,j2,na1,na2))
                       END DO
                    END DO
                 END DO
              end do  ! na2
           end do  ! na1
        end do   !  j2
     end do   ! j1
     close(filea2F)
     ENDIF  ! ionode

     resi = SUM ( ABS ( AIMAG( gamout ) ) )

     IF (resi > eps12) THEN
        WRITE (stdout,"(/5x,' fft-check warning: sum of imaginary terms = ',es12.6)") resi
     ELSE
        WRITE (stdout,"(/5x,' fft-check success (sum of imaginary terms < 10^-12)')")
     END IF

  ENDDO
  !
  DEALLOCATE (gaminp, gamout )
  !
END SUBROUTINE gammaq2r
!
!-----------------------------------------------------------------------
subroutine read_gamma (nqs, nat, ifn, xq, gaminp)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  implicit none
  !
  ! I/O variables
  integer, intent(in) :: nqs, nat, ifn
  real(DP), intent(out) :: xq(3,48)
  complex(DP), intent(out) :: gaminp(3,3,nat,nat,48)
  !
  logical :: lrigid
  integer :: i, j, na, nb, nt, iq
  real(DP) :: phir(3),phii(3)
  CHARACTER(LEN=75) :: line
  !
  !
  Do iq=1,nqs
     IF (ionode) THEN
        READ(ifn,*)
        READ(ifn,*)
        READ(ifn,*)
        READ(ifn,'(11X,3F14.9)')  (xq(i,iq),i=1,3)
     !     write(*,*) 'xq    ',iq,(xq(i,iq),i=1,3)
        READ(ifn,*)
     END IF
     CALL mp_bcast(xq(:,iq), ionode_id, world_comm)
     do na=1,nat
        do nb=1,nat
           IF (ionode) read(ifn,*) i,j
           CALL mp_bcast(i, ionode_id, world_comm)
           CALL mp_bcast(j, ionode_id, world_comm)
           if (i.ne.na) call errore('read_gamma','wrong na read',na)
           if (j.ne.nb) call errore('read_gamma','wrong nb read',nb)
           do i=1,3
              IF (ionode) read (ifn,*) (phir(j),phii(j),j=1,3)
              CALL mp_bcast(phir, ionode_id, world_comm)
              CALL mp_bcast(phii, ionode_id, world_comm)
              do j = 1,3
                 gaminp(i,j,na,nb,iq) = CMPLX(phir(j),phii(j),kind=DP)
              end do
              !           write(*,*) 'gaminp  ',(gaminp(i,j,na,nb,iq),j=1,3)
           end do
        end do
     end do
     !
  ENDDO
  RETURN
  !
end subroutine read_gamma
!
!----------------------------------------------------------------------------
SUBROUTINE trasl( phid, phiq, nq, nr1, nr2, nr3, nat, m1, m2, m3 )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, intent(in) ::  nr1, nr2, nr3, m1, m2, m3, nat, nq
  COMPLEX(DP), intent(in) :: phiq(3,3,nat,nat,48)
  COMPLEX(DP), intent(out) :: phid(nr1,nr2,nr3,3,3,nat,nat)
  !
  INTEGER :: j1,j2,  na1, na2
  !
  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           DO na2=1,nat
              phid(m1,m2,m3,j1,j2,na1,na2) = &
                   0.5d0 * (      phiq(j1,j2,na1,na2,nq) +  &
                          CONJG(phiq(j2,j1,na2,na1,nq)))
           END DO
        END DO
     END DO
  END DO
  !
  RETURN
END SUBROUTINE trasl
!----------------------------------------------------------------------
subroutine set_zasr ( zasr, nr1,nr2,nr3, nat, ibrav, tau, zeu)
  !-----------------------------------------------------------------------
  !
  ! Impose ASR - refined version by Nicolas Mounet
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  implicit none
  character(len=10) :: zasr
  integer ibrav,nr1,nr2,nr3,nr,m,p,k,l,q,r
  integer n,i,j,n1,n2,n3,na,nb,nat,axis,i1,j1,na1
  !
  real(DP) sum, zeu(3,3,nat)
  real(DP) tau(3,nat), zeu_new(3,3,nat)
  !
  real(DP) zeu_u(6*3,3,3,nat)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  integer zeu_less(6*3),nzeu_less,izeu_less
  ! indices of vectors zeu_u that are not independent to the preceding ones,
  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  real(DP) zeu_w(3,3,nat), zeu_x(3,3,nat),scal,norm2
  ! temporary vectors and parameters

  ! Initialization.
  ! n is the number of sum rules to be considered (if zasr.ne.'simple')
  ! and 'axis' is the rotation axis in the case of a 1D system
  ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2'
  ! and (Oz) if axis='3')
  !
  if((zasr.ne.'simple').and.(zasr.ne.'crystal').and.(zasr.ne.'one-dim') &
                       .and.(zasr.ne.'zero-dim')) then
      call errore('q2r','invalid Acoustic Sum Rulei for Z*:' // zasr, 1)
  endif
  if(zasr.eq.'crystal') n=3
  if(zasr.eq.'one-dim') then
     ! the direction of periodicity is the rotation axis
     ! It will work only if the crystal axis considered is one of
     ! the cartesian axis (typically, ibrav=1, 6 or 8, or 4 along the
     ! z-direction)
     if (nr1*nr2*nr3.eq.1) axis=3
     if ((nr1.ne.1).and.(nr2*nr3.eq.1)) axis=1
     if ((nr2.ne.1).and.(nr1*nr3.eq.1)) axis=2
     if ((nr3.ne.1).and.(nr1*nr2.eq.1)) axis=3
     if (((nr1.ne.1).and.(nr2.ne.1)).or.((nr2.ne.1).and. &
          (nr3.ne.1)).or.((nr1.ne.1).and.(nr3.ne.1))) then
        call errore('q2r','too many directions of &
             &   periodicity in 1D system',axis)
     endif
     if ((ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
          ((ibrav.ne.4).or.(axis.ne.3)) ) then
        write(stdout,*) 'zasr: rotational axis may be wrong'
     endif
     write(stdout,'("zasr rotation axis in 1D system= ",I4)') axis
     n=4
  endif
  if(zasr.eq.'zero-dim') n=6

  ! Acoustic Sum Rule on effective charges
  !
  if(zasr.eq.'simple') then
     do i=1,3
        do j=1,3
           sum=0.0d0
           do na=1,nat
               sum = sum + zeu(i,j,na)
            end do
            do na=1,nat
               zeu(i,j,na) = zeu(i,j,na) - sum/nat
            end do
         end do
      end do
   else
      ! generating the vectors of the orthogonal of the subspace to project
      ! the effective charges matrix on
      !
      zeu_u(:,:,:,:)=0.0d0
      do i=1,3
         do j=1,3
            do na=1,nat
               zeu_new(i,j,na)=zeu(i,j,na)
            enddo
         enddo
      enddo
      !
      p=0
      do i=1,3
         do j=1,3
            ! These are the 3*3 vectors associated with the
            ! translational acoustic sum rules
            p=p+1
            zeu_u(p,i,j,:)=1.0d0
            !
         enddo
      enddo
      !
      if (n.eq.4) then
         do i=1,3
            ! These are the 3 vectors associated with the
            ! single rotational sum rule (1D system)
            p=p+1
            do na=1,nat
               zeu_u(p,i,MOD(axis,3)+1,na)=-tau(MOD(axis+1,3)+1,na)
               zeu_u(p,i,MOD(axis+1,3)+1,na)=tau(MOD(axis,3)+1,na)
            enddo
            !
         enddo
      endif
      !
      if (n.eq.6) then
         do i=1,3
            do j=1,3
               ! These are the 3*3 vectors associated with the
               ! three rotational sum rules (0D system - typ. molecule)
               p=p+1
               do na=1,nat
                  zeu_u(p,i,MOD(j,3)+1,na)=-tau(MOD(j+1,3)+1,na)
                  zeu_u(p,i,MOD(j+1,3)+1,na)=tau(MOD(j,3)+1,na)
               enddo
               !
            enddo
         enddo
      endif
      !
      ! Gram-Schmidt orthonormalization of the set of vectors created.
      !
      nzeu_less=0
      do k=1,p
         zeu_w(:,:,:)=zeu_u(k,:,:,:)
         zeu_x(:,:,:)=zeu_u(k,:,:,:)
         do q=1,k-1
            r=1
            do izeu_less=1,nzeu_less
               if (zeu_less(izeu_less).eq.q) r=0
            enddo
            if (r.ne.0) then
               call sp_zeu(zeu_x,zeu_u(q,:,:,:),nat,scal)
               zeu_w(:,:,:) = zeu_w(:,:,:) - scal* zeu_u(q,:,:,:)
            endif
         enddo
         call sp_zeu(zeu_w,zeu_w,nat,norm2)
         if (norm2.gt.1.0d-16) then
            zeu_u(k,:,:,:) = zeu_w(:,:,:) / DSQRT(norm2)
         else
            nzeu_less=nzeu_less+1
            zeu_less(nzeu_less)=k
         endif
      enddo
      !
      ! Projection of the effective charge "vector" on the orthogonal of the
      ! subspace of the vectors verifying the sum rules
      !
      zeu_w(:,:,:)=0.0d0
      do k=1,p
         r=1
         do izeu_less=1,nzeu_less
            if (zeu_less(izeu_less).eq.k) r=0
         enddo
         if (r.ne.0) then
            zeu_x(:,:,:)=zeu_u(k,:,:,:)
            call sp_zeu(zeu_x,zeu_new,nat,scal)
            zeu_w(:,:,:) = zeu_w(:,:,:) + scal*zeu_u(k,:,:,:)
         endif
      enddo
      !
      ! Final substraction of the former projection to the initial zeu, to get
      ! the new "projected" zeu
      !
      zeu_new(:,:,:)=zeu_new(:,:,:) - zeu_w(:,:,:)
      call sp_zeu(zeu_w,zeu_w,nat,norm2)
      write(stdout,'("Norm of the difference between old and new effective ", &
           &  "charges: " , F25.20)') SQRT(norm2)
      !
      ! Check projection
      !
      !write(6,'("Check projection of zeu")')
      !do k=1,p
      !  zeu_x(:,:,:)=zeu_u(k,:,:,:)
      !  call sp_zeu(zeu_x,zeu_new,nat,scal)
      !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," zeu_new|zeu_u(k)= ",F15.10)') k,scal
      !enddo
      !
      do i=1,3
         do j=1,3
            do na=1,nat
               zeu(i,j,na)=zeu_new(i,j,na)
            enddo
         enddo
      enddo
   endif
   !
   !
   return
 end subroutine set_zasr
!
!----------------------------------------------------------------------
subroutine sp_zeu(zeu_u,zeu_v,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two effective charges matrices zeu_u and zeu_v
  ! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
  !
  USE kinds, ONLY : DP
  implicit none
  integer i,j,na,nat
  real(DP) zeu_u(3,3,nat)
  real(DP) zeu_v(3,3,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        scal=scal+zeu_u(i,j,na)*zeu_v(i,j,na)
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp_zeu
