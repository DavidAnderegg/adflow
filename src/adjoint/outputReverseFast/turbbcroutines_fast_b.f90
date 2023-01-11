!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
module turbbcroutines_fast_b
  implicit none

contains
!      ==================================================================
  subroutine applyallturbbcthisblock(secondhalo)
!
!       applyallturbbcthisblock sets the halo values of the
!       turbulent variables and eddy viscosity for the block the
!       variables in blockpointers currently point to.
!
    use constants
    use blockpointers
    use flowvarrefstate
    use inputphysics
    implicit none
!
!      subroutine arguments.
!
    logical, intent(in) :: secondhalo
!
!      local variables.
!
    integer(kind=inttype) :: nn, i, j, l, m
    real(kind=realtype), dimension(:, :, :, :), pointer :: bmt
    real(kind=realtype), dimension(:, :, :), pointer :: bvt, ww1, ww2
! loop over the boundary condition subfaces of this block.
bocos:do nn=1,nbocos
! loop over the faces and set the state in
! the turbulent halo cells.
      if (.not.wallfunctions) then
        select case  (bcfaceid(nn)) 
        case (imin) 
          do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
            do i=bcdata(nn)%icbeg,bcdata(nn)%icend
              do l=nt1,nt2
                w(1, i, j, l) = bvti1(i, j, l)
                do m=nt1,nt2
                  w(1, i, j, l) = w(1, i, j, l) - bmti1(i, j, l, m)*w(2&
&                   , i, j, m)
                end do
              end do
            end do
          end do
        case (imax) 
          do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
            do i=bcdata(nn)%icbeg,bcdata(nn)%icend
              do l=nt1,nt2
                w(ie, i, j, l) = bvti2(i, j, l)
                do m=nt1,nt2
                  w(ie, i, j, l) = w(ie, i, j, l) - bmti2(i, j, l, m)*w(&
&                   il, i, j, m)
                end do
              end do
            end do
          end do
        case (jmin) 
          do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
            do i=bcdata(nn)%icbeg,bcdata(nn)%icend
              do l=nt1,nt2
                w(i, 1, j, l) = bvtj1(i, j, l)
                do m=nt1,nt2
                  w(i, 1, j, l) = w(i, 1, j, l) - bmtj1(i, j, l, m)*w(i&
&                   , 2, j, m)
                end do
              end do
            end do
          end do
        case (jmax) 
          do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
            do i=bcdata(nn)%icbeg,bcdata(nn)%icend
              do l=nt1,nt2
                w(i, je, j, l) = bvtj2(i, j, l)
                do m=nt1,nt2
                  w(i, je, j, l) = w(i, je, j, l) - bmtj2(i, j, l, m)*w(&
&                   i, jl, j, m)
                end do
              end do
            end do
          end do
        case (kmin) 
          do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
            do i=bcdata(nn)%icbeg,bcdata(nn)%icend
              do l=nt1,nt2
                w(i, j, 1, l) = bvtk1(i, j, l)
                do m=nt1,nt2
                  w(i, j, 1, l) = w(i, j, 1, l) - bmtk1(i, j, l, m)*w(i&
&                   , j, 2, m)
                end do
              end do
            end do
          end do
        case (kmax) 
          do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
            do i=bcdata(nn)%icbeg,bcdata(nn)%icend
              do l=nt1,nt2
                w(i, j, ke, l) = bvtk2(i, j, l)
                do m=nt1,nt2
                  w(i, j, ke, l) = w(i, j, ke, l) - bmtk2(i, j, l, m)*w(&
&                   i, j, kl, m)
                end do
              end do
            end do
          end do
        end select
      end if
! set the value of the eddy viscosity, depending on the type of
! boundary condition. only if the turbulence model is an eddy
! viscosity model of course.
      if (eddymodel) then
        if (bctype(nn) .eq. nswalladiabatic .or. bctype(nn) .eq. &
&           nswallisothermal) then
! viscous wall boundary condition. eddy viscosity is
! zero at the wall.
          call bceddywall(nn)
        else
! any boundary condition but viscous wall. a homogeneous
! neumann condition is applied to the eddy viscosity.
          call bceddynowall(nn)
        end if
      end if
! extrapolate the turbulent variables in case a second halo
! is needed.
      if (secondhalo) call turb2ndhalo(nn)
    end do bocos
  end subroutine applyallturbbcthisblock
  subroutine bceddynowall(nn)
!
!       bceddynowall sets the eddy viscosity in the halo cells of
!       subface nn of the block given in blockpointers. the boundary
!       condition on the subface can be anything but a viscous wall.
!       a homogeneous neumann condition is applied, which means that
!       the eddy viscosity is simply copied from the interior cell.
!
    use constants
    use blockpointers
    implicit none
!
!      subroutine arguments.
!
    integer(kind=inttype), intent(in) :: nn
!
!      local variables.
!
    integer(kind=inttype) :: i, j
! determine the face id on which the subface and copy
    select case  (bcfaceid(nn)) 
    case (imin) 
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          rev(1, i, j) = rev(2, i, j)
        end do
      end do
    case (imax) 
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          rev(ie, i, j) = rev(il, i, j)
        end do
      end do
    case (jmin) 
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          rev(i, 1, j) = rev(i, 2, j)
        end do
      end do
    case (jmax) 
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          rev(i, je, j) = rev(i, jl, j)
        end do
      end do
    case (kmin) 
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          rev(i, j, 1) = rev(i, j, 2)
        end do
      end do
    case (kmax) 
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          rev(i, j, ke) = rev(i, j, kl)
        end do
      end do
    end select
  end subroutine bceddynowall
  subroutine bceddywall(nn)
!
!       bceddywall sets the eddy viscosity in the halo cells of
!       viscous subface nn of the block given in blockpointers.
!       as the eddy viscosity is zero at the wall, the value in the
!       halo is simply the negative value of the first interior cell.
!
    use constants
    use blockpointers
    implicit none
!
!      subroutine arguments.
!
    integer(kind=inttype), intent(in) :: nn
!
!      local variables.
!
    integer(kind=inttype) :: i, j
    real(kind=realtype) :: result1
! determine the face id on which the subface is located and
! loop over the faces of the subface and set the eddy viscosity
! in the halo cells.
    select case  (bcfaceid(nn)) 
    case (imin) 
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          result1 = saroughfact(2, i, j)
          rev(1, i, j) = result1*rev(2, i, j)
        end do
      end do
    case (imax) 
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          result1 = saroughfact(il, i, j)
          rev(ie, i, j) = result1*rev(il, i, j)
        end do
      end do
    case (jmin) 
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          result1 = saroughfact(i, 2, j)
          rev(i, 1, j) = result1*rev(i, 2, j)
        end do
      end do
    case (jmax) 
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          result1 = saroughfact(i, jl, j)
          rev(i, je, j) = result1*rev(i, jl, j)
        end do
      end do
    case (kmin) 
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          result1 = saroughfact(i, j, 2)
          rev(i, j, 1) = result1*rev(i, j, 2)
        end do
      end do
    case (kmax) 
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          result1 = saroughfact(i, j, kl)
          rev(i, j, ke) = result1*rev(i, j, kl)
        end do
      end do
    end select
  end subroutine bceddywall
  subroutine bcturbinflow(nn)
!
!       bcturbinflow applies the implicit treatment of the inflow
!       boundary conditions to subface nn. as the inflow boundary
!       condition is independent of the turbulence model, this routine
!       is valid for all models. it is assumed that the pointers in
!       blockpointers are already set to the correct block on the
!       correct grid level.
!
    use constants
    use blockpointers
    use flowvarrefstate
    implicit none
!
!      subroutine arguments.
!
    integer(kind=inttype), intent(in) :: nn
!
!      local variables.
!
    integer(kind=inttype) :: i, j, l
! loop over the faces of the subfaces and set the values of
! bvt and bmt such that the inflow state is linearly extrapolated
! with a fixed state at the face.
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
! loop over the number of turbulent variables.
        do l=nt1,nt2
          select case  (bcfaceid(nn)) 
          case (imin) 
            bvti1(i, j, l) = two*bcdata(nn)%turbinlet(i, j, l)
            bmti1(i, j, l, l) = one
          case (imax) 
            bvti2(i, j, l) = two*bcdata(nn)%turbinlet(i, j, l)
            bmti2(i, j, l, l) = one
          case (jmin) 
            bvtj1(i, j, l) = two*bcdata(nn)%turbinlet(i, j, l)
            bmtj1(i, j, l, l) = one
          case (jmax) 
            bvtj2(i, j, l) = two*bcdata(nn)%turbinlet(i, j, l)
            bmtj2(i, j, l, l) = one
          case (kmin) 
            bvtk1(i, j, l) = two*bcdata(nn)%turbinlet(i, j, l)
            bmtk1(i, j, l, l) = one
          case (kmax) 
            bvtk2(i, j, l) = two*bcdata(nn)%turbinlet(i, j, l)
            bmtk2(i, j, l, l) = one
          end select
        end do
      end do
    end do
  end subroutine bcturbinflow
  subroutine bcturboutflow(nn)
!
!       bcturboutflow applies the implicit treatment of the outflow
!       boundary conditions to subface nn. as the outflow boundary
!       condition is independent of the turbulence model, either
!       extrapolation or zero neumann, this routine is valid for all
!       models. it is assumed that the pointers in blockpointers are
!       already set to the correct block on the correct grid level.
!
    use constants
    use blockpointers
    use flowvarrefstate
    implicit none
!
!      subroutine arguments.
!
    integer(kind=inttype), intent(in) :: nn
!
!      local variables.
!
    integer(kind=inttype) :: i, j, l
! loop over the faces of the subfaces and set the values of bmt
! for an implicit treatment. for an outflow the turbulent variable
! variable is either extrapolated or zero neumann. as constant
! extrapolation is used this leads to an identical treatment, i.e.
! the halo value is identical to the value of the internal cell.
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        do l=nt1,nt2
          select case  (bcfaceid(nn)) 
          case (imin) 
            bmti1(i, j, l, l) = -one
          case (imax) 
            bmti2(i, j, l, l) = -one
          case (jmin) 
            bmtj1(i, j, l, l) = -one
          case (jmax) 
            bmtj2(i, j, l, l) = -one
          case (kmin) 
            bmtk1(i, j, l, l) = -one
          case (kmax) 
            bmtk2(i, j, l, l) = -one
          end select
        end do
      end do
    end do
  end subroutine bcturboutflow
  subroutine bcturbtreatment()
!
!       bcturbtreatment sets the arrays bmti1, bvti1, etc, such that
!       the physical boundary conditions are treated correctly.
!       it is assumed that the variables in blockpointers already
!       point to the correct block.
!       the turbulent variable in the halo is computed as follows:
!       whalo = -bmt*winternal + bvt for every block facer. as it is
!       possible to have a coupling in the boundary conditions bmt
!       actually are matrices. if there is no coupling between the
!       boundary conditions of the turbulence equations bmt is a
!       diagonal matrix.
!
    use constants
    use blockpointers
    use flowvarrefstate
    implicit none
!
!      local variable.
!
    integer(kind=inttype) :: nn, i, j, k, l, m
! initialize the arrays for the boundary condition treatment
! to zero, such that internal block boundaries are solved
! correctly (i.e. explicitly).
    do k=1,ke
      do j=1,je
        do l=nt1,nt2
          do m=nt1,nt2
            bmti1(j, k, l, m) = zero
            bmti2(j, k, l, m) = zero
          end do
          bvti1(j, k, l) = zero
          bvti2(j, k, l) = zero
        end do
      end do
    end do
    do k=1,ke
      do i=1,ie
        do l=nt1,nt2
          do m=nt1,nt2
            bmtj1(i, k, l, m) = zero
            bmtj2(i, k, l, m) = zero
          end do
          bvtj1(i, k, l) = zero
          bvtj2(i, k, l) = zero
        end do
      end do
    end do
    do j=1,je
      do i=1,ie
        do l=nt1,nt2
          do m=nt1,nt2
            bmtk1(i, j, l, m) = zero
            bmtk2(i, j, l, m) = zero
          end do
          bvtk1(i, j, l) = zero
          bvtk2(i, j, l) = zero
        end do
      end do
    end do
! loop over the boundary condition subfaces of this block.
bocos:do nn=1,nbocos
! determine the kind of boundary condition for this subface.
      select case  (bctype(nn)) 
      case (nswalladiabatic, nswallisothermal) 
! viscous wall. there is no difference between an adiabatic
! and an isothermal wall for the turbulent equations.
! set the implicit treatment of the wall boundary conditions.
        call bcturbwall(nn)
      case (symm, symmpolar, eulerwall) 
!=============================================================
!=============================================================
! symmetry, polar symmetry or inviscid wall. treatment of
! the turbulent equations is identical.
        call bcturbsymm(nn)
      case (farfield) 
!=============================================================
! farfield. the kind of boundary condition to be applied,
! inflow or outflow, depends on the local conditions.
        call bcturbfarfield(nn)
      case (slidinginterface, oversetouterbound, domaininterfaceall, &
&     domaininterfacerhouvw, domaininterfacep, domaininterfacerho, &
&     domaininterfacetotal) 
!=============================================================
! sliding mesh interface, overset outer boudaries, and
! domain interface with another code are not really boundary
! condition and therefore the values are kept.
        call bcturbinterface(nn)
      end select
    end do bocos
  end subroutine bcturbtreatment
  subroutine bcturbfarfield(nn)
!
!       bcturbfarfield applies the implicit treatment of the
!       farfield boundary condition to subface nn. as the farfield
!       boundary condition is independent of the turbulence model,
!       this routine is valid for all models. it is assumed that the
!       pointers in blockpointers are already set to the correct
!       block on the correct grid level.
!
    use constants
    use blockpointers
    use flowvarrefstate
    implicit none
!
!      subroutine arguments.
!
    integer(kind=inttype), intent(in) :: nn
!
!      local variables.
!
    integer(kind=inttype) :: i, j, l
    real(kind=realtype) :: nnx, nny, nnz, dot
! loop over the faces of the subfaces and set the values of
! bmt and bvt for an implicit treatment.
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
! determine the dot product between the outward pointing
! normal and the free stream velocity direction and add the
! possible grid velocity.
        dot = bcdata(nn)%norm(i, j, 1)*winf(ivx) + bcdata(nn)%norm(i, j&
&         , 2)*winf(ivy) + bcdata(nn)%norm(i, j, 3)*winf(ivz) - bcdata(&
&         nn)%rface(i, j)
! determine whether we are dealing with an inflow or
! outflow boundary here.
        if (dot .gt. zero) then
! outflow. simply extrapolation or zero neumann bc
! of the turbulent variables.
          do l=nt1,nt2
            select case  (bcfaceid(nn)) 
            case (imin) 
              bmti1(i, j, l, l) = -one
            case (imax) 
              bmti2(i, j, l, l) = -one
            case (jmin) 
              bmtj1(i, j, l, l) = -one
            case (jmax) 
              bmtj2(i, j, l, l) = -one
            case (kmin) 
              bmtk1(i, j, l, l) = -one
            case (kmax) 
              bmtk2(i, j, l, l) = -one
            end select
          end do
        else
! inflow. turbulent variables are prescribed.
          do l=nt1,nt2
            select case  (bcfaceid(nn)) 
            case (imin) 
              bvti1(i, j, l) = winf(l)
            case (imax) 
              bvti2(i, j, l) = winf(l)
            case (jmin) 
              bvtj1(i, j, l) = winf(l)
            case (jmax) 
              bvtj2(i, j, l) = winf(l)
            case (kmin) 
              bvtk1(i, j, l) = winf(l)
            case (kmax) 
              bvtk2(i, j, l) = winf(l)
            end select
          end do
        end if
      end do
    end do
  end subroutine bcturbfarfield
  subroutine bcturbinterface(nn)
!
!       bcturbinterface applies the halo treatment for interface halo
!       cells, sliding mesh interface and domain interface. as these
!       are not really boundary conditions, the variable bvt is simply
!       set to keep the current value.
!
    use constants
    use blockpointers
    use flowvarrefstate
    implicit none
! note that the original code had an error in the pointers...they
! were pointing to {il,jl,kl} and not {ie, je, ke}.
!
!      subroutine arguments.
!
    integer(kind=inttype), intent(in) :: nn
!
!      local variables.
!
    integer(kind=inttype) :: i, j, l
! loop over the faces of the subfaces and set the values of
! bvt to keep the current value.
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        do l=nt1,nt2
          select case  (bcfaceid(nn)) 
          case (imin) 
            bvti1(i, j, l) = w(1, i, j, l)
          case (imax) 
            bvti2(i, j, l) = w(ie, i, j, l)
          case (jmin) 
            bvtj1(i, j, l) = w(i, 1, j, l)
          case (jmax) 
            bvtj2(i, j, l) = w(i, je, j, l)
          case (kmin) 
            bvtk1(i, j, l) = w(i, j, 1, l)
          case (kmax) 
            bvtk2(i, j, l) = w(i, j, ke, l)
          end select
        end do
      end do
    end do
  end subroutine bcturbinterface
  subroutine bcturbsymm(nn)
!
!       bcturbsymm applies the implicit treatment of the symmetry
!       boundary condition (or inviscid wall) to subface nn. as the
!       symmetry boundary condition is independent of the turbulence
!       model, this routine is valid for all models. it is assumed
!       that the pointers in blockpointers are already set to the
!       correct block on the correct grid level.
!
    use constants
    use blockpointers
    use flowvarrefstate
    implicit none
!
!      subroutine arguments.
!
    integer(kind=inttype), intent(in) :: nn
!
!      local variables.
!
    integer(kind=inttype) :: i, j, l
! loop over the faces of the subfaces and set the values of bmt
! for an implicit treatment. for a symmetry face this means
! that the halo value is set to the internal value.
    do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
      do i=bcdata(nn)%icbeg,bcdata(nn)%icend
        do l=nt1,nt2
          select case  (bcfaceid(nn)) 
          case (imin) 
            bmti1(i, j, l, l) = -one
          case (imax) 
            bmti2(i, j, l, l) = -one
          case (jmin) 
            bmtj1(i, j, l, l) = -one
          case (jmax) 
            bmtj2(i, j, l, l) = -one
          case (kmin) 
            bmtk1(i, j, l, l) = -one
          case (kmax) 
            bmtk2(i, j, l, l) = -one
          end select
        end do
      end do
    end do
  end subroutine bcturbsymm
  subroutine turb2ndhalo(nn)
!
!       turb2ndhalo sets the turbulent variables in the second halo
!       cell for the given subface. simple constant extrapolation is
!       used to avoid problems.
!
    use constants
    use blockpointers
    use flowvarrefstate
    implicit none
!
!      subroutine arguments.
!
    integer(kind=inttype), intent(in) :: nn
!
!      local variables.
!
    integer(kind=inttype) :: i, j, l
! determine the face on which this subface is located and set
! some pointers accordingly.
! loop over the turbulent variables and set the second halo
! value. if this is an eddy model, also set the eddy viscosity.
    select case  (bcfaceid(nn)) 
    case (imin) 
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          do l=nt1,nt2
            w(0, i, j, l) = w(1, i, j, l)
          end do
          if (eddymodel) rev(0, i, j) = rev(1, i, j)
        end do
      end do
    case (imax) 
!===============================================================
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          do l=nt1,nt2
            w(ib, i, j, l) = w(ie, i, j, l)
          end do
          if (eddymodel) rev(ib, i, j) = rev(ie, i, j)
        end do
      end do
    case (jmin) 
!===============================================================
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          do l=nt1,nt2
            w(i, 0, j, l) = w(i, 1, j, l)
          end do
          if (eddymodel) rev(i, 0, j) = rev(i, 1, j)
        end do
      end do
    case (jmax) 
!===============================================================
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          do l=nt1,nt2
            w(i, jb, j, l) = w(i, je, j, l)
          end do
          if (eddymodel) rev(i, jb, j) = rev(i, je, j)
        end do
      end do
    case (kmin) 
!===============================================================
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          do l=nt1,nt2
            w(i, j, 0, l) = w(i, j, 1, l)
          end do
          if (eddymodel) rev(i, j, 0) = rev(i, j, 1)
        end do
      end do
    case (kmax) 
!===============================================================
      do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
        do i=bcdata(nn)%icbeg,bcdata(nn)%icend
          do l=nt1,nt2
            w(i, j, kb, l) = w(i, j, ke, l)
          end do
          if (eddymodel) rev(i, j, kb) = rev(i, j, ke)
        end do
      end do
    end select
  end subroutine turb2ndhalo
  subroutine turbbcnswall(secondhalo)
!
!       turbbcnswall applies the viscous wall boundary conditions
!       of the turbulent transport equations to a block. it is assumed
!       that the pointers in blockpointers are already set to the
!       correct block on the correct grid level.
!
    use constants
    use blockpointers
    use flowvarrefstate
    implicit none
!
!      subroutine argument.
!
    logical, intent(in) :: secondhalo
!
!      local variables.
!
    integer(kind=inttype) :: nn, i, j, l, m
! loop over the viscous subfaces of this block.
bocos:do nn=1,nviscbocos
! set the corresponding arrays.
      call bcturbwall(nn)
! loop over the faces and set the state in
! the turbulent halo cells.
      select case  (bcfaceid(nn)) 
      case (imin) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            do l=nt1,nt2
              w(1, i, j, l) = bvti1(i, j, l)
              do m=nt1,nt2
                w(1, i, j, l) = w(1, i, j, l) - bmti1(i, j, l, m)*w(2, i&
&                 , j, m)
              end do
              if (secondhalo) w(0, i, j, l) = w(1, i, j, l)
            end do
            if (eddymodel) then
              rev(1, i, j) = -rev(2, i, j)
              if (secondhalo) rev(0, i, j) = rev(1, i, j)
            end if
          end do
        end do
      case (imax) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            do l=nt1,nt2
              w(ie, i, j, l) = bvti2(i, j, l)
              do m=nt1,nt2
                w(ie, i, j, l) = w(ie, i, j, l) - bmti2(i, j, l, m)*w(il&
&                 , i, j, m)
              end do
              if (secondhalo) w(ib, i, j, l) = w(ie, i, j, l)
            end do
            if (eddymodel) then
              rev(ie, i, j) = -rev(il, i, j)
              if (secondhalo) rev(ib, i, j) = rev(ie, i, j)
            end if
          end do
        end do
      case (jmin) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            do l=nt1,nt2
              w(i, 1, j, l) = bvtj1(i, j, l)
              do m=nt1,nt2
                w(i, 1, j, l) = w(i, 1, j, l) - bmtj1(i, j, l, m)*w(i, 2&
&                 , j, m)
              end do
              if (secondhalo) w(i, 0, j, l) = w(i, 1, j, l)
            end do
            if (eddymodel) then
              rev(i, 1, j) = -rev(i, 2, j)
              if (secondhalo) rev(i, 0, j) = rev(i, 1, j)
            end if
          end do
        end do
      case (jmax) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            do l=nt1,nt2
              w(i, je, j, l) = bvtj2(i, j, l)
              do m=nt1,nt2
                w(i, je, j, l) = w(i, je, j, l) - bmtj2(i, j, l, m)*w(i&
&                 , jl, j, m)
              end do
              if (secondhalo) w(i, jb, j, l) = w(i, je, j, l)
            end do
            if (eddymodel) then
              rev(i, je, j) = -rev(i, jl, j)
              if (secondhalo) rev(i, jb, j) = rev(i, je, j)
            end if
          end do
        end do
      case (kmin) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            do l=nt1,nt2
              w(i, j, 1, l) = bvtk1(i, j, l)
              do m=nt1,nt2
                w(i, j, 1, l) = w(i, j, 1, l) - bmtk1(i, j, l, m)*w(i, j&
&                 , 2, m)
              end do
              if (secondhalo) w(i, j, 0, l) = w(i, j, 1, l)
            end do
            if (eddymodel) then
              rev(i, j, 1) = -rev(i, j, 2)
              if (secondhalo) rev(i, j, 0) = rev(i, j, 1)
            end if
          end do
        end do
      case (kmax) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            do l=nt1,nt2
              w(i, j, ke, l) = bvtk2(i, j, l)
              do m=nt1,nt2
                w(i, j, ke, l) = w(i, j, ke, l) - bmtk2(i, j, l, m)*w(i&
&                 , j, kl, m)
              end do
              if (secondhalo) w(i, j, kb, l) = w(i, j, ke, l)
            end do
            if (eddymodel) then
              rev(i, j, ke) = -rev(i, j, kl)
              if (secondhalo) rev(i, j, kb) = rev(i, j, ke)
            end if
          end do
        end do
      end select
    end do bocos
  end subroutine turbbcnswall
  subroutine bcturbwall(nn)
!
!       bcturbwall applies the implicit treatment of the viscous
!       wall boundary condition for the turbulence model used to the
!       given subface nn.
!       it is assumed that the pointers in blockpointers are
!       already set to the correct block.
!
    use blockpointers
    use flowvarrefstate
    use inputphysics
    use constants
    use paramturb
    implicit none
!
!      subroutine arguments.
!
    integer(kind=inttype), intent(in) :: nn
!
!      local variables.
!
    integer(kind=inttype) :: i, j, ii, jj, iimax, jjmax
    real(kind=realtype) :: tmpd, tmpe, tmpf, nu
    real(kind=realtype), dimension(:, :, :, :), pointer :: bmt
    real(kind=realtype), dimension(:, :, :), pointer :: bvt, ww2
    real(kind=realtype), dimension(:, :), pointer :: rlv2, dd2wall
    intrinsic min
    intrinsic max
    real(kind=realtype) :: result1
    integer(kind=inttype) :: y12
    integer(kind=inttype) :: y11
    integer(kind=inttype) :: y10
    integer(kind=inttype) :: y9
    integer(kind=inttype) :: y8
    integer(kind=inttype) :: y7
    integer(kind=inttype) :: y6
    integer(kind=inttype) :: y5
    integer(kind=inttype) :: y4
    integer(kind=inttype) :: y3
    integer(kind=inttype) :: y2
    integer(kind=inttype) :: y1
!        ================================================================
! determine the turbulence model used and loop over the faces
! of the subface and set the values of bmt and bvt for an
! implicit treatment.
    select case  (turbmodel) 
    case (spalartallmaras, spalartallmarasedwards) 
! spalart-allmaras type of model. value at the wall is zero,
! so simply negate the internal value.
      select case  (bcfaceid(nn)) 
      case (imin) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            result1 = saroughfact(2, i, j)
            bmti1(i, j, itu1, itu1) = -result1
          end do
        end do
      case (imax) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            result1 = saroughfact(il, i, j)
            bmti2(i, j, itu1, itu1) = -result1
          end do
        end do
      case (jmin) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            result1 = saroughfact(i, 2, j)
            bmtj1(i, j, itu1, itu1) = -result1
          end do
        end do
      case (jmax) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            result1 = saroughfact(i, jl, j)
            bmtj2(i, j, itu1, itu1) = -result1
          end do
        end do
      case (kmin) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            result1 = saroughfact(i, j, 2)
            bmtk1(i, j, itu1, itu1) = -result1
          end do
        end do
      case (kmax) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            result1 = saroughfact(i, j, kl)
            bmtk2(i, j, itu1, itu1) = -result1
          end do
        end do
      end select
    case (komegawilcox, komegamodified, mentersst) 
!        ================================================================
! k-omega type of models. k is zero on the wall and thus the
! halo value is the negative of the first internal cell.
! for omega the situation is a bit more complicated.
! theoretically omega is infinity, but it is set to a large
! value, see menter's paper. the halo value is constructed
! such that the wall value is correct. make sure that i and j
! are limited to physical dimensions of the face for the wall
! distance. due to the usage of the dd2wall pointer and the
! fact that the original d2wall array starts at 2, there is
! an offset of -1 present in dd2wall.
      select case  (bcfaceid(nn)) 
      case (imin) 
        iimax = jl
        jjmax = kl
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          if (j .gt. jjmax) then
            y1 = jjmax
          else
            y1 = j
          end if
          if (2 .lt. y1) then
            jj = y1
          else
            jj = 2
          end if
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            if (i .gt. iimax) then
              y2 = iimax
            else
              y2 = i
            end if
            if (2 .lt. y2) then
              ii = y2
            else
              ii = 2
            end if
            nu = rlv(2, i, j)/w(2, i, j, irho)
            tmpd = one/(rkwbeta1*d2wall(2, ii, jj)**2)
            bmti1(i, j, itu1, itu1) = one
            bmti1(i, j, itu2, itu2) = one
            bvti1(i, j, itu2) = two*60.0_realtype*nu*tmpd
          end do
        end do
      case (imax) 
        iimax = jl
        jjmax = kl
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          if (j .gt. jjmax) then
            y3 = jjmax
          else
            y3 = j
          end if
          if (2 .lt. y3) then
            jj = y3
          else
            jj = 2
          end if
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            if (i .gt. iimax) then
              y4 = iimax
            else
              y4 = i
            end if
            if (2 .lt. y4) then
              ii = y4
            else
              ii = 2
            end if
            nu = rlv(jl, i, j)/w(il, i, j, irho)
            tmpd = one/(rkwbeta1*d2wall(il, ii, jj)**2)
            bmti2(i, j, itu1, itu1) = one
            bmti2(i, j, itu2, itu2) = one
            bvti2(i, j, itu2) = two*60.0_realtype*nu*tmpd
          end do
        end do
      case (jmin) 
        iimax = il
        jjmax = kl
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          if (j .gt. jjmax) then
            y5 = jjmax
          else
            y5 = j
          end if
          if (2 .lt. y5) then
            jj = y5
          else
            jj = 2
          end if
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            if (i .gt. iimax) then
              y6 = iimax
            else
              y6 = i
            end if
            if (2 .lt. y6) then
              ii = y6
            else
              ii = 2
            end if
            nu = rlv(i, 2, j)/w(i, 2, j, irho)
            tmpd = one/(rkwbeta1*d2wall(ii, 2, jj)**2)
            bmtj1(i, j, itu1, itu1) = one
            bmtj1(i, j, itu2, itu2) = one
            bvtj1(i, j, itu2) = two*60.0_realtype*nu*tmpd
          end do
        end do
      case (jmax) 
        iimax = il
        jjmax = kl
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          if (j .gt. jjmax) then
            y7 = jjmax
          else
            y7 = j
          end if
          if (2 .lt. y7) then
            jj = y7
          else
            jj = 2
          end if
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            if (i .gt. iimax) then
              y8 = iimax
            else
              y8 = i
            end if
            if (2 .lt. y8) then
              ii = y8
            else
              ii = 2
            end if
            nu = rlv(i, jl, j)/w(i, jl, j, irho)
            tmpd = one/(rkwbeta1*d2wall(ii, jl, jj)**2)
            bmtj2(i, j, itu1, itu1) = one
            bmtj2(i, j, itu2, itu2) = one
            bvtj2(i, j, itu2) = two*60.0_realtype*nu*tmpd
          end do
        end do
      case (kmin) 
        iimax = il
        jjmax = jl
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          if (j .gt. jjmax) then
            y9 = jjmax
          else
            y9 = j
          end if
          if (2 .lt. y9) then
            jj = y9
          else
            jj = 2
          end if
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            if (i .gt. iimax) then
              y10 = iimax
            else
              y10 = i
            end if
            if (2 .lt. y10) then
              ii = y10
            else
              ii = 2
            end if
            nu = rlv(i, j, 2)/w(i, j, 2, irho)
            tmpd = one/(rkwbeta1*d2wall(ii, jj, 2)**2)
            bmtk1(i, j, itu1, itu1) = one
            bmtk1(i, j, itu2, itu2) = one
            bvtk1(i, j, itu2) = two*60.0_realtype*nu*tmpd
          end do
        end do
      case (kmax) 
        iimax = il
        jjmax = jl
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          if (j .gt. jjmax) then
            y11 = jjmax
          else
            y11 = j
          end if
          if (2 .lt. y11) then
            jj = y11
          else
            jj = 2
          end if
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            if (i .gt. iimax) then
              y12 = iimax
            else
              y12 = i
            end if
            if (2 .lt. y12) then
              ii = y12
            else
              ii = 2
            end if
            nu = rlv(i, j, kl)/w(i, j, kl, irho)
            tmpd = one/(rkwbeta1*d2wall(ii, jj, kl)**2)
            bmtk2(i, j, itu1, itu1) = one
            bmtk2(i, j, itu2, itu2) = one
            bvtk2(i, j, itu2) = two*60.0_realtype*nu*tmpd
          end do
        end do
      end select
    case (ktau) 
!        ================================================================
! k-tau model. both k and tau are zero at the wall, so the
! negative value of the internal cell is taken for the halo.
      select case  (bcfaceid(nn)) 
      case (imin) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            bmti1(i, j, itu1, itu1) = one
            bmti1(i, j, itu2, itu2) = one
          end do
        end do
      case (imax) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            bmti2(i, j, itu1, itu1) = one
            bmti2(i, j, itu2, itu2) = one
          end do
        end do
      case (jmin) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            bmtj1(i, j, itu1, itu1) = one
            bmtj1(i, j, itu2, itu2) = one
          end do
        end do
      case (jmax) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            bmtj2(i, j, itu1, itu1) = one
            bmtj2(i, j, itu2, itu2) = one
          end do
        end do
      case (kmin) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            bmtk1(i, j, itu1, itu1) = one
            bmtk1(i, j, itu2, itu2) = one
          end do
        end do
      case (kmax) 
        do j=bcdata(nn)%jcbeg,bcdata(nn)%jcend
          do i=bcdata(nn)%icbeg,bcdata(nn)%icend
            bmtk2(i, j, itu1, itu1) = one
            bmtk2(i, j, itu2, itu2) = one
          end do
        end do
      end select
    end select
  end subroutine bcturbwall
  function saroughfact(i, j, k)
! returns either the regular sa-boundary condition
! or the modified roughness-boundary condition
    use constants
    use inputphysics, only : useroughsa
    use blockpointers, only : ks, d2wall
    implicit none
! dummy arguments
    real(kind=realtype) :: saroughfact
! local variablse
    integer(kind=inttype) :: i, j, k
    if (.not.useroughsa) then
      saroughfact = -one
      return
    else
      saroughfact = (ks(i, j, k)-d2wall(i, j, k)/0.03_realtype)/(ks(i, j&
&       , k)+d2wall(i, j, k)/0.03_realtype)
    end if
  end function saroughfact
end module turbbcroutines_fast_b
