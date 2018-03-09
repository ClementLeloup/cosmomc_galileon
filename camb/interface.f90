module interface
  !use, intrinsic :: iso_c_binding
  use ISO_C_BINDING
  implicit none

  interface
     !from galileon.cc file
     function arrays(omegar, omegam, H0in, c2in, c3in, c4in, cGin, grhormass, nu_masses, nu_mass_eigenstates) bind(C, name='arrays_')
       import C_DOUBLE, C_INT, C_PTR ! Make iso c binding visible here
       real(kind=C_DOUBLE) :: omegar, omegam, H0in, c2in, c3in, c4in, cGin
       integer(C_INT) :: nu_mass_eigenstates
       type(C_PTR), value :: grhormass, nu_masses
       integer :: arrays
    end function arrays

    function GetX(point) bind(C, name='GetX_')
      import C_DOUBLE, C_PTR
      real(kind=C_DOUBLE) :: point, GetX
    end function GetX

    function GetH(point) bind(C, name='GetH_')
      import C_DOUBLE, C_PTR
      real(kind=C_DOUBLE) :: point, GetH
    end function GetH

    subroutine GetdHdX(point, hcamb, xcamb, dh, dx, grhormass, nu_masses, nu_mass_eigenstates) bind(C, name='GetdHdX_')
      import C_PTR, C_DOUBLE, C_INT
      real(kind=C_DOUBLE) :: point, hcamb, xcamb, dh, dx
      integer(C_INT) :: nu_mass_eigenstates
      type(C_PTR), value :: grhormass, nu_masses
    end subroutine GetdHdX

    function ct(point, grhormass, nu_masses, nu_mass_eigenstates) bind(C, name='ct_')
      import C_DOUBLE, C_PTR, C_INT
      real(kind=C_DOUBLE) :: point, ct
      integer(C_INT) :: nu_mass_eigenstates
      type(C_PTR), value :: grhormass, nu_masses
    end function ct

    function grhogal(point, hcamb, xcamb) bind(C, name='grhogal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: point, grhogal, hcamb, xcamb
    end function grhogal

    function gpresgal(point, hcamb, xcamb, dhcamb, dxcamb) bind(C, name='gpresgal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: gpresgal, point, hcamb, xcamb, dhcamb, dxcamb
      end function gpresgal

    function Chigal(point, hcamb, xcamb, dgrho, eta, dphi, dphiprime, k) bind(C, name='Chigal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: point, hcamb, xcamb
      real(kind=C_DOUBLE) :: dgrho, eta, dphi, dphiprime, k, Chigal
    end function Chigal

    function qgal(point, hcamb, xcamb, dgq, eta, dphi, dphiprime, k) bind(C, name='qgal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: point, hcamb, xcamb
      real(kind=C_DOUBLE) :: dgq, eta, dphi, dphiprime, k, qgal
    end function qgal

    function Pigal(point, hcamb, xcamb, dhcamb, dxcamb, dgrho, dgq, dgpi, eta, dphi, k) bind(C, name='Pigal_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: point, hcamb, xcamb, dhcamb, dxcamb
      real(kind=C_DOUBLE) :: dgrho, dgq, dgpi, eta, dphi, k, Pigal
    end function Pigal

    function dphisecond(point, hcamb, xcamb, dhcamb, dxcamb, dgrho, dgq, eta, dphi, dphiprime, k, deltafprime) bind(C, name='dphisecond_')
      import C_DOUBLE
      real(kind=C_DOUBLE) :: point, hcamb, xcamb, dhcamb, dxcamb
      real(kind=C_DOUBLE) :: dgrho, dgq, eta, dphi, dphiprime, k, dphisecond, deltafprime
    end function dphisecond

    function pigalprime(point, hcamb, xcamb, dhcamb, dxcamb, dgrho, dgq, dgpi, pidot, eta, dphi, dphiprime, k, grho, gpres, grhormass, nu_masses, nu_mass_eigenstates) bind(C, name='pigalprime_')
      import C_DOUBLE, C_INT, C_PTR
      real(kind=C_DOUBLE) :: point, hcamb, xcamb, dhcamb, dxcamb
      real(kind=C_DOUBLE) :: dgrho, dgq, dgpi, pidot, eta, dphi, dphiprime, k, grho, gpres, pigalprime
      integer(C_INT) :: nu_mass_eigenstates
      type(C_PTR), value :: grhormass, nu_masses
    end function pigalprime

    function crosschecks(point, hcamb, xcamb, dhcamb, dxcamb, dgrho, dgq, dgpi, eta, dphi, dphiprime, dphiprimeprime, k, grho, gpres, deltafprime) bind(C, name='crosschecks_')
      import C_PTR, C_DOUBLE
      real(kind=C_DOUBLE) :: point, hcamb, xcamb, dhcamb, dxcamb
      real(kind=C_DOUBLE) :: dgrho, dgq, dgpi, eta, dphi, dphiprime, dphiprimeprime, k, grho, gpres, deltafprime
      type(C_PTR) :: crosschecks
    end function crosschecks

    subroutine freegal() bind(C, name='freegal_')
      !subroutine to free non-necessary memory after perturbation calculations
    end subroutine freegal
    
 end interface
end module interface
