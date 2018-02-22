module interface_cosmomc
  !use, intrinsic :: iso_c_binding
  use ISO_C_BINDING
  implicit none

  interface
    function theta(omegam, orad, ombh2, h0, zstar, c2, c3, c4, cG, grhormass, nu_masses, nu_mass_eigenstates) bind(C, name='theta_')
      import C_DOUBLE, C_INT, C_PTR ! Make iso c binding visible here
      real(kind=C_DOUBLE) :: omegam, orad, ombh2, h0, zstar, c2, c3, c4, cG
      integer(C_INT) :: nu_mass_eigenstates
      type(C_PTR), value :: grhormass, nu_masses
      real(C_DOUBLE) :: theta
    end function theta
    
 end interface
end module interface_cosmomc
