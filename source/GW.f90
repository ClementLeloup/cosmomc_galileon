! Module created by Clement Leloup to use Gravitational Waves likelihood

    module GW
    use MatrixUtils
    use settings
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    use IniObjects
    implicit none
    private

    type, extends(TCosmoCalcLikelihood) :: TGWLikelihood
        integer :: num_gw ! total number of points used
        integer, allocatable :: type_bao(:) !one of the constants defined above
        real(mcp) :: rs_rescale = 1._mcp !if not generated with numerical CAMB rs
        real(mcp), allocatable, dimension(:) :: gw_z, gw_obs, gw_err
        real(mcp), allocatable, dimension(:,:) :: gw_invcov
        real(mcp) :: Hrd_fid,DA_rd_fid
    contains
    procedure :: LogLike => GW_LnLike
    procedure :: ReadIni => GW_ReadIni
    procedure :: InitProbDist => GW_InitProbDist
    procedure, private :: deltat
    end type TGWLikelihood

    public GWLikelihood_Add
    contains

    subroutine GWLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    class(TGWLikelihood), pointer :: this
    integer i
    Type(TSettingIni) :: DataSets, OverrideSettings

    if (.not. Ini%Read_Logical('use_GW',.false.)) return

    call Ini%TagValuesForName('gw_dataset', DataSets, filename=.true.)
    if (DataSets%Count==0) call MpiStop('Use_GW but no gw_dataset[NAMETAG] defined')

    do i= 1, DataSets%Count
        call Ini%SettingValuesForTagName('gw_dataset',DataSets%Name(i),OverrideSettings)
        allocate(TGWLikelihood::this)
        call this%ReadDatasetFile(Datasets%Value(i),OverrideSettings)
        this%tag = Datasets%Name(i)
        this%LikelihoodType = 'GW'
        this%needs_background_functions = .true.
        call LikeList%Add(this)
    end do

    call this%loadParamNames(trim(DataDir)//'GW.paramnames')

    if (Feedback>1) write(*,*) 'read GW data sets'

    end subroutine GWLikelihood_Add

    subroutine GW_ReadIni(this, Ini)
    class(TGWLikelihood) this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: gw_measurement, gw_measurements_file
    integer i,iopb
    Type(TTextFile) :: F
    logical :: hastype, haserror
    integer :: status
    character(LEN=:), allocatable :: InLine

    if (Feedback > 0 .and. MpiRank==0) write (*,*) 'reading GW data set: '//trim(this%name)
    this%num_gw = Ini%Read_Int('num_gw',1)

    allocate(this%gw_z(this%num_gw))
    allocate(this%gw_obs(this%num_gw))
    allocate(this%gw_err(this%num_gw))

    gw_measurements_file = Ini%ReadRelativeFileName('gw_measurements_file')
    haserror = Ini%Read_Logical('gw_measurements_file_has_error',.true.)
    call F%Open(gw_measurements_file)
    do i=1,this%num_gw
       if (F%ReadLineSkipEmptyAndComments(InLine)) then
          if (haserror) then
             read (InLine,*, iostat=iopb) this%gw_z(i),this%gw_obs(i),this%gw_err(i)
          else
             read (InLine,*, iostat=iopb) this%gw_z(i),this%gw_obs(i)
          end if
          if (iopb /= 0) call MpiStop('GW_ReadIni: Error reading gw_measurements_file: ' &
               //trim(this%name))
       else
          call MpiStop('GW_ReadIni: Missing line in gw_measurements_file: ' &
               //trim(this%name))
       end if
    end do
    call F%Close()
        
    !Not sure about this    
    !if (any(this%gw_z< 0.0001)) call MpiStop('Error reading GW measurements')

    call this%InitProbDist(Ini)

    end subroutine GW_ReadIni

    subroutine GW_InitProbDist(this, Ini)
    class(TGWLikelihood) this
    class(TSettingIni) :: Ini
    integer i
    character(LEN=:), allocatable ::gw_invcov_file

    allocate(this%gw_invcov(this%num_gw,this%num_gw))
    this%gw_invcov=0

    if (Ini%HasKey('gw_invcov_file')) then
        gw_invcov_file  = Ini%ReadRelativeFileName('gw_invcov_file')
        call File%ReadTextMatrix(gw_invcov_file, this%gw_invcov)
    else if (Ini%HasKey('gw_cov_file')) then
        gw_invcov_file  = Ini%ReadRelativeFileName('gw_cov_file')
        call File%ReadTextMatrix(gw_invcov_file, this%gw_invcov)
        call Matrix_Inverse(this%gw_invcov)
    else
        do i=1,this%num_gw
            !diagonal, or actually just 1..
            this%gw_invcov(i,i) = 1/this%gw_err(i)**2
        end do
    end if

    end subroutine GW_InitProbDist

    function deltat(this, z)
    Class(TGWLikelihood) :: this
    real(mcp) deltat, z

    if(CosmoSettings%use_galileon)then
       deltat = this%Calculator%GW_light_dt(z)
    else
       deltat = 0
    end if

    end function deltat

    function GW_LnLike(this, CMB, Theory, DataParams)
    Class(TGWLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    integer j
    real(mcp) GW_LnLike, z
    real(mcp)  :: GW_theory(this%num_gw), delta_t

    delta_t=DataParams(1)

    do j=1, this%num_gw
        z= this%gw_z(j)
        GW_theory(j) = this%deltat(z)+delta_t
        print *, "dt=", GW_theory(j), "dt_obs=", this%gw_obs
    end do

    GW_theory = GW_theory - this%gw_obs
    GW_LnLike = Matrix_QuadForm(this%gw_invcov,GW_theory) / 2

    if(feedback>1) write(*,*) trim(this%name)//' GW likelihood = ', GW_LnLike

    end function GW_LnLike

    end module GW
