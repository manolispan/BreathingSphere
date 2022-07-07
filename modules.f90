   module params
      INTEGER, PARAMETER :: DBL= selected_real_kind(14,200)
      INTEGER, PARAMETER :: DP =  kind(0.0d0)
      INTEGER, PARAMETER :: QP =  kind(0.0d0)  ! 8 
      INTEGER, PARAMETER :: CQP = kind(0.0d0) ! 8
      integer, parameter :: LMAXD = 20
      complex(kind=dp),allocatable :: han(:,:),bes(:,:),dhan(:,:),dbes(:,:) 
      complex(kind=dp),allocatable :: SCSbrsph(:)  
      complex(kind=dp),allocatable :: dynte(:,:,:),dyntm(:,:,:)  
      complex(kind=dp),allocatable :: dyntmat(:,:)      
    end module params

    module constants
      use params,only: dp
      complex(dp):: cone = (1.d0,0.d0)
      complex(dp):: ctwo = (2.d0,0.d0)
      complex(dp):: ci = (0.d0,1.d0)
      complex(dp):: czero = (0.d0,0.d0)
    end module constants
