      module neklab
      ! --> Abstract Krylov methods.
         use LightKrylov
      ! --> Definition of the abstract vectors in the Nek framework.
         use neklab_vectors
      ! --> Utility functions for Nek vectors
         use neklab_utils
      ! --> Definitions of the abstract linops in the Nek framework.
         use neklab_linops
      ! --> Definitions of the abstract systems in the Nek framework.
         use neklab_systems
      ! --> Stability analysis routines
         use neklab_analysis

         ! skern
         use neklab_pvectors
         use neklab_plinops
      
         private
      
      !------------------------------------------
      !-----     LIGHTKRYLOV RE-EXPORTS     -----
      !------------------------------------------
      
      ! Global variables.
         public :: dp, rtol_dp, atol_dp
      
      ! Krylov factorizations.
         public :: initialize_krylov_subspace
         public :: arnoldi_factorization, krylov_schur_restart
         public :: lanczos_bidiagonalization, lanczos_tridiagonalization
      
      ! Matrix factorizations.
         public :: eigs, eighs, svds
         public :: save_eigenspectrum

      ! Linear solvers.
         public :: cg
      ! Auxiliary exports
         public :: cg_dp_opts, cg_dp_metadata
      
      !----------------------------------
      !-----     NEKLAB EXPORTS     -----
      !----------------------------------
      
      ! Definition of the abstract vectors in the Nek framework.
         public :: nek_dvector
         public :: nek_pdvector
         public :: nek_ext_dvector
      
      ! Implementation of the standard linear operators.
         public :: exptA_linop, LNS_linop, DTD_linop
      
      ! Implementation of the abstract systems and Jacobians
         public :: nek_system, nek_system_upo
         public :: nek_jacobian, nek_jacobian_upo
      
      ! Baseflow computation
         public :: newton_fixed_point_iteration
         public :: newton_periodic_orbit
      
      ! Stability analysis exports.
         public :: linear_stability_analysis_fixed_point
         public :: transient_growth_analysis_fixed_point

      ! Miscellaneous
         public :: compare_nek_arnoldi
      
      ! Various utilities.
         public :: nek2vec, vec2nek
         public :: nek2pvec, pvec2nek
         public :: nek2ext_vec, ext_vec2nek
         public :: setup_nonlinear_solver, setup_linear_solver
         public :: apply_L

      end module neklab
