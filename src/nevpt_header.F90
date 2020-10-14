module nevpt_header

implicit none

#include "version.h" 

public print_nevpt_header

contains

subroutine print_nevpt_header(rel_ham)

  logical, intent(in) :: rel_ham
!
!     This program calculates the second order correction to the energy
!     using the n-electron valence state perturbation theory (NEVPT)
!     in the partially contracted (PC) ad strongly contracted (SC)
!     variants using the Dyall's Hamiltonian for the definition of the
!     zero order energies. All equations are based on the spin free
!     formulation of NEVPT.
!
!     Relevant papers:
!     1) C. Angeli, R. Cimiraglia, J-P Malrieu, 
!        Chem. Phys. Lett.,  317, 472-480, (2000)
!        Definition of the inactive (core+virtual) orbital energies.
!     2) C. Angeli, R. Cimiraglia, S. Evangelisti, T. Leininger, J-P Malrieu,
!        J. Chem. Phys., 114, 10252, (2001)
!        Introduction of NEVPT
!     3) C. Angeli, R. Cimiraglia, J-P Malrieu,
!        Chem. Phys. Lett., 350, 297-305, (2001)
!        Introduction of the formalism based on the "Koopmans' operators"
!     4) C. Angeli, R. Cimiraglia, J-P Malrieu,
!        J. Chem. Phys., 117, 10525-10533, (2002)
!        Spin free version of NEVPT

print '(/a)', ' ---------------------------------------------------------------------- '
print '(a )',' This program calculates the second order correction to the energy       '
print '(a )',' using the n-electron valence state perturbation theory (NEVPT)          '
print '(a )',' in the partially contracted (PC) ad strongly contracted (SC)            '
print '(a )',' variants using the Dyalls Hamiltonian for the definition of the        '
print '(a )',' zero order energies. All equations are based on the spin free           '
print '(a/)',' formulation of NEVPT.                                                   '
print '(a )',' Relevant papers:                                                        '
print '(a )',' 1) C. Angeli, R. Cimiraglia, J-P Malrieu,                               '
print '(a )','    Chem. Phys. Lett.,  317, 472-480, (2000)                             '
print '(a )','    Definition of the inactive (core+virtual) orbital energies.          '
print '(a )',' 2) C. Angeli, R. Cimiraglia, S. Evangelisti, T. Leininger, J-P Malrieu, '
print '(a )','    J. Chem. Phys., 114, 10252, (2001)                                   '
print '(a )','    Introduction of NEVPT                                                '
print '(a )',' 3) C. Angeli, R. Cimiraglia, J-P Malrieu,                               '
print '(a )','    Chem. Phys. Lett., 350, 297-305, (2001)                              '
print '(a )','    Introduction of the formalism based on the Koopmans operators        '
print '(a )',' 4) C. Angeli, R. Cimiraglia, J-P Malrieu,                               '
print '(a )','    J. Chem. Phys., 117, 10525-10533, (2002)                             '
print '(a/)','    Spin free version of NEVPT                                           '
print '(a )',' Main authors of the original NEVPT code:                                '
print '(a )','                  C. Angeli     (U Ferrara, Italy)                       '
print '(a )','                  R. Cimiraglia (U Ferrara, Italy)                       '
print '(a )','                  J-P Malrieu   (U Toulouse, France)                     '
!print '(a )',' Contributing author:                                                   '
!print '(a )','                  L. Tenti      (U Ferrara, Italy)                       '
print '(a )',' DMRG-driven NEVPT code:                                                 '
print '(a )','                  S. Knecht     (ETH Zurich)                             '
print '(a )','                  L. Freitag    (ETH Zurich)                             '
print '(a )','                  M. Reiher     (ETH Zurich)                             '
print '(a )',' DMRG-driven/CASSCF-driven relativistic NEVPT code:                      '
print '(a )','                  S. Knecht     (ETH Zurich)                             '
print '(/a,a )','                NEVPT version: ',NEVPT_VERSION
print '( a,a/)','                 Git revision: ',NEVPT_GIT_VERSION
if(rel_ham)then
print '(a/)','  Running the relativistic version                                       '
else
print '(a/)','  Running the non- and scalar-relativistic version                       '
print '( a)','  Please cite the following work in a publication that reports data      '
print '(a/)','  obtained with this code:                                               '
print '(a )','  L. Freitag, S. Knecht, C. Angeli, and M. Reiher,                       '
print '(a )','  J. Chem. Theory Comput., 13, 451-459 (2017).                           '
end if
print '(a/)', ' -----------------------------------------------------------------------'

end subroutine print_nevpt_header

end module nevpt_header
