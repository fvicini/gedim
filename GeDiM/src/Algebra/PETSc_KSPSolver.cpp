#include "PETSc_KSPSolver.hpp"

#include "PETSc_Array.hpp"
#include "PETSc_SparseArray.hpp"

using namespace Eigen;
using namespace std;

namespace Gedim
{
  // ***************************************************************************
  template class PETSc_KSPSolver<Vec, Mat, PETSc_SolverTypes::PETSc_KSPCG, PETSc_Preconditioners::PETSc_DEFAULT>;
  template class PETSc_KSPSolver<Vec, Mat, PETSc_SolverTypes::PETSc_KSPCG, PETSc_Preconditioners::PETSc_PCJACOBI>;
  template class PETSc_KSPSolver<Vec, Mat, PETSc_SolverTypes::PETSc_KSPBICG, PETSc_Preconditioners::PETSc_DEFAULT>;
  template class PETSc_KSPSolver<Vec, Mat, PETSc_SolverTypes::PETSc_KSPBICG, PETSc_Preconditioners::PETSc_PCJACOBI>;
  template class PETSc_KSPSolver<Vec, Mat, PETSc_SolverTypes::PETSc_KSPGMRES, PETSc_Preconditioners::PETSc_DEFAULT>;
  template class PETSc_KSPSolver<Vec, Mat, PETSc_SolverTypes::PETSc_KSPGMRES, PETSc_Preconditioners::PETSc_PCJACOBI>;
  // ***************************************************************************
  template<typename PETSc_ArrayType,
           typename PETSc_SparseArrayType,
           PETSc_SolverTypes PETSc_SolverType,
           PETSc_Preconditioners PETSc_Preconditioner>
  PETSc_KSPSolver<
  PETSc_ArrayType,
  PETSc_SparseArrayType,
  PETSc_SolverType,
  PETSc_Preconditioner>::PETSc_KSPSolver()
  {
    _rightHandSide = nullptr;
    _solution = nullptr;
  }
  template<typename PETSc_ArrayType,
           typename PETSc_SparseArrayType,
           PETSc_SolverTypes PETSc_SolverType,
           PETSc_Preconditioners PETSc_Preconditioner>
  PETSc_KSPSolver<
  PETSc_ArrayType,
  PETSc_SparseArrayType,
  PETSc_SolverType,
  PETSc_Preconditioner>::~PETSc_KSPSolver()
  {
    _rightHandSide = nullptr;
    _solution = nullptr;
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType,
           typename PETSc_SparseArrayType,
           PETSc_SolverTypes PETSc_SolverType,
           PETSc_Preconditioners PETSc_Preconditioner>
  void PETSc_KSPSolver<
  PETSc_ArrayType,
  PETSc_SparseArrayType,
  PETSc_SolverType,
  PETSc_Preconditioner>::Initialize(const ISparseArray& matrix,
                                    const IArray& rightHandSide,
                                    IArray& solution,
                                    const Configuration& configuration)
  {
    _rightHandSide = &rightHandSide;
    _solution = &solution;

    Initialize(matrix,
               configuration);
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType,
           typename PETSc_SparseArrayType,
           PETSc_SolverTypes PETSc_SolverType,
           PETSc_Preconditioners PETSc_Preconditioner>
  ILinearSolver::SolutionInfo PETSc_KSPSolver<
  PETSc_ArrayType,
  PETSc_SparseArrayType,
  PETSc_SolverType,
  PETSc_Preconditioner>::Solve() const
  {
    if (_rightHandSide == nullptr ||
        _solution == nullptr)
      throw runtime_error("No initialization found");

    return Solve(*_rightHandSide,
                 *_solution);
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType,
           typename PETSc_SparseArrayType,
           PETSc_SolverTypes PETSc_SolverType,
           PETSc_Preconditioners PETSc_Preconditioner>
  void PETSc_KSPSolver<
  PETSc_ArrayType,
  PETSc_SparseArrayType,
  PETSc_SolverType,
  PETSc_Preconditioner>::Initialize(const ISparseArray& matrix,
                                    const Configuration& configuration)
  {
    _config = configuration;

    const auto& A_PETSc_Array = static_cast<const PETSc_SparseArray<PETSc_ArrayType, PETSc_SparseArrayType>&>(matrix);
    const PETSc_SparseArrayType& A_PETSc = static_cast<const PETSc_SparseArrayType&>(A_PETSc_Array);

    KSPCreate(PETSC_COMM_WORLD, &linearSolver);
    KSPSetOperators(linearSolver, A_PETSc, A_PETSc);

    switch (PETSc_SolverType)
    {
      case PETSc_SolverTypes::PETSc_KSPCG:
        KSPSetType(linearSolver, KSPCG);
        break;
      case PETSc_SolverTypes::PETSc_KSPBICG:
        KSPSetType(linearSolver, KSPBICG);
        break;
      case PETSc_SolverTypes::PETSc_KSPGMRES:
        KSPSetType(linearSolver, KSPGMRES);
        break;
      default:
        throw std::runtime_error("Unsopported PETSc solver type");
    }

    KSPGetPC(linearSolver, &preconditioner);
    switch (PETSc_Preconditioner)
    {
      case PETSc_Preconditioners::PETSc_DEFAULT:
        break;
      case PETSc_Preconditioners::PETSc_PCJACOBI:
        PCSetType(preconditioner, PCJACOBI);
        PCSetFromOptions(preconditioner);
        break;
      default:
        throw std::runtime_error("Unsopported PETSc preconditioner type");
    }

    KSPSetTolerances(linearSolver,
                     _config.Tolerance,
                     PETSC_DEFAULT,
                     PETSC_DEFAULT,
                     _config.MaxIterations);

    CHKERRABORT(PETSC_COMM_WORLD,
                KSPSetFromOptions(linearSolver));
  }
  // ***************************************************************************
  template<typename PETSc_ArrayType,
           typename PETSc_SparseArrayType,
           PETSc_SolverTypes PETSc_SolverType,
           PETSc_Preconditioners PETSc_Preconditioner>
  ILinearSolver::SolutionInfo PETSc_KSPSolver<
  PETSc_ArrayType,
  PETSc_SparseArrayType,
  PETSc_SolverType,
  PETSc_Preconditioner>::Solve(const IArray& rightHandSide,
                               IArray& solution) const
  {
    const auto& rhs_PETSc_Array = static_cast<const PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>&>(rightHandSide);
    const PETSc_ArrayType& rhs_PETSc = static_cast<const PETSc_ArrayType&>(rhs_PETSc_Array);

    const auto& sol_PETSc_Array = static_cast<const PETSc_Array<PETSc_ArrayType, PETSc_SparseArrayType>&>(solution);
    const PETSc_ArrayType& sol_PETSc = static_cast<const PETSc_ArrayType&>(sol_PETSc_Array);

    CHKERRABORT(PETSC_COMM_WORLD,
                KSPSolve(linearSolver, rhs_PETSc, sol_PETSc));

    KSPConvergedReason reason;
    KSPGetConvergedReason(linearSolver, &reason);
    if (reason < 0)
    {
      throw std::runtime_error("KSP did not converge.\n");
      PetscPrintf(PETSC_COMM_WORLD, "KSP did not converge.\n");
    }

    PetscInt iterations;
    KSPGetIterationNumber(linearSolver, &iterations);

    PetscReal residual;
    KSPGetResidualNorm(linearSolver, &residual);

    return
    {
      static_cast<unsigned int>(iterations),
      static_cast<double>(residual)
    };
  }
  // ***************************************************************************
}