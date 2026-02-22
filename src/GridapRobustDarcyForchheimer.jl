module GridapRobustDarcyForchheimer

    using Gridap
    using Gridap.ReferenceFEs
    using Gridap.FESpaces
    using LinearOperators
    using LinearAlgebra
    using Krylov
    using Printf

    include("Preconditioners.jl")
    include("DarcyForchheimerTools.jl")

    export assemble_perturbed_darcy, assemble_darcy, assemble_darcy_noExact, compute_errors_darcy, compute_errors_perturbed_darcy, compute_operator_norm
    export generate_model2d, generate_model3d, generate_rectangle, setup_model_labels_rectangle_channel!, setup_unit_cube_labels!, setup_model_labels_unit_square!, computation_error_darcy_minimization
    export GridapLinearSolverPreconditioner, custom_stopping_condition
end
