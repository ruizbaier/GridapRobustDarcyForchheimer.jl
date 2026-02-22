module DarcyForchheimerConvergence3D_mixedBCs
  using GridapRobustDarcyForchheimer
  using Gridap
  import Gridap: ∇
  using Printf
  using Test
  using LineSearches: BackTracking
  using LinearAlgebra

  # Arbitrary model parameters
  const κ = 1.0 # permeability
  const r = 3.0     # Forchheimer exponent
  const FF = 1.0  # Forchheimer coefficient

  function normu(v)
    (v⋅v>0)*(sqrt(v⋅v))^(r-2)
  end

  function deri(v) #derivative of |v| wrt v
    (v⋅v>0)*(r-2)*(sqrt(v⋅v))^(r-4)*v
  end

  u_ex(x)  = VectorValue(cos(π*x[1])*sin(π*x[2])*sin(π*x[3]),
                          -sin(π*x[1])*cos(π*x[2])*sin(π*x[3]),
                          sin(π*x[1])*sin(π*x[2])*cos(π*x[3]))

  # or a div-free velocity 
  # u_ex(x)  = VectorValue(sin(π*x[1])*cos(π*x[2])*cos(π*x[3]),
  #                       -2*cos(π*x[1])*sin(π*x[2])*cos(π*x[3]),
  #                       cos(π*x[1])*cos(π*x[2])*sin(π*x[3]))                        

  p_ex(x)  = sin(π*x[1])*cos(π*x[2])*cos(π*x[3])

  f_ex(x) = (1.0/κ)*u_ex(x) + ∇(p_ex)(x) + FF*u_ex(x)*(normu∘u_ex)(x) 
  g_ex(x) = (∇⋅u_ex)(x)

  function solve_Forchheimer(model; k=k, generate_output=false)
    reffe_u = ReferenceFE(raviart_thomas,Float64,k)
    reffe_p = ReferenceFE(lagrangian,Float64,k)

    # FESpaces
    Uh_ = TestFESpace(model,reffe_u,dirichlet_tags="Gamma_u",conformity=:HDiv)
    Qh_ = TestFESpace(model,reffe_p,conformity=:L2)

    Uh = TrialFESpace(Uh_,u_ex)
    Qh = TrialFESpace(Qh_)

    Yh = MultiFieldFESpace([Uh_,Qh_])
    Xh = MultiFieldFESpace([Uh,Qh])

    Ω = Triangulation(model)
    dΩ = Measure(Ω,2*(k+2))
    Γp = BoundaryTriangulation(model, tags = "Gamma_p")
    n_Γp = get_normal_vector(Γp)
    dΓp = Measure(Γp,2*(k+2)-1)

    a(u,v)   = ∫( (1/κ)*(u⋅v))*dΩ
    c(u,w,v) = ∫(FF*v⋅((normu∘(u))*w))dΩ
    d1(u,w,v) = ∫(FF*v⋅((normu∘(u))*w))dΩ
    d2(u,w,v) = ∫(FF*(u⋅v)*((deri∘(u))⋅w))dΩ 
    b(v,q)   = ∫(-(∇⋅v)*q)*dΩ
    F(v)     = ∫(f_ex⋅v)*dΩ -∫(p_ex*(v⋅n_Γp))dΓp
    G(q)     = ∫((-1.0)*(g_ex*q))dΩ

    resid((u,p),(v,q)) = a(u,v) + c(u,u,v) + b(v,p) + b(u,q) - F(v) - G(q) 

    jacob((u,p),(du,dp),(v,q)) = a(du,v) + d1(u,du,v) + d2(u,du,v) + b(v,dp) + b(du,q)  

    oper = FEOperator(resid,jacob,Xh,Yh)
    nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(), iterations=25, ftol=1.0E-7)
    solver = FESolver(nls)
    x0=0.0001*ones(num_free_dofs(Xh)) 
    xh=FEFunction(Xh,x0)
    (uh,ph), _ = solve!(xh,solver,oper)

    pa(p,q) = ∫(p*q)*dΩ
    pdivl(q) = ∫(q*(∇⋅uh - g_ex))*dΩ
    projection = AffineFEOperator(pa,pdivl,Qh,Qh)
    loss_div = solve(projection)
    loss_divh=norm(get_free_dof_values(loss_div),Inf)

    euh = u_ex - uh
    eph = p_ex - ph
    eu = (sum(∫((euh⋅euh).^(0.5*r))dΩ))^(1.0/r) + sqrt(sum(∫((∇⋅euh)*(∇⋅euh))*dΩ))
    ep = sqrt(sum(∫(eph*eph)*dΩ))

    if generate_output 
      writevtk(Ω,"DFconvergence3D_k0_ncells=$(num_cells(model))",order=1,cellfields=["uh"=>uh,"ph"=>ph])
   end

   eu, ep, loss_divh, Gridap.FESpaces.num_free_dofs(Xh)

  end


  function  convergence_test(; nkmax, k=k, generate_output=false)
    eu = Float64[]
    ep = Float64[]
    ru = Float64[]
    rp = Float64[]
    hh = Float64[]
    losses_divh = Float64[]
    nn = Int[]
    newton_iters = Int[]

    push!(ru,0.)
    push!(rp,0.)

    for nk in 1:nkmax
        println("******** Refinement step: $nk")
        model=generate_model3d(nk)
        setup_unit_cube_labels!(model)
        error_u,error_p,loss_divh,ndofs=solve_Forchheimer(model; k=k, generate_output=generate_output)

        push!(nn,ndofs)
        push!(hh,sqrt(3/2)/(2^(nk)))
        #push!(hh,sqrt(2)/(2^nk)) # for 2D
        println("******** Total DoFs: ", nn[nk])
        
        push!(eu,error_u)
        push!(ep,error_p)
        push!(losses_divh,loss_divh)

        if nk>1
            push!(ru, log(eu[nk]/eu[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rp, log(ep[nk]/ep[nk-1])/log(hh[nk]/hh[nk-1]))
        end
    end
    println("=======================================================")
    println("   DoF  &    h   &   e(u)   &  r(u) &   e(p)   &  r(p) & loss_divh ")
    println("=======================================================")
    for nk in 1:nkmax
        @printf("%7d & %.4f & %.2e & %.3f & %.2e & %.3f & %.2e\n", 
                nn[nk], hh[nk], eu[nk], ru[nk], ep[nk], rp[nk], losses_divh[nk]);
    end
    println("=======================================================")
  end
  convergence_test(;nkmax=4,k=0,generate_output=false)
end

#output table for k=0

#    144 & 0.6124 & 9.78e-01 & 0.000 & 1.80e-01 & 0.000 & 1.40e-15  4
#   1152 & 0.3062 & 5.30e-01 & 0.884 & 9.61e-02 & 0.907 & 2.24e-15  5
#   9216 & 0.1531 & 2.72e-01 & 0.965 & 4.88e-02 & 0.977 & 6.53e-15  5
#  73728 & 0.0765 & 1.37e-01 & 0.991 & 2.45e-02 & 0.994 & 1.59e-14  5
# 589824 & 0.0383 & 6.85e-02 & 0.998 & 1.23e-02 & 0.999 & 1.64e-13  6
# 4096512 & 0.0191 & 3.43e-02 & 0.999 & 6.15e-03 & 0.999 & 1.31e-13 6

# and for k=1
#     624 & 0.6124 & 3.53e-01 & 0.000 & 6.32e-02 & 0.000 & 1.48e-14 5
#   4992 & 0.3062 & 9.86e-02 & 1.840 & 1.73e-02 & 1.871 & 3.33e-14  6
#  39936 & 0.1531 & 2.55e-02 & 1.950 & 4.42e-03 & 1.968 & 9.14e-14  6
# 319488 & 0.0765 & 6.45e-03 & 1.986 & 1.11e-03 & 1.992 & 2.25e-13  6
# 2555904 & 0.0383 & 1.62e-03 & 1.993 & 2.78e-04 & 1.996 & 1.01e-13 7 
# 20447232 & 0.0191 & 4.07e-04 & 1.997 & 6.95e-05 & 1.998 & 2.01e-13 7