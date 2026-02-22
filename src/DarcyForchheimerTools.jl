# Computation of error using the minimization sub-problem
function computation_error_darcy_minimization(xh, op, yh, opmin, dΩ, dΩ1, κ, u_ex, p_ex, p_min_ex)
    uh, ph, _ = xh
    pminh, _ = yh
 
    error_u = sqrt(sum(∫(1.0/κ*(u_ex-uh)⋅(u_ex-uh))*dΩ + ∫((∇⋅(u_ex-uh))*(∇⋅(u_ex-uh)))*dΩ))
    error_p = sqrt(sum(∫((p_ex-ph-pminh)*(p_ex-ph-pminh))*dΩ + ∫(κ*((∇(pminh))⋅(∇(pminh))))*dΩ))
    error_pmin = sqrt(sum(∫((p_min_ex-pminh)*(p_min_ex-pminh))*dΩ1 + ∫(κ*(∇(p_min_ex-pminh))⋅(∇(p_min_ex-pminh)))*dΩ1))
    return error_u, error_p, error_pmin
end

function generate_model2d(nk)
  domain =(0,1,0,1)
  n      = 2^nk #+ 1
  # Discrete model
  partition = (n,n)
  CartesianDiscreteModel(domain, partition) |> simplexify
end

function generate_model3d(nk)
  domain =(0,1,0,1,0,1)
  n      = 2^nk 
  # Discrete model
  partition = (n,n,n)
  CartesianDiscreteModel(domain, partition) |> simplexify
end

function generate_rectangle(nk)
  domain =(0,2,0,1)
  n      = 2^nk + 1
  # Discrete model
  partition = (2*n,n)
  CartesianDiscreteModel(domain, partition) |> simplexify
end

function setup_model_labels_unit_square!(model)
  labels = get_face_labeling(model)
  add_tag!(labels,"Gamma_p",[6,2,3,4,8]) #6top, right, and other 3 corners
  add_tag!(labels,"Gamma_u",[1,2,3,5,7]) # bottom and 7left and low-left corner
end 

#
#       3 --- 10 --- 4
#     19 |         20 |
#   /   13       /    14
# 7 --- 12 --- 8      |
# |      |     |      |
#15      1 ----| 9 -- 2
# |    17      16   18              ^ y+
# |  /         |  /                 |
# 5 --- 11 --- 6                    |       x+
#                                   / ------>
#                                  /
#        3 ---------- 4           v z+
#      / |          / |          
#   /    | 24    /    |  <-- 21 (back)
# 7 ---------- 8      |
# |      |     |   26 |
# |  25  1 ----|----- 2
# |    /  22   |    /
# |  /         |  /
# 5 ---------- 6
#       ^
#       |
#      23 (bottom)


function setup_unit_cube_labels!(model)
    labels = get_face_labeling(model)
    #add_tag!(labels,"Gamma",[1,2,4,5,6,8,21,23,26]) # for P2 and P1 only
    add_tag!(labels,"Gamma_u",[1,2,4,5,6,8,21,24,26,9,10,11,13,14,16,17,18,20])
    add_tag!(labels,"Gamma_p",[22,23,25])
  end
  
function setup_model_labels_rectangle_channel!(model)
  labels = get_face_labeling(model)
  add_tag!(labels,"Gamma_out",[8]) #right 
  add_tag!(labels,"Gamma_in",[7,1,3]) # left and 2 corners
  add_tag!(labels,"Gamma_topbot",[1,2,3,4,5,6]) # corners bottom and top
end  






function compute_operator_norm(B,xh)
  x=get_free_dof_values(xh)
  sqrt(dot(x,B*x))
end





