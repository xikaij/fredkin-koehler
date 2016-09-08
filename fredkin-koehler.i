[Mesh]
  file = sphere.e
  dim = 3
[]

[Variables]
  [./phi2]
    order = FIRST
    family = LAGRANGE
    initial_condition = 1.0
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = phi2
  [../]
[]

[UserObjects]
  [./bifmm]
    type = BoundaryIntegralFMM
    variable = phi2
    cx = 0.0
    cy = 0.0
    cz = 0.0
    boxWidth = 2.1
    TreeHeight = 4
    execute_on = timestep_begin
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
[]

[Outputs]
  output_initial = true
  exodus = true
  print_perf_log = true
[]
