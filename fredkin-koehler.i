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

[UserObjects]
  [./bifmm]
    type = BoundaryIntegralFMM
    variable = phi2
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
