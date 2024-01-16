using EffectiveTMatrix

## define host_medium, particles (sound soft or sound hard) and big cylinder radius:
dimension=2;
host_medium = Acoustic(dimension; ρ=1.0, c=1.0);
sound_soft_particle = Particle(Acoustic(2; ρ=1e-2, c=1.0+0im),Circle(1.0));   
sound_hard_particle = Particle(Acoustic(2; ρ=Inf, c=Inf),Circle(1.0)); 
cylinder_radius = 20.0; # big cylinder radius


# Choose range of modes where to compute Tₙ
basis_field_order = 4; # sample Tₙ for n ∈ [0;basis_field_order] 

## Monte Carlo keywords
kws_MC = Dict(
    :basis_order => 10,                       # truncation of the Foldy-Lax equation 
    :nb_iterations_max => 3000,              # max number of iterations
    :nb_iterations_step => 200,              # convergence criteria is checked at every steps
    :prec => 5e-2                            # precision of the Monte Carlo simulation
);

############################ DATA 1 ############################
# frequency loop - sound soft particles - volume_fraction = 0.05
# -> In the paper: Figure 2 (left), Figure 8, Table 2 (row 1). 

# list of frequencies
Ω = collect(0.05:0.015:1.5);

# particulate microstructure of the effective cylinder
particle = sound_soft_particle
sp_MC=Specie(particle; volume_fraction = .05)

# Initialize data
MC_vec = [MonteCarloResult(ω,sp_MC,cylinder_radius;
     basis_field_order=basis_field_order) for ω ∈ Ω];

# Update MC_vec with Monte Carlo simulations
sample_effective_t_matrix!(MC_vec, host_medium; kws_MC...);

# Update MC_vec with the Effective Waves Method
kws_MC[:basis_order]=4 # cannot be too high for sound soft particles
run_MC_validation!(host_medium, MC_vec; kws_MC...);

# save data
folder = pwd() # Choose preferred path for saving data
comment = "frequency loop - sound soft particles - volume_fraction = 0.05";
save(MC_vec,comment,folder)

############################ DATA 2 ############################
# frequency loop - sound hard particles - volume_fraction = 0.05 
# -> In the paper: Figure 2 (right), Table 1 (row 1).

# list of frequencies
Ω = collect(0.05:0.015:1.5);

# particulate microstructure of the effective cylinder
particle = sound_hard_particle
sp_MC=Specie(particle; volume_fraction = .05)

# Initialize data
MC_vec = [MonteCarloResult(ω,sp_MC,cylinder_radius;
     basis_field_order=basis_field_order) for ω ∈ Ω];

# Update MC_vec with Monte Carlo simulations
sample_effective_t_matrix!(MC_vec, host_medium; kws_MC...);

# Update MC_vec with the Effective Waves Method
kws_MC[:basis_order]=10
run_MC_validation!(host_medium, MC_vec; kws_MC...);

# save data
folder = pwd() # Choose preferred path for saving data
comment = "frequency loop - sound hard particles - volume_fraction = 0.05";
save(MC_vec,comment,folder)


############################ DATA 3 ############################
# frequency loop - sound soft particles - volume_fraction = 0.1
# -> In the paper: Figure 9 (top left), Table 2 (row 2). 

# list of frequencies
Ω = collect(.1:.1:1.5);

# particulate microstructure of the effective cylinder
particle = sound_soft_particle
sp_MC=Specie(particle; volume_fraction = .1)

# Initialize data
MC_vec = [MonteCarloResult(ω,sp_MC,cylinder_radius;
     basis_field_order=basis_field_order) for ω ∈ Ω];

# Update MC_vec with Monte Carlo simulations
sample_effective_t_matrix!(MC_vec, host_medium; kws_MC...);

# Update MC_vec with the Effective Waves Method
kws_MC[:basis_order]=4 # cannot be too high for sound soft particles
run_MC_validation!(host_medium, MC_vec; kws_MC...);

# save data
folder = pwd() # Choose preferred path for saving data
comment = "frequency loop - sound soft particles - volume_fraction = 0.1";
save(MC_vec,comment,folder)


############################ DATA 4 ############################
# frequency loop - sound hard particles - volume_fraction = 0.1 
# -> In the paper: Figure 9 (top right), Table 1 (row 2). 

# list of frequencies
Ω = collect(0.05:0.015:1.5);

# particulate microstructure of the effective cylinder
particle = sound_hard_particle
sp_MC=Specie(particle; volume_fraction = .1)

# Initialize data
MC_vec = [MonteCarloResult(ω,sp_MC,cylinder_radius;
     basis_field_order=basis_field_order) for ω ∈ Ω];

# Update MC_vec with Monte Carlo simulations
sample_effective_t_matrix!(MC_vec, host_medium; kws_MC...);

# Update MC_vec with the Effective Waves Method
kws_MC[:basis_order]=10
run_MC_validation!(host_medium, MC_vec; kws_MC...);

# save data
folder = pwd() # Choose preferred path for saving data
comment = "frequency loop - sound hard particles - volume_fraction = 0.1";
save(MC_vec,comment,folder)


############################ DATA 5 ############################
# frequency loop - sound soft particles - volume_fraction = 0.2
# -> In the paper: Figure 9 (bottom left), Table 2 (row 3). 

# list of frequencies
Ω = collect(.1:.1:1.5);

# particulate microstructure of the effective cylinder
particle = sound_soft_particle
sp_MC=Specie(particle; volume_fraction = .2)

# Initialize data
MC_vec = [MonteCarloResult(ω,sp_MC,cylinder_radius;
     basis_field_order=basis_field_order) for ω ∈ Ω];

# Update MC_vec with Monte Carlo simulations
sample_effective_t_matrix!(MC_vec, host_medium; kws_MC...);

# Update MC_vec with the Effective Waves Method
kws_MC[:basis_order]=4 # cannot be too high for sound soft particles
run_MC_validation!(host_medium, MC_vec; kws_MC...);

# save data
folder = pwd() # Choose preferred path for saving data
comment = "frequency loop - sound soft particles - volume_fraction = 0.2";
save(MC_vec,comment,folder)


############################ DATA 6 ############################
# frequency loop - sound hard particles - volume_fraction = 0.2
# -> In the paper: Figure 9 (bottom right), Table 1 (row 3). 

# list of frequencies
Ω = collect(0.05:0.015:1.5);

# particulate microstructure of the effective cylinder
particle = sound_hard_particle
sp_MC=Specie(particle; volume_fraction = .2)

# Initialize data
MC_vec = [MonteCarloResult(ω,sp_MC,cylinder_radius;
     basis_field_order=basis_field_order) for ω ∈ Ω];

# Update MC_vec with Monte Carlo simulations
sample_effective_t_matrix!(MC_vec, host_medium; kws_MC...);

# Update MC_vec with the Effective Waves Method
kws_MC[:basis_order]=10
run_MC_validation!(host_medium, MC_vec; kws_MC...);

# save data
folder = pwd() # Choose preferred path for saving data
comment = "frequency loop - sound hard particles - volume_fraction = 0.2";
save(MC_vec,comment,folder)



############################ DATA 7 ############################
# frequency loop - monopole particles - volume_fraction = 0.1
# -> In the paper: Figure 10

# list of frequencies
Ω = collect(0.05:0.015:1.5);

# particulate microstructure of the effective cylinder
particle = sound_hard_particle
kws_MC[:basis_order]=0 # basis_order = 0 corresponds to monopole scattering
sp_MC=Specie(particle; volume_fraction = .2)

# Initialize data
MC_vec = [MonteCarloResult(ω,sp_MC,cylinder_radius;
     basis_field_order=basis_field_order) for ω ∈ Ω];

# Update MC_vec with Monte Carlo simulations
sample_effective_t_matrix!(MC_vec, host_medium; kws_MC...);

# Update MC_vec with the Effective Waves Method
run_MC_validation!(host_medium, MC_vec; kws_MC...);

# save data
folder = pwd() # Choose preferred path for saving data
comment = "frequency loop - monopole particles - volume_fraction = 0.1";
save(MC_vec,comment,folder)
