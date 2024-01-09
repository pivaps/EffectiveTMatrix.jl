using EffectiveTMatrix
using Plots
using LaTeXStrings

## Validation for fixed set of parameters
## define host_medium
dimension=2;
host_medium = Acoustic(dimension; ρ=1.0, c=1.0);

## particulate microstructure of the effective cylinder

particle = Particle(Acoustic(2; ρ=1e-2, c=1.0),Circle(1.0));   # sound soft particle of radius 1.0 
sp_MC = Specie(particle; volume_fraction = .05);

## radius of the effective cylinder
cylinder_radius = 20.0;

basis_field_order = 4; # sample Tₙ for n ∈ [0;basis_field_order] 
ω=.1

## Monte Carlo keywords
kws_MC = Dict(
    :basis_order => 5,                       # truncation of the Foldy-Lax equation 
    :nb_iterations_max => 6000,              # max number of iterations
    :nb_iterations_step => 100,              # convergence criteria is checked at every steps
    :prec => 1e-2                            # precision of the Monte Carlo simulation
);

# Initialize MC
MC = MonteCarloResult(ω,sp_MC,cylinder_radius; basis_field_order=basis_field_order);

# Update MC with the Monte Carlo simulation
sample_effective_t_matrix!(MC, host_medium; kws_MC...);

# Update MC with the effective method
run_MC_validation!(host_medium, MC; kws_MC...);

plot(MC)
plot!(xlabel="n",ylabel=L"$\mathrm{R\!e}\,(\mathrm{T}_n)$")
savefig("MC_fixed_params.png")



## Initialize MC_vec
Ω = collect(.1:.05:1.5);
MC_vec = [MonteCarloResult(ω,sp_MC,cylinder_radius; basis_field_order=basis_field_order) for ω ∈ Ω];

# Update MC_vec with Monte Carlo simulations
sample_effective_t_matrix!(MC_vec, host_medium; kws_MC...);

# Update MC_vec with the effective method
run_MC_validation!(host_medium, MC_vec; kws_MC...);


# save data
folder = "/home"; # choose where to save Data
comment = "simulation description";
save(MC_vec,comment,folder)

# plot results
plot(MC_vec;MC_color=:black,MA_color=:steelblue1,exact_color=:brown1)
# savefig("MC_multi_frq.png")



