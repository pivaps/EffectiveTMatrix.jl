using Test, EffectiveTMatrix

## define host_medium
dimension=2;
host_medium = Acoustic(dimension; ρ=1.0, c=1.0);

## modal source
mode = 0;
source = mode_source(host_medium,mode);

## parameters for particles configurations
radius_big_cylinder = 20.0;
ϕ = 0.05; # volume fraction (density of particles)
separation_ratio = 1.001; # minimum distance between particles
particle = Particle(Acoustic(2; ρ=Inf, c=Inf),Circle(1.0)); # sound hard particles
sp_MC, sp_EF = generate_species(radius_big_cylinder,particle,ϕ,separation_ratio);

## Modes of the averaged scattered field
ω = .1;
basis_order = 5;
basis_field_order = 4;
nb_of_configurations = 5000;
F2 = mode_analysis(mode, ω, host_medium, sp_MC;
                radius_big_cylinder=radius_big_cylinder, 
                basis_order=basis_order, 
                basis_field_order=basis_field_order,
                nb_iterations=nb_of_configurations);

@test abs(F2[1]) > 10*abs(F2[2])
@test abs(F2[1]) > 10*abs(F2[3])
@test abs(F2[1]) > 10*abs(F2[4])
@test abs(F2[1]) > 10*abs(F2[5])