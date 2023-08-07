using EffectiveTMatrix
using Plots

## define host_medium
dimension=2;
host_medium = Acoustic(dimension; ρ=1.0, c=1.0);

## modal source
mode = 0;
source = mode_source(host_medium,mode);

## set frequency and define where to plot 
ω = .1;
M=N=40;
res=100;
bottomleft = [-M;-N]; topright = [M;N];
region = Box([bottomleft, topright]);
plot(source,ω;bounds=region,res=res)
# savefig("source_mode_"*string(mode)*".png")

## parameters for particles configurations
radius_big_cylinder = 20.0;
ϕ = 0.05; # volume fraction (density of particles)
separation_ratio = 1.001; # minimum distance between particles
particle = Particle(Acoustic(2; ρ=Inf, c=Inf),Circle(1.0)); # sound hard particles
sp_MC, sp_EF = generate_species(radius_big_cylinder,particle,ϕ,separation_ratio);

## draw one particles configuration and compute scattered field
particles_realisation = renew_particle_configurations(sp_MC,radius_big_cylinder);
plot(particles_realisation)
# savefig("particles_configuration.png")

sim = FrequencySimulation(particles_realisation,source);
basis_order=5;
scattered_field = run(sim,region,[ω];only_scattered_waves=true,basis_order=basis_order,res=res);
plot(scattered_field,ω; field_apply=real,seriestype = :contour,c=:balance) 
# savefig("scattering_mode_"*string(mode)*".png")

############ modes amplitudes ###############
basis_field_order = 4; # number of modes to compute
F1 = mode_analysis(mode, ω, host_medium, sp_MC;
                radius_big_cylinder=radius_big_cylinder, 
                basis_order=basis_order, 
                basis_field_order=basis_field_order,
                nb_iterations=1);

scatter(0:basis_field_order,abs.(F1),label="M=1",markerstrokewidth=.5,markersize=7,markershape=:dtriangle)
