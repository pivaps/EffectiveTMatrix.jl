using EffectiveTMatrix
using Plots
using LaTeXStrings

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

scatter(0:basis_field_order,abs.(F1),label=false,markerstrokewidth=.5,markersize=7,markershape=:dtriangle)
scatter!(xlabel="n",ylabel=L"$\mathfrak{F}_n$")
scatter!(title="modes for one realisation")
# savefig("deterministic_modal_decomposition_mode_"*string(mode)*".png")


## Average field
x_vec, _ = points_in_shape(region;resolution=res); # space discretization 
nb_of_configurations = 200;
A = complex(zeros(length(x_vec),nb_of_configurations)); # store the fields of each configurations
@time Threads.@threads for i=1:nb_of_configurations
    particles = renew_particle_configurations(sp_MC,radius_big_cylinder);
    sim = FrequencySimulation(particles,source);
    us = run(sim,x_vec,[ω];only_scattered_waves=true,basis_order=basis_order)
    A[:,i] = mean.(us.field[:,1]) # trick to get rid of SVector
end 
mean_A  = mean(A,dims=2);
mean_us = FrequencySimulationResult(mean_A,x_vec,[ω]);
plot(mean_us,ω; field_apply=real,seriestype = :contour,c=:balance) 
# savefig("mean_us_mode_"*string(mode)*".png")

## Modes of the averaged scattered field
F2 = mode_analysis(mode, ω, host_medium, sp_MC;
                radius_big_cylinder=radius_big_cylinder, 
                basis_order=basis_order, 
                basis_field_order=basis_field_order,
                nb_iterations=nb_of_configurations);

scatter(0:basis_field_order,abs.(F2),label=false,markerstrokewidth=.5,markersize=7,markershape=:dtriangle)
scatter!(xlabel="n",ylabel=L"$\langle\mathfrak{F}_n\rangle$")
scatter!(title="modes for one $(nb_of_configurations) realisations")
# savefig("average_modal_decomposition_mode_"*string(mode)*".png")