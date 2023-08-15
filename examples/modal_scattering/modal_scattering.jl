using EffectiveTMatrix
using Plots
using LaTeXStrings
using Statistics

## define host_medium
dimension=2;
host_medium = Acoustic(dimension; ρ=1.0, c=1.0);

## modal source
input_mode = 3;
source = mode_source(host_medium,input_mode);

## set frequency and define where to plot 
ω = .1;
    M=N=40;
    res=100;
    bottomleft = [-M;-N]; topright = [M;N];
    region = Box([bottomleft, topright]);
    plot(source,ω;bounds=region,res=res)
    plot!(title="radial source")
# savefig("source_mode_"*string(input_mode)*".png")

## parameters for particles configurations
radius_big_cylinder = 20.0;
particle = Particle(Acoustic(2; ρ=Inf, c=Inf),Circle(1.0)); # sound hard particles
sp_MC = Specie(particle; volume_fraction = .05);

## draw one particles configuration and compute scattered field
particles_realisation = renew_particle_configurations(sp_MC,radius_big_cylinder);
plot(particles_realisation)
# savefig("particles_configuration.png")

sim = FrequencySimulation(particles_realisation,source);
basis_order=5;
scattered_field = run(sim,region,[ω];only_scattered_waves=true,basis_order=basis_order,res=res);
plot(scattered_field,ω; field_apply=real,seriestype = :contour,c=:balance) 
plot!(title="scattered field")
# savefig("scattering_mode_"*string(input_mode)*".png")

############ modes amplitudes ###############
basis_field_order = 4; # number of modes to compute
F = mode_analysis(input_mode, ω, host_medium, sp_MC;
                radius_big_cylinder=radius_big_cylinder, 
                basis_order=basis_order, 
                basis_field_order=basis_field_order,
                nb_iterations=1);

scatter(0:basis_field_order,abs.(F),label=false,markerstrokewidth=.5,markersize=7,markershape=:dtriangle)
scatter!(xlabel="n",ylabel=L"$|\mathfrak{F}_n|$")
scatter!(title="modes for one realisation")
# savefig("deterministic_modal_decomposition_mode_"*string(input_mode)*".png")


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
plot!(title="average scattered field")
# savefig("mean_us_mode_"*string(input_mode)*".png")

## Modes of the averaged scattered field
F_average = mode_analysis(input_mode, ω, host_medium, sp_MC;
                radius_big_cylinder=radius_big_cylinder, 
                basis_order=basis_order, 
                basis_field_order=basis_field_order,
                nb_iterations=nb_of_configurations);

scatter(0:basis_field_order,abs.(F_average),label=false,markerstrokewidth=.5,markersize=7,markershape=:dtriangle)
scatter!(xlabel="n",ylabel=L"$|\langle\mathfrak{F}_n\rangle|$")
scatter!(title="modes for $(nb_of_configurations) realisations")
# savefig("average_modal_decomposition_mode_"*string(input_mode)*".png")