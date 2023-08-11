using EffectiveTMatrix
using Plots
using LaTeXStrings

## define host_medium
dimension=2;
host_medium = Acoustic(dimension; ρ=1.0, c=1.0);

## particulate microstructure of the effective cylinder
particle = Particle(Acoustic(2; ρ=Inf, c=Inf),Circle(1.0));   # sound hard particles of radius 1.0 
ϕ = 0.4;                                                      # density of particles
sp_MC = Specie(particle; volume_fraction = ϕ) 
micro = Microstructure(host_medium,[sp_MC]);

## define the effective cylinder
cylinder_radius = 20.0;
cylinder = Material(Circle(cylinder_radius),micro);

## compute the coefficients of the T-matrix 
ω = .2; # frequency
N = 5;  # Tₙ for n ∈ [-N,N]
T = t_matrix(ω, host_medium, cylinder, basis_field_order=N);

## plot Tₙ for n ∈ [0,N]
scatter(0:N,real.(T[N+1:2N+1]),label=L"$\mathrm{Re\,(T}_n)$")
scatter!(0:N,imag.(T[N+1:2N+1]),label=L"$\mathrm{Im\,(T}_n)$")
scatter!(title="Effective T-matrix diagonal terms", xlabel=L"$n$", ylabel=L"$\mathrm{T}_n$")
savefig("T-matrix.png")

## average scattered field from incident plane wave
# define the plane wave
psource = plane_source(host_medium; direction = [1.0,0.0]);

# define box where to compute the average scattered field 
P=Q=4*cylinder_radius;                   # bounding box size
bottomleft = [-P;-Q]; topright = [P;Q];
region = Box([bottomleft, topright]);    # bounding box

us = average_scattered_field(ω, region, psource, cylinder; basis_field_order=N);
plot(us,ω; field_apply=real,seriestype = :contour,c=:balance) 
plot!(title="average pressure field")
savefig("plane_wave_scattering.png")



