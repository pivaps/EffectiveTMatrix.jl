using EffectiveTMatrix
using Plots
using LaTeXStrings

## define host_medium
dimension=2;
host_medium = Acoustic(dimension; ρ=1.0, c=1.0);

## particulate microstructure of the effective cylinder
particle = Particle(Acoustic(2; ρ=Inf, c=Inf),Circle(1.0));   # sound hard particles of radius 1.0 
sp_MC = Specie(particle; volume_fraction = .4) 
microstructure = Microstructure(host_medium,[sp_MC]);

## define the effective cylinder
cylinder_radius = 20.0;
cylinder = Material(Circle(cylinder_radius),microstructure);

## compute the coefficients of the T-matrix 
ω = .2; # frequency
N = 5;  # Tₙ for n ∈ [-N,N]
T = t_matrix(ω, host_medium, cylinder, basis_field_order=N);

## plot Tₙ for n ∈ [0,N]
scatter(0:N,real.(T[N+1:2N+1]),label=L"$\mathrm{Re\,(T}_n)$")
scatter!(0:N,imag.(T[N+1:2N+1]),label=L"$\mathrm{Im\,(T}_n)$")
scatter!(title="Effective T-matrix diagonal terms", xlabel=L"$n$", ylabel=L"$\mathrm{T}_n$")
# savefig("T-matrix.png")

## Define incident planar wave and point source
# define the plane wave
plane_wave = plane_source(host_medium; direction = [1.0,0.0]);

# define the point source
x0 = [-4.5cylinder_radius,0.0];                    # origin of the point source
A = 1.0;                                           # amplitude of the point source
point_wave = point_source(host_medium, x0, A);


# define box where to compute the average scattered field 
P=Q=4*cylinder_radius;                   # bounding box size
bottomleft = [-P;-Q]; topright = [P;Q];
region = Box([bottomleft, topright]);    # bounding box

## scattering from incident plane wave
us_plane_wave = average_scattered_field(ω, region, plane_wave, cylinder; basis_field_order=N);
plot(us_plane_wave, ω; field_apply=real,seriestype = :contour,c=:balance) 
plot!(title="Scattering from incident plave wave")
# savefig("plane_wave_scattering.png")

## scattering from point source
us_point_source = average_scattered_field(ω, region, point_wave, cylinder; basis_field_order=N);
plot(us_point_source,ω; field_apply=real,seriestype = :contour,c=:balance) 
plot!(title="Scattering from incident point source")
# savefig("point_source_scattering.png")

