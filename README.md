# EffectiveTMatrix.jl
A Julia library for computing the effective T-matrix of a random particulate sphere or cylinder.

The maths is briefly explained here through the examples. More details will be available soon in an upcoming publication. 

## I) Introduction: The mode to mode scattering in a 2D setting
### Modal source

For a 2D setting, a modal source is of the form 

$$\mathrm V_n(\mathbf{kr}) := \mathrm J_n(kr)\mathrm e^{\mathrm in\theta}$$

where $\mathbf{r}=(r\cos\theta,r\sin\theta)$. It is defined by 

```julia
dimension=2;                                       
host_medium = Acoustic(dimension; ρ=1.0, c=1.0);   # 2D acoustic problem
mode = 0;                                          # mode n=0 
source = mode_source(host_medium,mode);            # V_0
```

To plot the source field, we need to set a frequency $\omega$, define a box where to plot the field and choose a resolution

```julia
ω=0.1;                                    # frequency
M=N=20;                                   # sizes of the rectangle where to plot
bottomleft = [-M;-N]; topright = [M;N];
region = Box([bottomleft, topright]);
plot(source,ω;bounds=region,res=100)
```

<p align="center">
    <img
    src="examples/source_mode_0.png"
    alt="Alt text"
    title=""
    style="display: inline-block; margin: 0 auto; max-width: 200px">
</p>

### Scattering from one configuration of particles

We can draw one configuration of particles confined in a disc of choosen radius "radius_big_cylinder", then compute the scattered field as follows,

We first need to define the parameters for particles configurations:

```julia
# parameters for particles configurations
radius_big_cylinder = 20.0;                                 # radius of cylinder where particles are confined 
ϕ = 0.05;                                                   # volume fraction (density of particles)
separation_ratio = 1.001;                                   # minimum distance between particles
particle = Particle(Acoustic(2; ρ=Inf, c=Inf),Circle(1.0)); # particles type: sound hard particles
sp_MC, sp_EF = generate_species(radius_big_cylinder,
                                particle,
                                ϕ,
                                separation_ratio);          # defines specie sp_MC (sp_EF will be explained and used later)
```

sp_MC contains the statistics of the particles configuration, particles configurations are then drawn as follows:

```julia
particles_realisation = renew_particle_configurations(sp_MC,radius_big_cylinder);
plot(particles_realisation)
```

<p align="center">
    <img
    src="examples/particles_configuration.png"
    alt="Alt text"
    title="Particles configuration"
    style="display: inline-block; margin: 0 auto; max-width: 200px">
</p>



We can now compute the scattered field

```julia
basis_order=5;
sim = FrequencySimulation(particles_realisation,source);
scattered_field = run(sim,region,[ω];only_scattered_waves=true,basis_order=basis_order,res=res);

plot(scattered_field,ω; field_apply=real,seriestype = :contour,c=:balance) 
```

<p align="center">
    <img
    src="examples/scattering_mode_0.png"
    alt="Alt text"
    title=""
    style="display: inline-block; margin: 0 auto; max-width: 200px">
</p>

The Scattered field can be decomposed in modes $\mathrm U_n$ defined by 
$$\mathrm U_n(k\mathbf{r}) := \mathrm H_n(kr)\mathrm e^{\mathrm in\theta}$$

We have 

$$u_s(\mathbf{r}) =  \sum_{n=-\infty}^{+\infty} \mathfrak{F}_n\mathrm U_n(k\mathbf{r})$$

The modal amplitudes are computed as follows:

```julia
basis_field_order = 4;                                    # number of modes to compute
F1 = mode_analysis(mode, ω, host_medium, sp_MC;
                radius_big_cylinder=radius_big_cylinder, 
                basis_order=basis_order, 
                basis_field_order=basis_field_order,
                nb_iterations=1);

scatter(0:basis_field_order,abs.(F1),label=false,markerstrokewidth=.5,markersize=7,markershape=:dtriangle)
```

<p align="center">
    <img
    src="examples/deterministic_modal_decomposition_mode_0.png"
    alt="Alt text"
    title=""
    style="display: inline-block; margin: 0 auto; max-width: 200px">
</p>

### Average scattered field over particles configurations

The previous steps can be repeated in order to compute the averaged scattered field over several particles configurations, denoted $\langle u_s \rangle(\mathbf{r})$. 

$$\langle u_s \rangle(\mathbf{r}) =  \sum_{n=-\infty}^{+\infty} \langle \mathfrak{F}_n\rangle \mathrm U_n(k\mathbf{r})$$

The averaged scattered field is computed as follows:

```julia 
x_vec, _ = points_in_shape(region;resolution=res);                            # space discretization 
nb_of_configurations = 200;
A = complex(zeros(length(x_vec),nb_of_configurations));                       # store the fields of each configurations
@time Threads.@threads for i=1:nb_of_configurations
    particles = renew_particle_configurations(sp_MC,radius_big_cylinder);
    sim = FrequencySimulation(particles,source);
    us = run(sim,x_vec,[ω];only_scattered_waves=true,basis_order=basis_order)
    A[:,i] = mean.(us.field[:,1]) 
end 
mean_A  = mean(A,dims=2);
mean_us = FrequencySimulationResult(mean_A,x_vec,[ω]);
plot(mean_us,ω; field_apply=real,seriestype = :contour,c=:balance) 
```

<p align="center">
    <img
    src="examples/mean_us_mode_0.png"
    alt="Alt text"
    title=""
    style="display: inline-block; margin: 0 auto; max-width: 200px">
</p>


The empirical average of $\langle \mathfrak{F}_n\rangle$ obtained with 200 configurations is computed by

```julia 
F2 = mode_analysis(mode, ω, host_medium, sp_MC;
                radius_big_cylinder=radius_big_cylinder, 
                basis_order=basis_order, 
                basis_field_order=basis_field_order,
                nb_iterations=nb_of_configurations);

scatter(0:basis_field_order,abs.(F2),label=false,markerstrokewidth=.5,markersize=7,markershape=:dtriangle)
scatter!(xlabel="n",ylabel=L"$\langle\mathfrak{F}_n\rangle$")
scatter!(title="modes for one $(nb_of_configurations) realisations")
```


<p align="center">
    <img
    src="examples/average_modal_decomposition_mode_0.png"
    alt="Alt text"
    title=""
    style="display: inline-block; margin: 0 auto; max-width: 200px">
</p>



## II) The effective T-matrix 

After ensemble averaging over particle configurations within a cylindrical region the system inherits cylindrical symmetry. 
For example, if the source has radial symmetry then the average scattered field will also have radial symmetry. This is also true for sources 
with more general rotational symmetry resulting in scattered fields of the same rotational symmetry.

If we respectively denote by $\mathrm V_n$ and $\mathrm U_n$ the incident and scattered modes of order $n$, then they are linearly related by some 
complex number $\mathrm T_n$ as follows:

$$
   \mathrm  U_n = \mathrm T_n \mathrm V_n,\quad \mathrm T_n\in \mathbb{C},
$$

where  only depends on the properties of the random material, such as it's radius and constitutive particles. Since the scattering 
problem is linear, the knowledge of the $\mathrm T_n$ allows to describe the scattering from any incident field, after decomposing the latter into 
the modes $\mathrm V_n$.

This library computes $\mathrm T_n$, provided the particles types and the material radius by using a 
We do so by using the effective wave method approach \cite{gower2021effective}. Finally, we show that it is possible to
improve the rate of convergence of Monte Carlo simulations thanks to the cylindrical symmetry. We use the latter to validate our results.