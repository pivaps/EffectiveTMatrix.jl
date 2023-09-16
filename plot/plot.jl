@recipe function plot(MC::MonteCarloResult; field_apply=real,
    MC_color=:black, MA_color=:steelblue2, exact_color=:brown2)

    basis_field_order = length(MC.μ)-1
    η = uncertainty(MC)
    seriestype --> :scatter
    

    @series begin
        label --> "MC"
        yerror --> field_apply(η)
        color --> MC_color
        0:basis_field_order, field_apply.(MC.μ)
    end

    @series begin
        label --> "EWM"
        color --> exact_color
        0:basis_field_order, field_apply.(MC.μeff)
    end

    @series begin
        label --> "EWM-MA"
        color --> MA_color
        0:basis_field_order, field_apply.(MC.μeff0)
    end
end


@recipe function plot(MC_vec::Vector{MonteCarloResult{N}};
    mode=0, xproperty::Symbol=:ω, x_field_apply=real,
    MC_filter=1,MC_color=:black, MA_color=:steelblue2, exact_color=:brown2) where N
    
    basis_field_order = length(MC_vec[1].μ)
    if mode > basis_field_order
        @error("mode must range between 0 and basis_field_order = $basis_field_order")
    end
    η = [1.96*MC.σ[mode+1]/sqrt(MC.nb_iterations[mode+1]) for MC in MC_vec] 
    x_vec = x_field_apply([getproperty(MC,xproperty) for MC in MC_vec])
    layout := (1,2)
    
    ϕ = MC_vec[1].sp_MC.volume_fraction
    c = MC_vec[1].sp_MC.particle.medium.c
    ρ = MC_vec[1].sp_MC.particle.medium.ρ
    plot_title --> "ϕ = $ϕ, c=$c, ρ=$ρ\nmode N=$mode"

    # subplots attributes
    title --> ["\nReal\n" "\nImage\n"]
    xlabel --> "ka"
    ylabel --> ["T-matrix" ""]
    legend --> [true false]

    for (i,field_apply) in enumerate((real,imag))
        @series begin
            label --> "EWM"
            linewidth --> 3
            linecolor --> exact_color
            # linestyle --> :dash
            subplot := i
            x_vec,field_apply.([MC.μeff[mode+1] for MC in MC_vec])
        end

        @series begin
            label --> "EWM-MA"
            linewidth --> 3
            linecolor --> MA_color
            linestyle --> :dash
            subplot := i
            x_vec,field_apply.([MC.μeff0[mode+1] for MC in MC_vec])
        end

        @series begin
            label --> "MC"
            line --> false
            markershape --> :circle
            markersize --> 2.5
            markercolor --> MC_color
            markeralpha --> .7
            markerstrokecolor --> MC_color
            fillcolor := MC_color
            fillalpha := 0.3
            ribbon := field_apply(η)
            # yerror --> field_apply(η)
            subplot := i
            x_vec[1:MC_filter:end],field_apply.([MC.μ[mode+1] for MC in MC_vec[1:MC_filter:end]])
        end
    end
end



@recipe function plot(MC_vec_vec::Vector{Vector{MonteCarloResult{N}}};
    mode=0, xproperty::Symbol=:ω, field_apply=real, x_field_apply=real,
    subtitles=["\nDirichlet\n" "\nNeumann\n"],
    MC_filter=1,MC_color=:black, MA_color=:steelblue2, exact_color=:brown2) where N
    
    MC_vec = MC_vec_vec[1];
    basis_field_order = length(MC_vec[1].μ)
    if mode > basis_field_order
        @error("mode must range between 0 and basis_field_order = $basis_field_order")
    end

    layout := (1,2)
    
    ϕ = MC_vec[1].sp_MC.volume_fraction
    c = MC_vec[1].sp_MC.particle.medium.c
    ρ = MC_vec[1].sp_MC.particle.medium.ρ
    plot_title --> "ϕ = $ϕ, c=$c, ρ=$ρ\nmode N=$mode"

    # subplots attributes
    title --> subtitles
    xlabel --> "ka"
    ylabel --> ["T-matrix" ""]
    legend --> [true false]

    for i=1:2

        η = [1.96*MC.σ[mode+1]/sqrt(MC.nb_iterations[mode+1]) for MC in MC_vec_vec[i]] 
        x_vec = x_field_apply([getproperty(MC,xproperty) for MC in MC_vec_vec[i]])


        @series begin
            label --> "EWM"
            linewidth --> 3
            linecolor --> exact_color
            subplot := i
            x_vec,field_apply.([MC.μeff[mode+1] for MC in MC_vec_vec[i]])
        end

        @series begin
            label --> "EWM-MA"
            linewidth --> 3
            linecolor --> MA_color
            linestyle --> :dash
            subplot := i
            x_vec,field_apply.([MC.μeff0[mode+1] for MC in MC_vec_vec[i]])
        end

        @series begin
            label --> "MC"
            line --> false
            markershape --> :circle
            markersize --> 2.5
            markercolor --> MC_color
            markeralpha --> .7
            markerstrokecolor --> MC_color
            fillcolor := MC_color
            fillalpha := 0.3
            ribbon := field_apply(η)
            # yerror --> field_apply(η)
            subplot := i
            x_vec[1:MC_filter:end],field_apply.([MC.μ[mode+1] for MC in MC_vec_vec[i][1:MC_filter:end]])
            #x_vec[1:MC_filter:end],field_apply.([MC.μ[mode+1] for MC in MC_vec[1:MC_filter:end]])
    end
end
end

