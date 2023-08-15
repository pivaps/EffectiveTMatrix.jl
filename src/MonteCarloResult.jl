# The object MonteCarloResult contains the parameters and result of the the MC simulations
# N is basis_field_order+1
 mutable struct MonteCarloResult{N} 
    basis_order::Int64
    ω::Float64
    sp_MC::Specie
    R::Float64
    μ::SVector{N,ComplexF64} 
    σ::SVector{N,ComplexF64}
    nb_iterations::SVector{N,Int64}
    μeff::SVector{N,ComplexF64}
    μeff0::SVector{N,ComplexF64}
end

function MonteCarloResult(ω::Float64,sp_MC::Specie,R::Float64; basis_order=5::Int, basis_field_order=0::Int)
    N = basis_field_order+1
    a = Array{Int,1}(undef,N)
    b = Array{ComplexF64,1}(undef,N)
    return MonteCarloResult{N}(basis_order, ω, sp_MC, R, b, b, a, b, b)
end

function uncertainty(V::Vector{Float64})
    return 1.96*std(V)/sqrt(length(V))    
end

function uncertainty(V::Vector{ComplexF64})
    return uncertainty(real.(V))+im*uncertainty(imag.(V))
end

function uncertainty(MC::MonteCarloResult)
    return 1.96*MC.σ./sqrt.(MC.nb_iterations)
end

function relative_error(MC::MonteCarloResult)
    return abs.(MC.μ-MC.μeff)./abs.(MC.μ), abs.(MC.μ-MC.μeff0)./abs.(MC.μ)    
end

function CSV_string_vec(v::Vector)
    string_v = string()
    for x in v
        string_v *= "$x,"
    end
    return "\""*chop(string_v)*"\""
end

function MC_write(MC_vec::Vector{MonteCarloResult},file_path::String) 
    basis_field_order = length(MC_vec[1].μ)-1
    header = [
        "basis_order",
        "basis_field_order",

        "ω",
        "R",
        
        "particle_radius",
        "c_particle",
        "ρ_particle",
        "separation_ratio",
        "volume_fraction",

        "μ",
        "σ",
        "nb_iterations",

        "μeff",
        "μeff0"
        ]

        open(file_path, "w") do f
            H=string()
            for h in header
                H*=h*","
            end 
            write(f, chop(H))

            write(f,"\n")
            for MC in MC_vec
                write(f,
                    "$(MC.basis_order),",
                    "$(basis_field_order),",

                    "$(MC.ω),",
                    "$(MC.R),",

                    "$(outer_radius(MC.sp_MC.particle)),",
                    "$(MC.sp_MC.particle.medium.c),",
                    "$(MC.sp_MC.particle.medium.ρ),",
                    "$(MC.sp_MC.separation_ratio),",
                    "$(MC.sp_MC.volume_fraction),"
                )

                write(f, CSV_string_vec(MC.μ) ,",")
                write(f, CSV_string_vec(MC.σ) ,",")
                write(f, CSV_string_vec(MC.nb_iterations) ,",")

                write(f, CSV_string_vec(MC.μeff) ,",")
                write(f, CSV_string_vec(MC.μeff0))

                write(f,"\n")
            end
            write(f,"\n")
        end
end


function MC_read(file_path::String)
    # file_path = "/home/kevish/Documents/Numeric/Julia/CylindersProject/Final/Data/Temp/4/MC.csv"
    file = CSV.File(file_path; types=Dict(:c_particle => ComplexF64)) 
    MC_vec=Vector{MonteCarloResult}()
    for i = 1:length(file)
        particle = Particle(Acoustic(2; ρ=file.ρ_particle[i], c=file.c_particle[i]),Circle(file.particle_radius[i]))
        sp_MC = Specie(particle; volume_fraction = file.volume_fraction[i],separation_ratio=file.separation_ratio[i]) 
        
        μ = parse.(ComplexF64,split(file.μ[i],","))
        σ = parse.(ComplexF64,split(file.σ[i],","))
        nb_iterations = parse.(Int64,split(file.nb_iterations[i],","))

        μeff = parse.(ComplexF64,split(file.μeff[i],","))
        μeff0 = parse.(ComplexF64,split(file.μeff0[i],","))


        push!(MC_vec,
                MonteCarloResult(file.basis_order[i],file.ω[i],sp_MC,file.R[i],μ,σ,nb_iterations,μeff,μeff0)
                )
    end
    return MC_vec
end

function MC_read(data_folder::Int=1)  MC_read(pwd()*"/Data/"*string(data_folder)*"/MC.csv") end

## commented below: function contains part of code to write realisations F
# function write_realisations(MC_vec::Vector{MonteCarloResultTemp},file_path::String)
#         open(file_path, "w") do f
#             # write header F0, F1 ...
#             FH = string()
#             for N = 0:basis_field_order
#                 FH*="F"*string(N)*","
#             end
#             write(f, chop(FH))

#             write(f,"\n")
#             for MC in MC_vec
#                 write(f,
#                     "$(MC.basis_order),",
#                     "$(MC.basis_field_order),",

#                     "$(MC.ω),",
#                     "$(MC.R),",

#                     "$(outer_radius(MC.sp_MC.particle)),",
#                     "$(MC.sp_MC.particle.medium.c),",
#                     "$(MC.sp_MC.particle.medium.ρ),",
#                     "$(MC.sp_MC.separation_ratio),",
#                     "$(MC.sp_MC.volume_fraction),"
#                 )

#                 write(f, CSV_string_vec(MC.μeff) ,",")
#                 write(f, CSV_string_vec(MC.μeff0),",")
                
#                 for N = 1:basis_field_order
#                     write(f, CSV_string_vec(MC.F[N]) ,",")
#                 end
#                 write(f, CSV_string_vec(MC.F[basis_field_order+1]))

#                 write(f,"\n")
#             end
#             write(f,"\n")
#         end
# end
# -----------------------------------------------------------------------------------------

## Commented below, the function contains part of the code to read realisations F
# function read_realisation(data_folder::Int=1)
#     file_path = pwd()*"/Data/"*string(data_folder)*"/realisations/MCtemp.csv"
#     file = CSV.File(file_path; types=Dict(:c_particle => ComplexF64)) 
#     MCtemp_vec=Vector{MonteCarloResultTemp}()
#     for i = 1:length(file)
#         particle = Particle(Acoustic(2; ρ=file.ρ_particle[i], c=file.c_particle[i]),Circle(file.particle_radius[i]))
#         sp_MC = Specie(particle; volume_fraction = file.volume_fraction[i],separation_ratio=file.separation_ratio[i]) 

#         μeff = parse.(ComplexF64,split(file.μeff[i],","))
#         μeff0 = parse.(ComplexF64,split(file.μeff0[i],","))

#         F = [ComplexF64[] for _ = 0:file.basis_field_order[i]]
#         for N=0:file.basis_field_order[i]
#             F[N+1] = parse.(ComplexF64,split(getproperty(file,Symbol("F"*string(N)))[i],","))
#         end

#         push!(MCtemp_vec,
#                 MonteCarloResultTemp(file.basis_order[i],file.basis_field_order[i],file.ω[i],sp_MC,file.R[i],F,μeff,μeff0)
#                 )
#     end
#     return MCtemp_vec
# end
# -----------------------------------------------------------------------------------------

# commented below: function which contains part of the code to save realisations F
# function save(MCtemp_vec::Vector{MonteCarloResultTemp},description::String,all_data_path::String=pwd())
#     new_data_path = save(MonteCarloResult(MCtemp_vec),description,all_data_path)
#     realisations_path = new_data_path*"/realisations"
#     mkdir(realisations_path)
#     MCtemp_write(MCtemp_vec,realisations_path*"/MCtemp.csv")
# end

function load_parameters(data_folder::Int=1)
    file_path = pwd()*"/Data/"*string(data_folder)*"/MC.csv"
    file = CSV.File(file_path; types=Dict(:c_particle => ComplexF64)) 
    i = 1 # parameters of first element of MC_vec in the file

    particle = Particle(Acoustic(2; ρ=file.ρ_particle[i], c=file.c_particle[i]),Circle(file.particle_radius[i]))
    sp_MC = Specie(particle; volume_fraction = file.volume_fraction[i],separation_ratio=file.separation_ratio[i]) 
        
    ω = file.ω[i]
    radius_big_cylinder = file.R[i]
    basis_order = file.basis_order[i]
    basis_field_order = file.basis_field_order[i]
    sp_EF = sp_MC_to_EF(sp_MC,radius_big_cylinder)
    
    return ω, radius_big_cylinder, basis_order, basis_field_order, sp_MC, sp_EF
end


# please provide folder path which CONTAINS the folder named "Data" inside
function save(MC_vec::Vector{MonteCarloResult},description::String,all_data_path::String=pwd())

    # Check all_data_path is a valid path 
    if !isdir(all_data_path)
        all_data_path = pwd()
        @warn("path "* all_data_path * "does not exist and is replaced with "* pwd())
    end

    # chop last "/" in pathname if any
    if all_data_path[end] == '/'
        all_data_path=chop(all_data_path)
    end

    # check if all_data_path contains folder Data
    # create Data folder otherwise.
    if !isdir(all_data_path*"/Data")
    # folder Data not found
        if all_data_path[end-4:end]!="/Data" 
        # path doesn't include Data neither
            mkdir(all_data_path*"/Data")
            @warn("Data Folder created at "* all_data_path)
        else 
        # path includes Data
            all_data_path = chop(all_data_path;tail=5)
        end
    end

    # update metadata file if it exists
    metadata_path = all_data_path*"/Data/metadata.csv"
    if isfile(metadata_path)
    # metadata file exists: update it
        global num_current_data
        open(metadata_path) do f
            num_total_data = -1 
            for _ in enumerate(eachline(f))
                num_total_data += 1
            end
            global num_current_data = num_total_data+1
        end
        open(metadata_path,"a") do f
            write(f,"$num_current_data,",description,",",string(Dates.now()),"\n")
        end
    else
    # metadata file doesn't exist: create, put headers, and new data info
        open(metadata_path,"w") do f
            num_current_data = 1
            write(f,"Folder,Description,Date\n")
            write(f,"$num_current_data,",description,",",string(Dates.now()),"\n")
        end
    end    

    # save new data
    new_data_path = all_data_path*"/Data/"*string(num_current_data)
    if isdir(new_data_path)
        @error("Folder $new_data_path already exist.")
    end
    mkdir(new_data_path)
    MC_write(MC_vec,new_data_path*"/MC.csv")
    return new_data_path
end



# delete saved data and update metadatafile
function delete()

end

function Base.show(io::IO,MC::MonteCarloResult)
    basis_field_order = length(MC.μ)-1
    stdm_r = real.(MC.σ)./ sqrt.(MC.nb_iterations)
    stdm_i = imag.(MC.σ)./ sqrt.(MC.nb_iterations)

    print(io,"----------------------------------------------------------------------------------------------------")
    print(io, 
        "
        ω = $(MC.ω)
        container radius: $(MC.R)
        Particle type: $(MC.sp_MC.particle) 
        Values:
     ")

     for N = 0:basis_field_order
        print(io, 
        "
        F_$(N) =  $(MC.μ[N+1]) ± $(stdm_r[N+1]+im*stdm_i[N+1]) ($(MC.nb_iterations[N+1]) iterations) 
        ")
     end

     print(io,"----------------------------------------------------------------------------------------------------")
end