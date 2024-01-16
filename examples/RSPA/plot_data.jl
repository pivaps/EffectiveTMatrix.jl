using EffectiveTMatrix
using Plots
using CSV
using Statistics
using DataFrames

############################ plots related to Monte Carlo validation ############################
# plots are stored in folders Data/1, Data/2, ..., Data/7
# as T0.pdf, ..., T4.pdf

root_folder = pwd() # This folder should contain Data
for p = 1:7
    data_path = root_folder * "/Data/"*string(p) *"/"
    MC_vec = MC_read(data_path * "MC.csv")
    for m = 0:4
        plot(MC_vec;mode=m)
        savefig(data_path * "T" * string(m) * ".pdf")
    end
end



############################ generates Table 1 and Table 2 ############################
# saved as table1_sound_hard_particles.csv and table2_sound_soft_particles.csv

root_folder = pwd() # This folder should contain Data

### Table 1 : Sound hard particles
errors_sound_hard = Array{Float64}(undef,3,5);
line = 0;
for p in [2,4,6]
    line+=1
    data_path = root_folder * "/Data/"*string(p) *"/"
    MC_vec = MC_read(data_path * "MC.csv")
    L = length(MC_vec)
    e = Array{Float64}(undef,L,5)
    for l=1:L
        x, _ = relative_error(MC_vec[l])
        e[l,:] = x;
    end
    mean_e = vec(mean(e,dims=1))
    errors_sound_hard[line,:] = mean_e
end

### Table 2 : Sound soft particles
errors_sound_soft = Array{Float64}(undef,3,5);
line = 0;
for p in [1,3,5]
    line+=1
    data_path = root_folder * "/Data/"*string(p) *"/"
    MC_vec = MC_read(data_path * "MC.csv")
    L = length(MC_vec)
    e = Array{Float64}(undef,L,5)
    for l=1:L
        x, _ = relative_error(MC_vec[l])
        e[l,:] = x;
    end
    mean_e = vec(mean(e,dims=1))
    errors_sound_soft[line,:] = mean_e
end     

volume_fractions = [0.05,0.1,0.2];  
df1 = DataFrame(hcat(volume_fractions,errors_sound_hard),:auto)
df2 = DataFrame(hcat(volume_fractions,errors_sound_soft),:auto)
rename!(df1,[:volume_fraction, :e0, :e1, :e2, :e3, :e4])     
rename!(df2,[:volume_fraction, :e0, :e1, :e2, :e3, :e4])  

CSV.write(root_folder * "/table1_sound_hard_particles.csv",df1)
CSV.write(root_folder * "/table2_sound_soft_particles.csv",df2)