using Random
using StatsBase
using Distributions
using DelimitedFiles

### parameter values ####
b1::Float64=1.0
b2::Float64=1.6
c1::Float64=0.96
c2::Float64=1.5
b3::Float64=-2.5
np::Int=2000 #number of patches
n::Int=10 #patch size at generation 0 (patch size later fluctuates)
p_mut::Float64=0.005 #mutation probability
sigma_mut::Float64=0.0025 #standard deviation for mutation size
gen::Int=100000 #number of generations to be iterated
sample_freq::Int=1000 #number of generations between each sampling of the entire population
gamma::Int=1
f0::Float64=1.05
delta::Float64=1*f0
chi::Float64=0.09;
lambda::Float64=0.7
m::Float64=0.5
############


### computes of the maximal value of z (eq. 3) ####
zmax::Float64=0
if b2>0
zmax=-b1/(2*b3);
else
zmax=-b1/(2*(b2+b3));
end
############


### initialization ###
iniz::Float64=0.1; #sets the initial value of z
pop=fill(fill(iniz,n),np) #creates a collection of n individuals carrying the trait z, each located in one of the np patch
e0::Float64=(lambda*n*iniz*(b1+(b2+b3)*iniz))/(1-lambda) #computes the ecological equilibrium (eq. C-12)
envt=fill(e0,np) #creates a collection of np patches, each at the ecological equilibrium


### prepares the output directory and files ####
dirname="island_model_fluctuating_patch_size" 
path=string(dirname,"_c1_",c1,"_c2_",c2,"_b1_",b1,"_b2_" , b2 ,"_b3_" , b3,"_gamma_",gamma ,"_lambda_",lambda ,"_m_",m) #full directory name (with parameter values for clarity)
mkpath(string(path,"/sampling")) #creates directories for the output 

text0=string(0," ",iniz, " ",0, " ", e0," ", 0, " ",n*np) #fills a vector with the intial values for: generation (i.e. iteration of the life-cycle, starting at 0); average trait value; variance in trait value; average environmental quality, variance in environmental quality; total population size
open(string(path,"/variables.txt"), "w") do io
    println(io, text0) #writes this vector in a file
end
############

##### saves the parameters used in the current simulation run in a file ####
namepara=string("b1 b2 b3 c1 c2 gamma lambda m f0 delta chi p_mut sigma_mut n_indiv n_patches n_generations sample_freq")
allpara=string(b1," ",b2," ",b3," ",c1," ",c2," ",gamma," ",lambda," ",m," ",f0," ",delta," ",chi," ",p_mut," ",sigma_mut," ",n," ",np," ",gen," ",sample_freq)
open(string(path,"/parameters.txt"), "w") do io
    println(io, namepara)
    println(io, allpara)
end
###############

#### defines a function computing the amount of environmental modifiations from the investment of individuals within a patch (from a vector of traits) ####
function investment(b1,b2,b3,indiv_patch)
    n_indiv::Int=length(indiv_patch)
    B::Float64=0 #if there is no individual in the patch, the total contribution is 0
    if n_indiv==1 
        B=b1*indiv_patch[1]+b2*indiv_patch[1]*indiv_patch[1]#if there is one individual, between-individuals interaction do not affect the environment
    elseif n_indiv>1
        B=b1*sum(indiv_patch)+(b2-b3/(n_indiv-1))*sum(indiv_patch.*indiv_patch)+b3/(n_indiv-1)*sum(indiv_patch*indiv_patch')  #else, B is computed from from eq. 2 / eq. D31
    end
    return B
end

### defines a function sampling the number of offspring produced from a Poisson distribution with parameter equal to individual fecundity (computed using eq. 4 / eq. C29)
function reproduction(f0,delta,gamma,B_patch,envt_patch,c1,c2,indiv)
    noff::Int=rand(Poisson(f0+delta*(gamma*B_patch+envt_patch-c1*indiv-c2*indiv*indiv))) #
    return noff
end

### defines a function that samples which offspring survive in a given patch, from the vector of offspring competing within that patch
function competition(popof_deme,chi)
    noff::Int=length(popof_deme)
    new_patch::Vector{Float64}=sample(popof_deme, rand(Binomial(noff,(1/(1+chi*noff)))),replace=false)
    return new_patch
end

B=fill(0.0,np) #### declares a vector collecting the contribution (B(z)) in each patches for ease of computation

#### computes the probabilities of for an offspring to move into each patch
p=fill(m/(np-1),np)
p[1]=1-m

@time for it in 1:gen
    popof= Vector{Float64}[] #declares the vector collecting the trait of all offspring produced by parents (after reproduction and migration, before competition)
    for i in 1:np
        push!(popof, Vector{Float64}[]) #specifies that this is a vector of vector (vector of individuals within each patch)
    end
    popnext= Vector{Float64}[] #declares the vector collecting the trait of the next generation (surviving offspring, after competiting)
    for i in 1:np #in each patch
        B[i]=investment(b1,b2,b3,pop[i]) #computes the change in environment quality
        babies=fill(0.0,0) #declares a vector collecting the trait of all offspring (after reproduction, before migration)
        for i2 in 1:length(pop[i]) #for each individual in the patch, samples the number of offspring produced
            n_off=reproduction(f0,delta,gamma,B[i],envt[i],c1,c2,pop[i][i2])
            append!(babies,fill(pop[i][i2],n_off)) #appends the trait of each offspring produced to the vector collecting the trait of all offspring
        end
        id_disp=sample(1:np, Weights(p), length(babies)) #samples the identity of the patch that each offspring competes in
        for ideme in 1:length(babies)
            push!(popof[id_disp[ideme]],babies[ideme]) #to each patch, append the offspring that will compete in in
        end
        #updates the probability to be assigned to each patch, according to which patch is considered next
        p[i]=m/(np-1)
        if i<np
            p[i+1]=1-m
        else
            p[1]=1-m
        end
    end


    for i in 1:np #in each patch
        push!(popnext,competition(popof[i],chi)) #samples which offspring survives
        id=sample(1:length(popnext[i]), rand(Binomial(length(popnext[i]),p_mut)),replace=false) #samples which of these surviving offspring mutates
        mutants=popnext[i][id].+rand(Normal(0,sigma_mut), length(id)) #updates the trait value of the mutants, sampling the mutation from a normal distribution
        mutants[mutants.<0].=0 #truncates the mutant trait values below 0
        mutants[mutants.>zmax].=zmax #truncates the mutant trait value above zmax
        popnext[i][id]=mutants #updates the trait value in the overall population
        envt[i]=(envt[i]+B[i]) #updates the environment quality (before degradation of the public good)
    end

    global pop=popnext  #updates the parent generation (adults die and the offspring replace them)

    vec_pop=collect(Iterators.flatten(pop)) #flattens the vector of patch made of vector of traits into a large vector of traits
    meanz::Float64=mean(vec_pop ) #computes the average trait value in the population
    varz::Float64=mean(vec_pop.*vec_pop)-meanz^2 #computes the variance in trait value in the population
    meane=mean(envt) #computes the average environment quality
    vare=mean(envt.*envt)-meane^2 #computes the variance in environment quality
    

    n_all_pop=length(vec_pop) #computes the total population size
    text=string(it," ", meanz, " ",varz, " ", meane," ", vare, " ", n_all_pop)
    open(string(path,"/variables.txt"), "a") do io #writes the averages, variances and total, population size in the output files
        println(io, text)
    end

    if it%sample_freq==0 #each sample_freq generations
        #samples the entire population (trait values, environment quality, patch size)
        println(it)
        writedlm(string(path,"/sampling/e_generation_",it,".txt"), envt)
        writedlm(string(path,"/sampling/z_generation_",it,".txt"), vec_pop)
        nindiv=fill(1,np) #creates a vector that collects all the patch sizes
        for i in 1:np 
            deme::Vector{Float64}=pop[i] #computes patch size
            nindiv[i]=length(deme)
        end
        writedlm(string(path,"/sampling/i_generation_",it,".txt"), nindiv)
    end
  end
  global envt=envt.*lambda #updates the environment quality (after degradation of the public good)
end

