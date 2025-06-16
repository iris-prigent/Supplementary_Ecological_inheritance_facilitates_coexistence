using Random
using StatsBase
using Distributions
using DelimitedFiles

### parameter values ###
b1::Float64=1.0
b2::Float64=2.0
c1::Float64=1.07
c2::Float64=1.5
b3::Float64=-2.5
xnp::Int=50 #number of the patches on the x axis
ynp::Int=50 #number of patches on the y axis
n::Int = 10 #patch size
p_mut::Float64=0.005 #mutation probability
sigma_mut::Float64=0.0025 #standard deviation for mutation size
gen::Int=100000 #number of generations to be iterated
sample_freq::Int=1000 #number of generations between each sampling of the entire population
gamma::Int=1
delta::Float64=1
f0::Float64=1
lambda::Float64=0.7
m::Float64=0.5


### computes of the maximal value of z (eq. 3) ####
zmax::Float64=0
if b2>0
zmax=-b1/(2*b3);
else
zmax=-b1/(2*(b2+b3));
end
############


### initialization ###
iniz::Float64=0.1 #sets the initial value of z
pop=fill(iniz,(n,xnp,ynp)) #creates a collection of xnp by ynp individuals carrying the trait z, each located on the x and y axis
e0::Float64=(lambda*n*iniz*(b1+(b2+b3)*iniz))/(1-lambda) #computes the ecological equilibrium (eq. C-12)
envt=fill(e0,(xnp,ynp)) #creates a collection of  xnp by ynp patches, each at the ecological equilibrium

### prepares the output directory and files ####
dirname="lattice_model_fixed_patch_size" 
path=string(dirname,"_c1_",c1,"_c2_",c2,"_b1_",b1,"_b2_" , b2 ,"_b3_" , b3,"_gamma_",gamma ,"_lambda_",lambda ,"_m_",m) #full directory name (with parameter values for clarity)
mkpath(string(path,"/sampling")) #creates directories for the output 


text0=string(0," ",iniz, " ",0, " ", e0," ", 0) #fills a vector with the intial values for: generation (i.e. iteration of the life-cycle, starting at 0); average trait value; variance in trait value; average environmental quality, variance in environmental quality)
open(string(path,"/variables.txt"), "w") do io
    println(io, text0) #writes this vector in a file
end
############


##### saves the parameters used in the current simulation run in a file ####
namepara=string("b1 b2 b3 c1 c2 gamma lambda m f0 delta p_mut sigma_mut n_indiv n_patches_x_axis n_patches_y_axis n_generations sample_freq")
allpara=string(b1," ",b2," ",b3," ",c1," ",c2," ",gamma," ",lambda," ",m," ",f0," ",delta," ",p_mut," ",sigma_mut," ",n," ",xnp," ",ynp," ",gen," ",sample_freq)
open(string(path,"/parameters.txt"), "w") do io
    println(io, namepara)
    println(io, allpara)
end
#######################

#### defines a function computing the amount of environmental modifiations from the investment of individuals within a patch (from a vector of n traits) ####
function investment(b1::Float64,b2::Float64,b3::Float64,indiv_patch::Vector{Float64})
    B::Float64=b1*sum(indiv_patch)+(b2-b3/(n-1))*sum(indiv_patch.*indiv_patch)+b3/(n-1)*sum(indiv_patch*indiv_patch')
    return B
end

#### defines a function computing the fecundity of all individuals in one patch ####
function fecundity(c1::Float64,c2::Float64,gamma::Int,delta::Float64,f0::Float64,envt_patch::Float64,indiv_patch::Vector{Float64},B::Float64)
    fec::Vector{Float64}=((-c1.*indiv_patch-c2.*indiv_patch.*indiv_patch.+gamma.*B.+envt_patch).*delta.+f0)
    return fec
end

#### defines a function weighting  each individual in each neighbouring patch according to their fecundity and the probability that their offspring migrates into the patch considered, eq.E41 ####
function migration(fec_9,m)
    p::Array{Float64}=fec_9*(m/8)
    p[:,2,2]=fec_9[:,2,2]*(1-m)
    return p
end

#declares the vector collecting the fecundities of all individuals in the population, for ease of computation
fec=fill(0.0,(n,xnp,ynp))

@time for it in 1:gen
    popnext=fill(0.0,(n,xnp,ynp)) #declares the vector collecting the trait of the next generation
    for i in 1:xnp, j in 1:ynp #patch per patch
        patch_pop=pop[:,i,j]
        envt_pop=envt[i,j]
        B::Float64=investment(b1,b2,b3,patch_pop) # computes the individuals' contribution to their environment
        fec[:,i,j]=fecundity(c1,c2,gamma,delta,f0,envt_pop,patch_pop,B)  #computes the individuals' fecundity
        envt[i,j]=(envt_pop+B) #updates the environment quality according to the current generation's contribution
    end

    for i in 1:xnp, j in 1:ynp #patch per patch
        xi=[i-1,i,i+1] #selects patches at less than 1 distance (i.e. patch) apart in the x direction
        yi=[j-1,j,j+1] #selects patches at less than 1 distance (i.e. patch) apart in the y direction
        xi[xi.<1].=xnp #corrects for the first patches in the x axis
        xi[xi.>xnp].=1 #corrects for the last patches in the x axis
        yi[yi.<1].=ynp #corrects for the first patches in the y axis
        yi[yi.>ynp].=1 #corrects for the last patches in the y axis
        p=migration(fec[:,xi,yi],m) #wcreates a vector weighting the parental generation by their fecundity and the probability that their offspring migrates into the patch considered
        popnext[:,i,j] = sample(vec(pop[:,xi,yi]), Weights(vec(p)), n) #samples the successful offspring from the weighted parental generation
    end
    id=sample(1:(xnp*ynp*n), rand(Binomial(xnp*ynp*n,p_mut)),replace=false)  #samples which individual mutates
    mutants=popnext[id].+rand(Normal(0,sigma_mut), length(id)) #updates the trait value of the mutants, sampling the mutation from a normal distribution
    mutants[mutants.<0].=0 #truncates the mutant trait values below 0
    mutants[mutants.>zmax].=zmax  #truncates the mutant trait values above zmax
    popnext[id]=mutants #updates the trait value in the overall population
    global pop=popnext  #updates the parent generation (adults die and the offspring replace them)
     meanz::Float64=sum(pop)/(xnp*ynp*n) #computes the average trait value in the population
     varz::Float64=sum(vec(pop).*vec(pop))/(xnp*ynp*n)-meanz^2 #computes the variance in trait value in the population
     meane::Float64=sum(envt)/(xnp*ynp) #computes the average environment quality
     vare::Float64=sum(envt.*envt)/(xnp*ynp)-meane^2 #computes the variance in environment quality
    text=string(it," ", meanz, " ",varz, " ", meane," ", vare)
    open(string(path,"/variables.txt"), "a") do io  #writes the averages and variances in the output files
        println(io, text)
    end
    if it%sample_freq==0#each sample_freq generations
        #samples the entire population (trait values and environment quality)
        println(it)
        writedlm(string(path,"/sampling/e_generation_",it,".txt"), vec(envt))
        writedlm(string(path,"/sampling/z_generation_",it,".txt"),  vec(pop))
    end
    global envt=envt*lambda #updates the environment quality (after degradation of the public good)
end
