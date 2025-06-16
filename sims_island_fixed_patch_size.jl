using Random
using StatsBase
using Distributions
using DelimitedFiles

### parameter values ####
b1::Float64=1.0 
b2::Float64=1.6
b3::Float64=-2.5
c1::Float64=0.96
c2::Float64=1.5
np::Int=2000 #number of patches
n::Int=10 #patch size
p_mut::Float64=0.005 #mutation probability
sigma_mut::Float64=0.0025 #standard deviation for mutation size
gen::Int=100000 #number of generations to be iterated
sample_freq::Int=1000 #number of generations between each sampling of the entire population
gamma::Float64=1
delta::Float64=1
f0::Float64=1
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
iniz=0.1 #sets the initial value of z
pop=fill(iniz,(n,np)) #creates a collection of n individuals carrying the trait z, each located in one of the np patch
e0::Float64=(lambda*n*iniz*(b1+(b2+b3)*iniz))/(1-lambda) #computes the ecological equilibrium (eq. C-12)
envt=fill(e0,np) #creates a collection of np patches, each at the ecological equilibrium


### prepares the output directory and files ####
dirname="island_model_fixed_patch_size" 
path=string(dirname,"_c1_",c1,"_c2_",c2,"_b1_",b1,"_b2_" , b2 ,"_b3_" , b3,"_gamma_",gamma ,"_lambda_",lambda ,"_m_",m) #full directory name (with parameter values for clarity)
mkpath(string(path,"/sampling")) #creates directories for the output 

text0=string(0," ",iniz, " ",0, " ", e0," ", 0) #fills a vector with the intial values for: generation (i.e. iteration of the life-cycle, starting at 0); average trait value; variance in trait value; average environmental quality, variance in environmental quality
open(string(path,"/variables.txt"), "w") do io
    println(io, text0) #writes this vector in a file
end
############


##### saves the parameters used in the current simulation run in a file ####
namepara=string("b1 b2 b3 c1 c2 gamma lambda m f0 delta p_mut sigma_mut n_indiv n_patches n_generations sample_freq")
allpara=string(b1," ",b2," ",b3," ",c1," ",c2," ",gamma," ",lambda," ",m," ",f0," ",delta," ",p_mut," ",sigma_mut," ",n," ",np," ",gen," ",sample_freq)
open(string(path,"/parameters.txt"), "w") do io
    println(io, namepara)
    println(io, allpara)
end
############


#### defines a function computing the amount of environmental modifiations from the investment of individuals within a patch (from a vector of traits) ####
function investment(b1,b2,b3,indiv_patch)
    B::Float64=b1*sum(indiv_patch)+(b2-b3/(n-1))*sum(indiv_patch.*indiv_patch)+b3/(n-1)*sum(indiv_patch*indiv_patch') #from eq. 2 / eq. C28
return B
end

#### defines a function computing the fecundity of all individuals in the population ####
function fecundity(c1,c2,gamma,delta,f0,envt,pop,B)
    f::Matrix{Float64}=((-c1.*pop-c2.*pop.*pop.+gamma.*B'.+envt').*delta.+f0) #from eq. 4 / eq. C29
    return f
end

#### computes the probabilities of moving patch and remaining in the focal patch ####
p_disp::Float64=(m/(np-1))
p_remain::Float64=(1-m)

#### declares vectors for ease of computation
B=fill(0.0,np) ### vector collecting the contribution (B(z)) in each patches

p=fill(0.0,(n,np)) ### vector collecting the probability of being sampled as a parent (collecting the elements computed using eq. C30,for each k)

@time  for it in 1:gen
    popnext=fill(0.0,(n,np))#declares the vector collecting the trait of the next generation
    for i in 1:np #patch per patch, compute the individuals' contribution to their environment
        B[i]=investment(b1,b2,b3,pop[:,i])
    end
    fec_all=fecundity(c1,c2,gamma,delta,f0,envt,pop,B)#computes all individuals' fecundity
    p[:,1]=fec_all[:,1].*p_remain #weights the fecundities by the probability of having an offspring competing in the patch being filled (first, patch 1)
    p[:,2:np]=fec_all[:,2:np].*p_disp
    popnext[:,1] = sample(vec(pop), Weights(vec(p)), n) #samples the next generation
    for i in 2:np #patch per patch
        p[:,i-1]=fec_all[:,i-1].*p_disp #updates the vector weighting each individual
        p[:,i]=fec_all[:,i].*p_remain
        popnext[:,i] = sample(vec(pop), Weights(vec(p)), n) #and samples the next generation
    end
    id=sample(1:(np*n), rand(Binomial(n*np,p_mut)),replace=false) #samples which individual mutates
    mutants::Vector{Float64}=popnext[id].+rand(Normal(0,sigma_mut), length(id)) #updates the trait value of the mutants, sampling the mutation from a normal distribution
    mutants[mutants.<0].=0 #truncates the mutant trait values below 0
    mutants[mutants.>zmax].=zmax #truncates the mutant trait value above zmax
    popnext[id]=mutants #updates the trait value in the overall population
    global pop=popnext #updates the parent generation (adults die and the offspring replace them)
    global envt=(envt.+B) #updates the environment quality (before degradation of the public good)

    meanz::Float64=sum(pop)/(np*n) #computes the average trait value in the population
    varz::Float64=sum(pop.*pop)/(np*n)-meanz^2 #computes the variance in trait value in the population
    meane::Float64=sum(envt)/(np) #computes the average environment quality
    vare::Float64=sum(envt.*envt)/(np)-meane^2 #computes the variance in environment quality
    
    text=string(it," ", meanz, " ",varz, " ", meane," ", vare)
    open(string(path,"/variables.txt"), "a") do io #writes the averages and variances in the output files
      println(io, text)
  end

  if it%sample_freq==0 #each sample_freq generations
    #sampless the entire population (trait values and environment quality)
      println(it)
      writedlm(string(path,"/sampling/e_generation_",it,".txt"), vec(envt))
      writedlm(string(path,"/sampling/z_generation_",it,".txt"), vec(pop))
    end
  end
  global envt=envt.*lambda #updates the environment quality (after degradation of the public good)
end

