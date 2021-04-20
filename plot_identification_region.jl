using Plots; pyplot()
println("Caution. The 3d plot works under Julia v.1.4.1 Plots v1.6.10")
using Polyhedra
using Makie
C = 1
F = 5
B=range(0.05,stop=0.5,length=100)
E=range(0.05,stop=0.5,length=100)
lb(B,E) = C/B +0*E
ub(B,E) = C/(-E) +0*B
f3(B,E) = 0#B + E
f4(B,E) = F
Plots.plot(camera=(-10,20),
           xlabel = "B",
           ylabel = "E",
           zlabel = "β",
           legend = :topleft,
           title = "UB, LB of separated (red, blue), LB of composite (green)")
Plots.plot!(B,E,lb,st=:surface,
           label = "Upper bound of β",
           color=:red,
           alpha = 0.1)
Plots.plot!(B,E,ub,st=:surface,
           label = "Lower bound of β",
           #color=:blue,
           alpha = 0.1)
Plots.plot!(B,E,f3,st=:surface,
           label = "Lower bound of β",
           color=:green,
           alpha = 0.1)
Plots.plot!(B,E,f4,st=:surface,
           label = "Upper bound of β",
           color=:orange,
           alpha = 0.1)
savefig("julia_merger_figure/identification_region_3d_positive_B_positive_E")


C = 1
F = 5
B=range(0.05,stop=0.5,length=100)
E=range(-0.05,stop=-0.5,length=100)
lb(B,E) = C/B +0*E
ub(B,E) = C/(-E) +0*B
f3(B,E) = 0#B + E
Plots.plot(camera=(-20,20),
           xlabel = "B",
           ylabel = "E",
           zlabel = "β",
           zlim=(-20,20),
           legend = :topleft,
           title = "UB, LB of separated (red, blue), LB of composite (green)")
Plots.plot!(B,E,lb,st=:surface,
           label = "Upper bound of β",
           color=:red, alpha = 0.1)
Plots.plot!(B,E,ub,st=:surface,
           label = "Lower bound of β",
           color=:blue, alpha = 0.1)
Plots.plot!(B,E,f3,st=:surface,
           label = "Lower bound of β",
           color=:green, alpha = 0.1)
Plots.plot!(B,E,f4,st=:surface,
           label = "Upper bound of β",
           color=:orange, alpha = 0.1)
savefig("julia_merger_figure/identification_region_3d_positive_B_negative_E")

#=
using Polyhedra
using Makie
v = Polyhedra.convexhull([0, 0, 0]) +
    Polyhedra.conichull([1, 0, 0], [0, 1, 0], [0, 0, 1])
p = Polyhedra.polyhedron(v)
m = Polyhedra.Mesh(p)
Makie.mesh(m, alpha = 0.1)
Makie.wireframe!(m)
using Polyhedra
h = Polyhedra.HalfSpace([1, 2], 1) ∩
    Polyhedra.HalfSpace([-1, 0], 0) ∩
    Polyhedra.HalfSpace([0, -1], 0)
p = Polyhedra.polyhedron(h)
Plots.plot(p, ratio=:equal)
=#

println("Akkus et al 2015, mergermatches. Appendix. A.1.4. Performance in Table A1")
using CSV
#using Plots
#using JuMP, Ipopt
using Random
using Distributions
using LinearAlgebra
using Dates
using DataFrames
using Gurobi
using DataFramesMeta
using Combinatorics
using Optim
using BlackBoxOptim # on behalf of DEoptim
function expand_grid(args...)
    nargs= length(args)
    if nargs == 0
      error("expand_grid need at least one argument")
    end
    iArgs= 1:nargs
    nmc= "Var" .* string.(iArgs)
    nm= nmc
    d= map(length, args)
    orep= prod(d)
    rep_fac= [1]
    # cargs = []
    if orep == 0
        error("One or more argument(s) have a length of 0")
    end
    cargs= Array{Any}(undef,orep,nargs)
    for i in iArgs
        x= args[i]
        nx= length(x)
        orep= Int(orep/nx)
        mapped_nx= vcat(map((x,y) -> repeat([x],y), collect(1:nx), repeat(rep_fac,nx))...)
        cargs[:,i] .= x[repeat(mapped_nx,orep)]
        rep_fac= rep_fac * nx
    end
    convert(DataFrame,cargs)
end
## ------------------------------ ##
## A script to implement Fox's    ##
## matching estimator.            ##
## ------------------------------ ##

## First, set the parameters of our example
#=
set.seed(209)
library(MASS)
library(DEoptim)
library(lpSolve)
library(gregmisc)
=#
#=
matchval = function(Ab,At,Bb,Bt){
  val = Ab*At + 1.5*Bb*Bt
  return(val)
}
=#
function matchval(Ab,At,Bb,Bt)
  val = 1.0.*Ab.*At .+ 1.5.*Bb.*Bt
  return val
end

## ------------------------------ ##
## This function produces a data  ##
## set of simulated matches and   ##
## attributes.                    ##
## ------------------------------ ##
num_agents = 4
sd_err = 4.0
beta = Vector{Float64}([1.5])
#givemedata = function(num_agents, sd_err = 4){
function givemedata(num_agents::Int64, sd_err::Float64)
    #I = num_agents
    I = num_agents
    #means  = c(10,10)
    means  = [10.0, 10.0]
    #covars = matrix(c(1,0.25,0.25,1),nrow=2)
    covars = [1 0.25;0.25 1]
    #buydata = mvrnorm(I, mu=means, Sigma=covars)
    buydata = rand(Distributions.MvNormal(means, covars), I)
    #buydata = data.frame(buydata, buyid = 1:I)
    buyid = Array{Int64,1}(1:I)
    buydata = hcat(buyid, buydata')
    buydata = convert(DataFrame, buydata)
    #colnames(buydata) = c("Ab","Bb")
    names!(buydata, [:id, :Ab,  :Bb])

    #tardata = mvrnorm(I, mu=means, Sigma=covars)
    tardata = rand(Distributions.MvNormal(means, covars), I)
    #tardata = data.frame(tardata, tarid = I + 1:I)
    tarid = Array((1+I):(I+I))
    tardata = hcat(tarid, tardata')
    tardata = convert(DataFrame, tardata)
    #colnames(tardata) = c("At","Bt")
    names!(tardata, [:id, :At,  :Bt])

    #matchmaker = expand.grid(buyid = buydata$buyid, tarid = tardata$tarid)
    matchmaker = expand_grid(buyid, tarid)
    names!(matchmaker, [:buyid, :tarid])
    #matchdat = merge(matchmaker, tardata)
    matchdat = DataFrames.join(matchmaker, tardata, on = [(:tarid, :id)], kind = :left)
    #matchdat = merge(matchdat, buydata)
    matchdat = DataFrames.join(matchdat, buydata, on = [(:buyid, :id)], kind = :left)
    sort!(matchdat, (:buyid, :tarid));
    #matchdat = within(matchdat, mval <- matchval(Ab,At,Bb,Bt))
    mval = matchval(matchdat.Ab,matchdat.At,matchdat.Bb,matchdat.Bt)
    #matchdat = within(matchdat, mval <- mval + rnorm(length(matchdat$mval), mean = 0, sd_err) )
    mval = mval .+ rand(Distributions.Normal(0, sd_err), length(mval))
    matchdat = hcat(matchdat, mval)
    rename!(matchdat, :x1 => :mval)
    #matchdat = matchdat[order(matchdat$tarid),]
    #matchdat = matchdat[sort!(matchdat,:buyid),]
    #matchdat = matchdat[order(matchdat$buyid),]

    #Buy = I
    Buy = I
    #Tar = I
    Tar = I
    #eye = diag(rep(1,Tar))
    eye = convert(Array{Float64,2},LinearAlgebra.Diagonal(ones(Tar)) )
    #buy.const = matrix(rep(eye, Buy), nrow=Tar)
    global temp = eye
    for i in 1:(Tar-1)
        global temp
        temp = hcat(temp, eye)
    end
    buy_const = temp
    #tar_const = NULL
    tar_const = 0

    #for(i in 1:nrow(eye)){
    #  thisconst = rep(eye[i,], each = Tar)
    #  tar.const = rbind(tar.const, thisconst)
    #}
    tar_const = zeros(size(eye)[1],size(eye)[1]*size(eye)[2])
    for i in 1:size(eye)[1]
        for j in ((i-1)*size(eye)[2]+1):((i-1)*size(eye)[2]+size(eye)[2])
            tar_const[i,j] = 1.0
        end
    end
    #=
      obj = matchdat$mval
      all.const = rbind(tar.const, buy.const)
      rhs = matrix(rep(1, nrow(all.const)),ncol=1)
      dir = as.vector(rep(c("<="), nrow(all.const)))
      dir2 = as.vector(rep(c(">="), ncol(all.const)))
      mylp = lp("max", objective.in = obj, const.mat = all.const, const.dir =  dir, const.rhs = rhs)
    =#
    obj = matchdat.mval
    all_const = vcat(tar_const, buy_const)
    rhs = ones(size(all_const)[1])

    env = Gurobi.Env()
    setparam!(env, "Presolve", 0)
    println("Time for reading whole model:")
    @time model = gurobi_model(env;
    	name = "lp_01",
    	sense = :maximize,
    	f = obj,
    	A = all_const, # Aeq
    	b = rhs, # beq
    	lb = zeros(length(mval)),
    	ub = ones(length(mval))
    	)
    println("Time for optimizing model:")
    @time Gurobi.optimize(model)
    # show results
    objv = get_objval(model)
    println("objvalue　= ", objv)
    #matches = round(mylp$solution, 1)
    matches = get_solution(model)
    #matchdat = cbind(matchdat, matches)
    matchdat = hcat(matchdat, matches)
    rename!(matchdat, :x1 => :matches)
    #duals = lp("min", rhs, t(all.const), dir2, obj)$solution
    env = Gurobi.Env()
    setparam!(env, "Presolve", 0)
    println("Time for reading whole model:")
    @time model = gurobi_model(env;
    	name = "lp_02",
  	sense = :minimize,
    	f = rhs,
    	A = convert(Array{Float64,2}, (-all_const)'),
    	b = -obj,
        lb = zeros(length(rhs)) # dual is non-negative?
    	)
    println("Time for optimizing model:")
    @time Gurobi.optimize(model)
    objv = get_objval(model)
    println("objvalue　= ", objv)
    duals = get_solution(model) # tar price may be negative
    lo = I+1
    hi = 2*I
    #duals = data.frame(tarid=I+1:I, tarprice = duals[lo:hi])
    duals = DataFrame(tarid=Array((1+I):(I+I)), tarprice = duals[lo:hi])
    #matchdat = merge(matchdat, duals)
    matchdat = DataFrames.join(matchdat, duals, on = [(:tarid, :tarid)], kind = :left)
    #obsd = matchdat[matchdat$matches==1,]
    @linq obsd = matchdat |>
  	  where(:matches .== 1.0)
    return(obsd)
end
dat = givemedata(4, 4.0)

## ------------------------ ##
## Form the Inequalities    ##
## for both "with" and      ##
## "without" estimators.    ##
## ------------------------ ##

#=
ineq  = function(mat, idx){
  prin = mat[idx,idx]
  ineq = prin[1,1]+prin[2,2]-prin[1,2]-prin[2,1]
  return(ineq)
}
=#
function ineq(mat::Array{Float64,2}, idx::Vector{Int64})
    prin = mat[idx,idx]
    ineq = prin[1,1]+prin[2,2]-prin[1,2]-prin[2,1]
    return ineq
end

#=
with_ineq  = function(mat, vec, idx){
  iq1 = mat[idx[2], idx[2]] - mat[idx[2], idx[1]] - vec[idx[2]] + vec[idx[1]]
  iq2 = mat[idx[1], idx[1]] - mat[idx[1], idx[2]] - vec[idx[1]] + vec[idx[2]]
  ineq = c(iq1, iq2) ## ifelse( (iq1 > 0) & (iq2 > 0) , 1, 0)
  return(ineq)
}
=#
#idx = index_list[1,:]
function with_ineq(comper::Array{Float64,2}, prc::Vector{Float64}, idx::Vector{Int64})
    iq1 = comper[idx[2], idx[2]] - comper[idx[2], idx[1]] - prc[idx[2]] + prc[idx[1]]
    iq2 = comper[idx[1], idx[1]] - comper[idx[1], idx[2]] - prc[idx[1]] + prc[idx[2]]
    res_ineq = vcat(iq1, iq2) ## ifelse( (iq1 > 0) & (iq2 > 0) , 1, 0)
    return res_ineq
end

#=
score_b_with = function(beta, data){    ## Computes the maximum score function, given data and beta
  A   = data$Ab%o%data$At
  B   = data$Bb%o%data$Bt
  prc = data$tarprice
  index_list = combinations(nrow(A),2)
  ineqs = matrix(rep(-1000, 2*length(index_list[,1])), ncol = 2)

  comper = A + beta*B
  for(j in 1:length(index_list[,1])){
    ineqs[j,] <- with_ineq(comper, prc, index_list[j,])
  }
  #return(sum(ineqs[,1]>0 & ineqs[,2]>0))
   return(sum(ineqs>0))
}
=#
data = givemedata(4, 4.0)
function score_b_with(beta::Vector{Float64},data::DataFrame,num_agents::Int64)
    I = num_agents
    beta = beta[1] # for Optim
    #A   = data$Ab%o%data$At
    #B   = data$Bb%o%data$Bt
    A = kron(data.Ab, data.At') #Take care of row and column
    B = kron(data.Bb, data.Bt') #Take care of row and column
    prc = convert(Vector{Float64}, data.tarprice)
    temp = [Combinatorics.combinations(1:I,2)...]
    index_list = Array{Int64,2}(undef, length(temp), 2)
    for i in 1:length(temp)
        index_list[i,1] = temp[i][1]
        index_list[i,2] = temp[i][2]
    end
    #ineqs = matrix(rep(-1000, 2*length(index_list[,1])), ncol = 2)
    ineqs = Array{Float64,2}(undef,length(index_list[:,1]),2)
    comper = A + beta*B
    for j in 1:length(index_list[:,1])
        ineqs[j,:] = with_ineq(comper, prc, index_list[j,:])
    end
    # level of ineq is too big
    #return(sum(ineqs>0))
    #display(ineqs)
    res = sum(ineqs.>0)
    return res
end
#=
score_b = function(beta, data){    ## Computes the maximum score function, given data and beta
  A   = data$Ab%o%data$At
  B   = data$Bb%o%data$Bt
  index_list = combinations(nrow(A),2)
  ineqs = rep(-1000, length(index_list[,1]))

  comper = A + beta*B
  for(j in 1:length(index_list[,1])){
    ineqs[j] <- ineq(comper, index_list[j,])
  }
  return(sum(ineqs>0))
}
=#
data = givemedata(10, 4.0)
#data = obsdat
function score_b(beta::Vector{Float64}, data::DataFrame, num_agents::Int64)
    I = num_agents
    beta = beta[1] # for Optim
    #A   = data$Ab%o%data$At
    #B   = data$Bb%o%data$Bt
    A = kron(data.Ab, data.At') #Take care of row and column
    B = kron(data.Bb, data.Bt') #Take care of row and column
    temp = [Combinatorics.combinations(1:I,2)...]
    index_list = Array{Int64,2}(undef, length(temp), 2)
    for i in 1:length(temp)
        index_list[i,1] = temp[i][1]
        index_list[i,2] = temp[i][2]
    end
    ineqs = fill(-1000.0, length(index_list[:,1]))
    comper = A + beta*B
    for j in 1:length(index_list[:,1])
        ineqs[j] = ineq(comper, index_list[j,:])
    end
    #global ineqs
    #return(sum(ineqs>0))
    display(I)
    res = sum(ineqs.>0)
    return res
end

#score_b([10.1], data)
## ----------------------------------- ##
## Given a data set, findmin() will    ##
## find the minimum on the interval 0  ##
## to 50.                              ##
## ----------------------------------- ##
println("findmin = function(FUN){DEoptim(FUN, lower=c(0), upper = c(50)\n    ,control = DEoptim.control(storepopfreq=5, NP = 50, reltol=0.000001, steptol = 5, strategy=1, avWinner=TRUE, itermax = 200, trace = FALSE))}")

## ---------------------------- ##
## The remainder of the code    ##
## calls the previous functions ##
## to produce the estimates.    ##
## ---------------------------- ##

#=
sumfun = function(vec){
rw = nrow(vec)
co = ncol(vec)
meanmat = matrix(c(1, 1.5), ncol=co, nrow=rw, byrow=TRUE)
return(c(mean(vec-meanmat, na.rm=TRUE), sqrt(mean((vec-meanmat)^2, na.rm=TRUE))))
}
=#
function sumfun(vec)
    rw = size(vec)[1]
    co = size(vec)[2]
    #meanmat = matrix(c(1, 1.5), ncol=co, nrow=rw, byrow=TRUE)
    meanmat = [1 1.5] # true parameters
    mean = mean(vec-meanmat)
    sqrt = sqrt(mean((vec-meanmat)^2))
    return mean, sqrt
end
num_its =10
myests = Array{Float64,2}(undef, num_its*1, 1)
#mean, sqrt = sumfun(myests)

#=
myfun = function(beta, dat, i){
return(dat[i,"Ab"]*dat[,"At"]+beta*dat[i,"Bb"]*dat[,"Bt"])
}
=#
typeof([1.0])
beta = Vector{Float64}([1.5])
dat = givemedata(4, 4.0)
function myfun(beta::Vector{Float64}, dat::DataFrame, i::Int64)
    beta = beta[1] # for optim
    res = dat.Ab[i]*dat.At .+ beta.*dat.Bb[i]*dat.Bt
    return res
end
myfun(beta, dat, 1)
#=
mylogit = function(beta, dat, i){
  f = exp(myfun(beta, dat, i) - dat[, "tarprice"])
  den = sum(f)
  return(log(f[1]/den))
}
=#
function mylogit(beta::Vector{Float64}, dat::DataFrame, i::Int64)
    f = exp.(myfun(beta, dat, i) .- dat.tarprice)
    den = sum(f)
    res = log(f[1]/den) # log-likelihood
    return res
end
mylogit(beta, dat, 1)
#=
mylogiter = function(beta, dat){
logix = rep(NA, nrow(dat))
for(i in 1:nrow(dat)){
logix[i] = mylogit(beta, dat, i)
}
return(logix)
}
=#
typeof(beta[1])
function mylogiter(beta::Vector{Float64}, dat::DataFrame)
    logix = Array{Float64,1}(undef, size(dat)[1])
    for i in 1:(size(dat)[1])
        logix[i] = mylogit(beta, dat, i) # log-likelihood
    end
    llk = sum(logix) # typo for original code
    return -llk # minimizing llk by Optim
end
mylogiter(beta, dat)
obsdat = givemedata(4, 4.0)
@time @show est_res = Optim.optimize(beta -> mylogiter(beta, obsdat), zeros(length(beta)), LBFGS());
est = Optim.minimizer(est_res)
mylogiter([1.5], obsdat)
##########################
# maximum score estimator
##########################
num_agents = 4
sd_err = 4.0
beta = Vector{Float64}([1.5])
#obsdat = CSV.read("./testdata_10.csv")
num_agents = 10
#obsdat = CSV.read("./testdata_20.csv")
num_agents = 20
#delete!(obsdat,:Column1)
#maxscore_mc = function(num_agents, sd_err, num_its = 100, quietly = FALSE, withtrans = TRUE, ml = FALSE){
function maxscore_mc(num_agents::Int64, sd_err::Float64;
                     num_its::Int64 = 10,
                     quietly::Bool = false,
                     withtrans::Bool = true,
                     ml::Bool = false)
    #ptm <- proc.time()
    start = Dates.unix2datetime(time())
    #myests <<- matrix(rep(NA, num_its*1), ncol=1)
    myests = Array{Float64,2}(undef, num_its*1, 1)
    #for(i in 1:num_its){
    for i in 1:num_its
       #obsdat <<- givemedata(num_agents, sd_err)
       println("Create obsdat for iteration $i \n" )
       obsdat = givemedata(num_agents, sd_err)
       #if(ml){
       if ml == true
           #est = tryCatch(maxLik(function(beta){mylogiter(beta, obsdat)}, start = 2)$estimate, error = function(err) NA)
           @time @show est_res = Optim.optimize(beta -> mylogiter(beta, obsdat), zeros(length(beta)), LBFGS());
           est = Optim.minimizer(est_res)
           #if(!quietly){
           if quietly!== true
               #cat( "with transfer (low var) rep ", i, ": ",  as.numeric(est), " Time: ", proc.time()-ptm, "\n" )
               finish = convert(Int, Dates.value(Dates.unix2datetime(time())-start))/1000
               println("###############################################\n" )
               println("with transfer (low var) rep ", i, ": ", est, " Time: ", finish, "\n" )
               println("###############################################\n" )
               myests[i,:] = est
            end
       elseif ml!==true
           #if(withtrans){ score_bthis = function(x){ -1*score_b_with(x, obsdat) + 100000} }
           if withtrans == true
               function score_bthis_with(beta::Vector{Float64}, obsdat::DataFrame)
                   res = -1.0*score_b_with(beta, obsdat, num_agents) + 100000.0 # need to be Float64 for bboptimize
                   println("score: $(res) \n" )
                   return res
               end
               m_res = BlackBoxOptim.bboptimize(beta -> score_bthis_with(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta), Method = :de_rand_1_bin, MaxSteps=200)
               println("score: ",score_bthis_with(m_res.archive_output.best_candidate, obsdat))
           #else { score_bthis = function(x){ -1*score_b(x, obsdat) + 100000 }}
           else#if withtrans == false
               function score_bthis(beta::Vector{Float64}, obsdat::DataFrame)
                   res = -1.0*score_b(beta, obsdat, num_agents) + 100000.0 # need to be Float64 for bboptimize
                   println("score:$(res) \n" )
                   return res
               end
               println("******** test **********")
               m_res = BlackBoxOptim.bboptimize(beta -> score_bthis(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta), Method = :de_rand_1_bin,  MaxSteps=100)
               println("score: ", score_bthis(m_res.archive_output.best_candidate, obsdat))
               println("score of correct: ", score_b(m_res.archive_output.best_candidate, obsdat, num_agents))
               println("TRUE score of correct: ", score_b([1.5], obsdat, num_agents))
               #score_b([2.3], obsdat, 20)
           end
           #m = findmin(score_bthis)
           #@time @show m_res = Optim.optimize(beta -> score_bthis(beta, obsdat), zeros(length(beta)), LBFGS())
           #m = Optim.minimizer(m_res)
           #m_res = BlackBoxOptim.bboptimize(beta -> score_bthis(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta))
           #m_res = BlackBoxOptim.bboptimize(beta -> score_bthis(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta), Method = :de_rand_1_bin)
           m = m_res.archive_output.best_candidate
           #if(!quietly){
           if !quietly == true
              #cat( "with transfer (low var) rep ", i, ": ",  as.numeric(m$optim$bestmem), " Time: ", proc.time()-ptm, "\n" )
              finish = convert(Int, Dates.value(Dates.unix2datetime(time())-start))/1000
              println("###############################################\n" )
              println("with transfer (low var) rep ", i, ": ", m, " \nTime: ", finish, "\n" )
              println("###############################################\n" )
              myests[i,:] = m
           end
           #myests[i,] = as.numeric(m$optim$bestmem)
           myests[i,:] = m
       end
   end
   #runtime = proc.time() - ptm
   finish = convert(Int, Dates.value(Dates.unix2datetime(time())-start))/1000
   #mean = mean(vec-meanmat)
   #sqrt = sqrt(mean((vec-meanmat)^2))
   meanmat = 1.5 # true parameter
   res_mean = mean(myests.-meanmat)
   res_sqrt = sqrt(mean((myests.-meanmat).^2))
   #res = c(sumfun(myests), runtime[1])
   #res = sumfun(myests)
   #myests <<- myests
   global myests
   @show myests
   return res_mean, res_sqrt
end

@time @show res = maxscore_mc(10, 1.0, ml =true, withtrans = true, quietly = false)
@time @show res = maxscore_mc(10, 1.0, ml =false, withtrans = true, quietly = false)
@time @show res = maxscore_mc(20, 1.0, ml =false, withtrans = false, quietly = false)
@time @show res = maxscore_mc(100, 1.0, ml =false, withtrans = false, quietly = false)
@time @show res = maxscore_mc(100, 1.0, ml =false, withtrans = true, quietly = false)
#=
set.seed(209)
x10 =maxscore_mc(num_agents = 100, sd_err = 1, ml = TRUE)
x11 =maxscore_mc(num_agents = 100, sd_err = 5, ml = TRUE)
x12 =maxscore_mc(num_agents = 100, sd_err = 20, ml = TRUE)
=#
@time @show x10 = maxscore_mc(100, 1.0, ml =true, withtrans = true, quietly = false)
@time @show x11 = maxscore_mc(100, 5.0, ml =true, withtrans = true, quietly = false)
@time @show x12 = maxscore_mc(100, 20.0, ml =true, withtrans = true, quietly = false)

#=
set.seed(209)
x1 = maxscore_mc(num_agents = 100, sd_err = 1, withtrans=FALSE)
x2 = maxscore_mc(num_agents = 100, sd_err = 1)
x3 = maxscore_mc(num_agents = 100, sd_err = 5, withtrans=FALSE)
x4 = maxscore_mc(num_agents = 100, sd_err = 5)
x5 = maxscore_mc(num_agents = 100, sd_err = 20, withtrans = FALSE)
x6 = maxscore_mc(num_agents = 100, sd_err = 20)
tab1 = rbind(x1,x2,x3,x4, x5, x6)
colnames(tab1) = c("bias", "RMSE", "Run Time")
row.names(tab1) = c("without.low.var", "with.low.var", "without.med.var","with.med.var", "without.high.var", "with.high.var")
round(tab1,3)
=#
@time @show x1 = maxscore_mc(100, 1.0, withtrans = false, quietly = false)
@time @show x2 = maxscore_mc(100, 1.0, withtrans = true, quietly = false)
@time @show x3 = maxscore_mc(100, 5.0, withtrans = false, quietly = false)
@time @show x4 = maxscore_mc(100, 5.0, withtrans = true, quietly = false)
@time @show x3 = maxscore_mc(100, 20.0, withtrans = false, quietly = false)
@time @show x4 = maxscore_mc(100, 20.0, withtrans = true, quietly = false)

#=
x7 = maxscore_mc(num_agents = 100, sd_err = 1)
x8 = maxscore_mc(num_agents = 100, sd_err = 5)
x9 = maxscore_mc(num_agents = 100, sd_err = 20)
rbind(x7,x8,x9)
=#

println("## Change the functions to \"simultaneous functions(Multinomial Logit)\" for this")







println("Akkus et al 2015, mergermatchesmultiple. Appendix. A.1.5. Performance in Table A2")
obsdat = CSV.read("./testdata_10_non.csv")
num_agents = 10
obsdat = CSV.read("./testdata_20.csv")
num_agents = 20
delete!(obsdat,:Column1)
function matchval2(Ab,At,Bb,Bt,Ct)
  val = 1.0.*Ab.*At .+ 1.5.*Bb.*Bt .+ 2.0 .* Ct
  return val
end
function givemedata2(num_agents::Int64, sd_err::Float64)
    #I = num_agents
    I = num_agents
    #means  = c(10,10)
    means  = [10.0, 10.0]
    #covars = matrix(c(1,0.25,0.25,1),nrow=2)
    covars = [1 0.25;0.25 1]
    #buydata = mvrnorm(I, mu=means, Sigma=covars)
    buydata = rand(Distributions.MvNormal(means, covars), I)
    #buydata = data.frame(buydata, buyid = 1:I)
    buyid = Array{Int64,1}(1:I)
    buydata = hcat(buyid, buydata')
    buydata = convert(DataFrame, buydata)
    #colnames(buydata) = c("Ab","Bb")
    names!(buydata, [:id, :Ab,  :Bb])

    #tardata = mvrnorm(I, mu=means, Sigma=covars)
    tardata = rand(Distributions.MvNormal(means, covars), I)
    #tardata = data.frame(tardata, tarid = I + 1:I)
    tarid = Array((1+I):(I+I))
    #--------------------#
    # non-interactive term
    #--------------------#
    println("non-interactive Match specific term: Ct = rnorm(I, 10, 1)")
    Ct = rand(Distributions.Normal(10, 1), I)
    tardata = hcat(tarid, tardata', Ct)
    tardata = convert(DataFrame, tardata)
    #colnames(tardata) = c("At","Bt")
    names!(tardata, [:id, :At, :Bt, :Ct])

    #matchmaker = expand.grid(buyid = buydata$buyid, tarid = tardata$tarid)
    matchmaker = expand_grid(buyid, tarid)
    names!(matchmaker, [:buyid, :tarid])
    #matchdat = merge(matchmaker, tardata)
    matchdat = DataFrames.join(matchmaker, tardata, on = [(:tarid, :id)], kind = :left)
    #matchdat = merge(matchdat, buydata)
    matchdat = DataFrames.join(matchdat, buydata, on = [(:buyid, :id)], kind = :left)
    sort!(matchdat, (:buyid, :tarid));
    #matchdat = within(matchdat, mval <- matchval(Ab,At,Bb,Bt))
    mval = matchval2(matchdat.Ab,matchdat.At,matchdat.Bb,matchdat.Bt,matchdat.Ct)
    #matchdat = within(matchdat, mval <- mval + rnorm(length(matchdat$mval), mean = 0, sd_err) )
    mval = mval .+ rand(Distributions.Normal(0, sd_err), length(mval))
    matchdat = hcat(matchdat, mval)
    rename!(matchdat, :x1 => :mval)
    #matchdat = matchdat[order(matchdat$tarid),]
    #matchdat = matchdat[sort!(matchdat,:buyid),]
    #matchdat = matchdat[order(matchdat$buyid),]

    #Buy = I
    Buy = I
    #Tar = I
    Tar = I
    #eye = diag(rep(1,Tar))
    eye = convert(Array{Float64,2},LinearAlgebra.Diagonal(ones(Tar)) )
    #buy.const = matrix(rep(eye, Buy), nrow=Tar)
    global temp = eye
    for i in 1:(Tar-1)
        global temp
        temp = hcat(temp, eye)
    end
    buy_const = temp
    #tar_const = NULL
    tar_const = 0

    #for(i in 1:nrow(eye)){
    #  thisconst = rep(eye[i,], each = Tar)
    #  tar.const = rbind(tar.const, thisconst)
    #}
    tar_const = zeros(size(eye)[1],size(eye)[1]*size(eye)[2])
    for i in 1:size(eye)[1]
        for j in ((i-1)*size(eye)[2]+1):((i-1)*size(eye)[2]+size(eye)[2])
            tar_const[i,j] = 1.0
        end
    end
    #=
      obj = matchdat$mval
      all.const = rbind(tar.const, buy.const)
      rhs = matrix(rep(1, nrow(all.const)),ncol=1)
      dir = as.vector(rep(c("<="), nrow(all.const)))
      dir2 = as.vector(rep(c(">="), ncol(all.const)))
      mylp = lp("max", objective.in = obj, const.mat = all.const, const.dir =  dir, const.rhs = rhs)
    =#
    obj = matchdat.mval
    all_const = vcat(tar_const, buy_const)
    rhs = ones(size(all_const)[1])

    env = Gurobi.Env()
    setparam!(env, "Presolve", 0)
    println("Time for reading whole model:")
    @time model = gurobi_model(env;
    	name = "lp_01",
    	sense = :maximize,
    	f = obj,
    	A = all_const, # Aeq
    	b = rhs, # beq
    	lb = zeros(length(mval)),
    	ub = ones(length(mval))
    	)
    println("Time for optimizing model:")
    @time Gurobi.optimize(model)
    # show results
    objv = get_objval(model)
    println("objvalue　= ", objv)
    #matches = round(mylp$solution, 1)
    matches = get_solution(model)
    #matchdat = cbind(matchdat, matches)
    matchdat = hcat(matchdat, matches)
    rename!(matchdat, :x1 => :matches)
    #duals = lp("min", rhs, t(all.const), dir2, obj)$solution
    env = Gurobi.Env()
    setparam!(env, "Presolve", 0)
    println("Time for reading whole model:")
    @time model = gurobi_model(env;
    	name = "lp_02",
  	sense = :minimize,
    	f = rhs,
    	A = convert(Array{Float64,2}, (-all_const)'),
    	b = -obj,
        lb = zeros(length(rhs)) # dual is non-negative?
    	)
    println("Time for optimizing model:")
    @time Gurobi.optimize(model)
    objv = get_objval(model)
    println("objvalue　= ", objv)
    duals = get_solution(model) # tar price may be negative
    lo = I+1
    hi = 2*I
    #duals = data.frame(tarid=I+1:I, tarprice = duals[lo:hi])
    duals = DataFrame(tarid=Array((1+I):(I+I)), tarprice = duals[lo:hi])
    #matchdat = merge(matchdat, duals)
    matchdat = DataFrames.join(matchdat, duals, on = [(:tarid, :tarid)], kind = :left)
    #obsd = matchdat[matchdat$matches==1,]
    @linq obsd = matchdat |>
  	  where(:matches .== 1.0)
    return(obsd)
end
dat = givemedata2(4, 4.0)
#=
sdz = rep(NA, 100)
for(i in 1:100){
  dat = givemedata2(100)
  sdz[i] = sd(dat$mval)
}
mean(sdz)
=#
data = givemedata2(4, 4.0)
#beta = [1.0 1.5 2.0]
beta = Vector{Float64}([1.5, 2.0, 1.0])
println("score_b_with_non ")
function score_b_with_non(beta::Vector{Float64},data::DataFrame, num_agents::Int64)
    I = num_agents
    #beta = beta[1] # for Optim
    #A   = data$Ab%o%data$At
    #B   = data$Bb%o%data$Bt
    A = kron(data.Ab, data.At') #Take care of row and column
    B = kron(data.Bb, data.Bt') #Take care of row and column
    # C   = rep(1,length(data$Ct))%o%data$Ct
    C = kron(ones(length(data.Ct)),data.Ct')
    prc = convert(Vector{Float64}, data.tarprice)
    temp = [Combinatorics.combinations(1:I,2)...]
    index_list = Array{Int64,2}(undef, length(temp), 2)
    for i in 1:length(temp)
        index_list[i,1] = temp[i][1]
        index_list[i,2] = temp[i][2]
    end
    #ineqs = matrix(rep(-1000, 2*length(index_list[,1])), ncol = 2)
    ineqs = Array{Float64,2}(undef,length(index_list[:,1]),2)
    # comper = beta[3]*A + beta[1]*B + beta[2]*C
    comper = beta[3]*A + beta[1]*B .+ beta[2]*C
    for j in 1:length(index_list[:,1])
        ineqs[j,:] = with_ineq(comper, prc, index_list[j,:])
    end
    # level of ineq is too big
    #return(sum(ineqs>0))
    #display(ineqs)
    res = sum(ineqs.>0)
    return res
end
#=
score_b = function(beta, data){    ## Computes the maximum score function, given data and beta
  A   = data$Ab%o%data$At
  B   = data$Bb%o%data$Bt
  index_list = combinations(nrow(A),2)
  ineqs = rep(-1000, length(index_list[,1]))

  comper = A + beta*B
  for(j in 1:length(index_list[,1])){
    ineqs[j] <- ineq(comper, index_list[j,])
  }
  return(sum(ineqs>0))
}
=#
data = givemedata(10, 4.0)
#data = obsdat

println("score_b_non")
beta = Vector{Float64}([1.5, 2.0, 1.0])
beta = Vector{Float64}([2.0, 1.0]) # constant is not identified

function score_b_non(beta::Vector{Float64}, data::DataFrame, num_agents::Int64)
    I = num_agents
    #beta = beta[1] # for Optim
    #A   = data$Ab%o%data$At
    #B   = data$Bb%o%data$Bt
    A = kron(data.Ab, data.At') #Take care of row and column
    B = kron(data.Bb, data.Bt') #Take care of row and column
    # C   = rep(1,length(data$Ct))%o%data$Ct
    C = kron(ones(length(data.Ct)),data.Ct')
    temp = [Combinatorics.combinations(1:I,2)...]
    index_list = Array{Int64,2}(undef, length(temp), 2)
    for i in 1:length(temp)
        index_list[i,1] = temp[i][1]
        index_list[i,2] = temp[i][2]
    end
    ineqs = fill(-1000.0, length(index_list[:,1]))
    # comper = A + beta[1]*B + beta[2]*C # typo?!
    println("Need to fix constant coefficient = 1")
    comper = A + beta[1]*B .+ beta[2]*C
    for j in 1:length(index_list[:,1])
        ineqs[j] = ineq(comper, index_list[j,:])
    end
    #global ineqs
    #return(sum(ineqs>0))
    display(I)
    res = sum(ineqs.>0)
    return res
end
println("maxscore_mc2")
beta = Vector{Float64}([1.5, 2.0, 1.0])
beta = Vector{Float64}([2.0, 1.0]) # constant is not identified
obsdat = CSV.read("./testdata_100_non.csv")
num_agents = 100
#obsdat = CSV.read("./testdata_20.csv")
#num_agents = 20
delete!(obsdat,:Column1)
obsdat = givemedata2(100, sd_err)
function maxscore_mc2(num_agents::Int64, sd_err::Float64;
                     num_its::Int64 = 10,
                     quietly::Bool = false,
                     withtrans::Bool = true,
                     misspec::Bool = false)
    #ptm <- proc.time()
    start = Dates.unix2datetime(time())
    #myests <<- matrix(rep(NA, num_its*1), ncol=1)
    myests = zeros(num_its, 3)
    #for(i in 1:num_its){
    for i in 1:num_its
       #obsdat <<- givemedata(num_agents, sd_err)
       println("Create obsdat for iteration $i \n" )
       println("Use non-interactive term \n" )
       obsdat = givemedata2(num_agents, sd_err)
       #if(ml){
       if misspec == true
           #if(withtrans){ score_bthis = function(x){ -1*score_b_with(x, obsdat) + 100000} }
           if withtrans == true
               function score_bthis_with(beta::Vector{Float64}, obsdat::DataFrame)
                   res = -1.0*score_b_with(beta, obsdat, num_agents) + 100000.0 # need to be Float64 for bboptimize
                   println("score: $(res) \n" )
                   return res
               end
               m_res = BlackBoxOptim.bboptimize(beta -> score_bthis_with(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta), Method = :de_rand_1_bin, MaxSteps=200)
               println("score: ",score_bthis_with(m_res.archive_output.best_candidate, obsdat))
           #else { score_bthis = function(x){ -1*score_b(x, obsdat) + 100000 }}
           else#if withtrans == false
               function score_bthis(beta::Vector{Float64}, obsdat::DataFrame)
                   res = -1.0*score_b(beta, obsdat, num_agents) + 100000.0 # need to be Float64 for bboptimize
                   println("score:$(res) \n" )
                   return res
               end
               println("******** test **********")
               m_res = BlackBoxOptim.bboptimize(beta -> score_bthis(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta), Method = :de_rand_1_bin,  MaxSteps=100)
               println("score: ", score_bthis(m_res.archive_output.best_candidate, obsdat))
               println("score of correct: ", score_b(m_res.archive_output.best_candidate, obsdat, num_agents))
               println("TRUE score of correct: ", score_b([1.5], obsdat, num_agents))
               #score_b([2.3], obsdat, 20)
           end
           #m = findmin(score_bthis)
           #@time @show m_res = Optim.optimize(beta -> score_bthis(beta, obsdat), zeros(length(beta)), LBFGS())
           #m = Optim.minimizer(m_res)
           #m_res = BlackBoxOptim.bboptimize(beta -> score_bthis(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta))
           #m_res = BlackBoxOptim.bboptimize(beta -> score_bthis(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta), Method = :de_rand_1_bin)
           m = m_res.archive_output.best_candidate
           #if(!quietly){
           if !quietly == true
              #cat( "with transfer (low var) rep ", i, ": ",  as.numeric(m$optim$bestmem), " Time: ", proc.time()-ptm, "\n" )
              finish = convert(Int, Dates.value(Dates.unix2datetime(time())-start))/1000
              println("###############################################\n" )
              println("with transfer (low var) rep ", i, ": ", m, " \nTime: ", finish, "\n" )
              println("###############################################\n" )
              myests[i,:] = m
           end
           #myests[i,] = as.numeric(m$optim$bestmem)
           myests[i,:] = m
       elseif misspec!==true
           #if(withtrans){ score_bthis = function(x){ -1*score_b_with(x, obsdat) + 100000} }
           if withtrans == true
               function score_bthis_non_with(beta::Vector{Float64}, obsdat::DataFrame)
                   res = -1.0*score_b_with_non(beta, obsdat, num_agents) + 100000.0 # need to be Float64 for bboptimize
                   println("score: $(res) \n" )
                   return res
               end
               #findmin_3d(score_bthis)
               m_res = BlackBoxOptim.bboptimize(beta -> score_bthis_non_with(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta), Method = :de_rand_1_bin, MaxSteps=2000)
               println("score: ",score_bthis_non_with(m_res.archive_output.best_candidate, obsdat))
               println("TRUE score_b_non of correct: ", score_b_with_non([1.5,2.0,1.0], obsdat, num_agents))
               #score_b([2.3], obsdat, 20)
           #else { score_bthis = function(x){ -1*score_b(x, obsdat) + 100000 }}
           else#if withtrans == false
               function score_bthis_non(beta::Vector{Float64}, obsdat::DataFrame)
                   res = -1.0*score_b_non(beta, obsdat, num_agents) + 100000.0 # need to be Float64 for bboptimize
                   println("score:$(res) \n" )
                   return res
               end
               println("******** test **********")
               m_res = BlackBoxOptim.bboptimize(beta -> score_bthis_non(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta)-1, Method = :de_rand_1_bin,  MaxSteps=2000)
               println("score: ", score_bthis_non(m_res.archive_output.best_candidate, obsdat))
               println("score_b of correct: ", score_b(m_res.archive_output.best_candidate, obsdat, num_agents))
               println("TRUE score_b_non of correct: ", score_b_non([1.5,2.0], obsdat, num_agents))
               #score_b([2.3], obsdat, 20)
           end
           #m = findmin(score_bthis)
           #@time @show m_res = Optim.optimize(beta -> score_bthis(beta, obsdat), zeros(length(beta)), LBFGS())
           #m = Optim.minimizer(m_res)
           #m_res = BlackBoxOptim.bboptimize(beta -> score_bthis(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta))
           #m_res = BlackBoxOptim.bboptimize(beta -> score_bthis(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta), Method = :de_rand_1_bin)
           m = m_res.archive_output.best_candidate
           #if(!quietly){
           if !quietly == true
              #cat( "with transfer (low var) rep ", i, ": ",  as.numeric(m$optim$bestmem), " Time: ", proc.time()-ptm, "\n" )
              finish = convert(Int, Dates.value(Dates.unix2datetime(time())-start))/1000
              println("###############################################\n" )
              println("with transfer (low var) rep ", i, ": ", m, " \nTime: ", finish, "\n" )
              println("###############################################\n" )
              myests[i,1] = m[1]
              myests[i,2] = m[2]
           end
           #myests[i,] = as.numeric(m$optim$bestmem)
           myests[i,1] = m[1]
           myests[i,2] = m[2]
       end
   end
   #runtime = proc.time() - ptm
   finish = convert(Int, Dates.value(Dates.unix2datetime(time())-start))/1000
   #mean = mean(vec-meanmat)
   #sqrt = sqrt(mean((vec-meanmat)^2))
   meanmat = [1.5 2.0 1.0]
   println("true param: ", meanmat)
   res_mean = mean(myests.-meanmat)
   res_sqrt = sqrt(mean((myests.-meanmat).^2))
   #res = c(sumfun(myests), runtime[1])
   #res = sumfun(myests)
   #myests <<- myests
   global myests
   @show myests
   return res_mean, res_sqrt
end

@time @show res_mean, res_sqrt = maxscore_mc2(100, 1.0, misspec =false, withtrans = true, quietly = false)
@time @show res_mean, res_sqrt = maxscore_mc2(100, 1.0, misspec =false, withtrans = false, quietly = false)
myests
@time @show res_mean, res_sqrt = maxscore_mc2(100, 1.0, misspec =true, withtrans = true, quietly = false)
@time @show res_mean, res_sqrt = maxscore_mc2(100, 1.0, misspec =true, withtrans = false, quietly = false)













println("Akkus et al 2015, mergermatchesmultiple. Appendix. A.2.3 Performance in Table A5")

mdat = dat
#mvaluator = function(vec, mdat){
function mvaluator(vec,mdat)
    #  mdate  = mdat$mval[mdat$buyid == sort(vec)[1] & mdat$tarid == sort(vec)[2]]
    mdate  = mdat.mval[mdat.buyid == sort(vec)[1] & mdat.tarid == sort(vec)[2]]
    #  return(mdate)
    return mdate
end

#datagetter = function(num_agents, sd_err, num_its = 100){
function datagetter_two_sided(num_agents::Int64, sd_err, num_its)
  #sdz = rep(NA, num_its)
  sdz = Array{Float64,1}(undef,num_its)
  #for(i in 1:num_its){
  for i in 1:num_its
    obsdat = "not yet"
    ct = 1
    while obsdat == "not yet"
      println("Generating Data.  Try ", ct, "\n")
      obsdat = givemedata(num_agents, sd_err)
      ct = ct +1
    end
    println("\n Data Generation Success!\n")
    sdz[i] = std(obsdat.mval)
  end
  return sdz
end
#sdzo2 = datagetter(100, 20, num_its=10)
#sdzo22 = datagetter(100, 5, num_its=100)
sdzo2 = datagetter_two_sided(100, 20.0, 10)
sdzo22 = datagetter_two_sided(100, 5.0, 10)
## ------------------------------ ##
## This function produces a data  ##
## set of simulated matches and   ##
## attributes.                    ##
## ------------------------------ ##
num_agents = 10
function givemedata_one_sided(num_agents::Int64, sd_err::Float64)
    #I = num_agents
    I = num_agents
    #means  = c(10,10)
    means  = [10.0, 10.0]
    #covars = matrix(c(1,0.25,0.25,1),nrow=2)
    covars = [1 0.25;0.25 1]
    #buydata = mvrnorm(I, mu=means, Sigma=covars)
    #buydata = rand(Distributions.MvNormal(means, covars), I)
    buydata = rand(Distributions.MvNormal(means, covars), 2*I)
    #buydata = data.frame(buydata, buyid = 1:I)
    #buyid = Array{Int64,1}(1:I)
    buyid = Array{Int64,1}(1:(2*I))
    buydata = hcat(buyid, buydata')
    buydata = convert(DataFrame, buydata)
    #colnames(buydata) = c("Ab","Bb")
    names!(buydata, [:id, :Ab,  :Bb])

    #tardata = mvrnorm(I, mu=means, Sigma=covars)
    #tardata = rand(Distributions.MvNormal(means, covars), I)
    tardata = rand(Distributions.MvNormal(means, covars), 2*I)
    #tardata = data.frame(tardata, tarid = I + 1:I)
    #tarid = Array((1+I):(I+I))
    tarid = Array{Int64,1}(1:(2*I))
    tardata = hcat(tarid, tardata')
    tardata = convert(DataFrame, tardata)
    #colnames(tardata) = c("At","Bt")
    names!(tardata, [:id, :At,  :Bt])

    #matchmaker = data.frame(combinations(2*I, 2))
    temp = collect(combinations(1:2*I,2))
    matchmaker = Array{Int64,2}(undef,size(temp)[1],2)
    for i in 1: size(temp)[1]
        matchmaker[i,1] = temp[i][1]
        matchmaker[i,2] = temp[i][2]
    end
    matchmaker = convert(DataFrame, matchmaker)
    names!(matchmaker, [:buyid, :tarid])
    #matchdat = merge(matchmaker, tardata)
    matchdat = DataFrames.join(matchmaker, tardata, on = [(:tarid, :id)], kind = :left)
    #matchdat = merge(matchdat, buydata)
    matchdat = DataFrames.join(matchdat, buydata, on = [(:buyid, :id)], kind = :left)
    sort!(matchdat, (:buyid, :tarid));
    #matchdat = within(matchdat, mval <- matchval(Ab,At,Bb,Bt))
    mval = matchval(matchdat.Ab,matchdat.At,matchdat.Bb,matchdat.Bt)
    #matchdat = within(matchdat, mval <- mval + rnorm(length(matchdat$mval), mean = 0, sd_err) )
    mval = mval .+ rand(Distributions.Normal(0, sd_err), length(mval))
    matchdat = hcat(matchdat, mval)
    #matchdat = matchdat[order(matchdat$tarid),]
    #matchdat = matchdat[sort!(matchdat,:buyid),]
    #matchdat = matchdat[order(matchdat$buyid),]
    rename!(matchdat, :x1 => :mval)

    #buy.const =  t(model.matrix(?0+as.factor(matchdat$tarid)))
    buy_const = zeros(2*I-1, size(matchdat)[1])
    for i in 1:2*I-1
        diff = 2*I-1
        head_id = i
        for j in 1:i
            if j == 1
                temp_id = head_id
            else
                temp_id = head_id + diff
            end
            #println("id", i, "  ",temp_id)
            buy_const[i, temp_id] = 1.0
            head_id = temp_id
            diff = diff - 1
        end
    end
    buy_const[10,136]
    buy_const[9,44]
    #tar.const = t(model.matrix(?0+as.factor(matchdat$buyid)))
    tar_const = zeros(2*I-1, size(matchdat)[1])
    head_id = 1
    diff = 2*I-1-1
    for i = 1:2*I-1
        #global diff, head_id
        tar_const[i,head_id:(head_id+diff)] = ones(diff+1)
        head_id = (head_id+diff+1)
        diff = diff - 1
    end
    tar_const
    #obuy.const = buy.const[-(nrow(buy.const)),]
    obuy_const = buy_const[1:(size(buy_const)[1])-1,:]
    #otar.const = tar.const[-1,]
    otar_const = tar_const[2:(size(buy_const)[1]),:]
    #o.const = obuy.const + otar.const
    o_const = obuy_const .+ otar_const

    #=
    #Buy = I
    Buy = I
    #Tar = I
    Tar = I
    #eye = diag(rep(1,Tar))
    eye = convert(Array{Float64,2},LinearAlgebra.Diagonal(ones(Tar)) )
    #buy.const = matrix(rep(eye, Buy), nrow=Tar)
    global temp = eye
    for i in 1:(Tar-1)
        global temp
        temp = hcat(temp, eye)
    end
    buy_const = temp
    #tar_const = NULL
    tar_const = 0

    #for(i in 1:nrow(eye)){
    #  thisconst = rep(eye[i,], each = Tar)
    #  tar.const = rbind(tar.const, thisconst)
    #}
    tar_const = zeros(size(eye)[1],size(eye)[1]*size(eye)[2])
    for i in 1:size(eye)[1]
        for j in ((i-1)*size(eye)[2]+1):((i-1)*size(eye)[2]+size(eye)[2])
            tar_const[i,j] = 1.0
        end
    end
    =#

    #=
    obj = matchdat$mval
    all.const = rbind(buy.const, tar.const, o.const)
    rhs = matrix(rep(1, nrow(all.const)),ncol=1)
    dir = as.vector(rep(c("<="), nrow(all.const)))

    mylp = lp("max", objective.in = obj, const.mat = all.const, const.dir =  dir, const.rhs = rhs)
    =#
    obj = matchdat.mval
    all_const = vcat(buy_const, tar_const, o_const)
    rhs = ones(size(all_const)[1])
    all_const[5,71]
    env = Gurobi.Env()
    setparam!(env, "Presolve", 0)
    println("Time for reading whole model:")
    @time model = gurobi_model(env;
    	name = "lp_01",
    	sense = :maximize,
    	f = obj,
    	A = all_const, # Aeq
    	b = rhs, # beq
    	lb = zeros(length(mval)),
    	ub = ones(length(mval))
    	)
    println("Time for optimizing model:")
    @time Gurobi.optimize(model)
    # show results
    objv = get_objval(model)
    println("objvalue　= ", objv)
    #matches = round(mylp$solution, 1)
    matches = get_solution(model)
    #matchdat = cbind(matchdat, matches)
    matchdat = hcat(matchdat, matches)
    rename!(matchdat, :x1 => :matches)
    matchdat[matchdat.matches .>=0.001 ,:]
    obsd = matchdat[matchdat.matches .== 1,:]
    dat = "not yet"
    #=
    row.names(all.const) = gsub("as.factor\\(matchdat\\$buyid\\)", "", gsub("as.factor\\(matchdat\\$tarid\\)", "", row.names(all.const)))
    dir2 = as.vector(rep(c(">="), ncol(all.const)))
    rhs = matrix(rep(1, nrow(all.const)),ncol=1)

    duals = lp("min", rhs, t(all.const), dir2, obj)$solution
    prices = t(model.matrix(~0+as.factor(row.names(all.const))))%*%duals
    row.names(prices) = gsub("as.factor\\(row.names\\(all.const\\)\\)", "", row.names(prices))
    pricesdf = data.frame(tarprice = prices, tarid = row.names(prices))

    dat = merge(obsd, pricesdf)
    }=#
    if size(obsd)[1]==I
        #duals = lp("min", rhs, t(all.const), dir2, obj)$solution
        env = Gurobi.Env()
        setparam!(env, "Presolve", 0)
        println("Time for reading whole model:")
        @time model = gurobi_model(env;
            name = "lp_02",
        sense = :minimize,
            f = rhs,
            A = convert(Array{Float64,2}, (-all_const)'),
            b = -obj,
            lb = zeros(length(rhs)) # dual is non-negative?
            )
        println("Time for optimizing model:")
        @time Gurobi.optimize(model)
        objv = get_objval(model)
        println("objvalue　= ", objv)
        #duals = lp("min", rhs, t(all.const), dir2, obj)$solution
        duals = get_solution(model) # tar price may be negative
        #prices = t(model.matrix(~0+as.factor(row.names(all.const))))%*%duals
        lo = I+1
        hi = 2*I
        #duals = data.frame(tarid=I+1:I, tarprice = duals[lo:hi])
        duals = DataFrame(tarid=Array((1+I):(I+I)), tarprice = duals[lo:hi])
        #matchdat = merge(matchdat, duals)
        matchdat = DataFrames.join(matchdat, duals, on = [(:tarid, :tarid)], kind = :left)
        #obsd = matchdat[matchdat$matches==1,]
        @linq obsd = matchdat |>
          where(:matches .== 1.0)
        dat = obsd
    end
    return(dat)
end



function datagetter_one_sided(num_agents::Int64, sd_err, num_its)
  #sdz = rep(NA, num_its)
  sdz = Array{Float64,1}(undef,num_its)
  #for(i in 1:num_its){
  for i in 1:num_its
    obsdat = "not yet"
    ct = 1
    while obsdat == "not yet"
      println("Generating Data.  Try ", ct, "\n")
      obsdat = givemedata_one_sided(num_agents, sd_err)
      ct = ct +1
    end
    println("\n Data Generation Success!\n")
    sdz[i] = std(obsdat.mval)
  end
  return sdz
end

sdzo = datagetter_one_sided(100, 20.0, num_its)
sdzo1 = datagetter_one_sided(100, 5.0, num_its)

#c(mean(sdzo1), mean(sdzo))   ## 1-sided matching
[mean(sdzo1), mean(sdzo)]
#c(mean(sdzo22) , mean( sdzo2))  ## 2-sided matching
[mean(sdzo22), mean(sdzo2)]







#=
with_ineq_onesided  = function(matBT, matBB, matTT, vec, idx){  ## matBT is the outerproduct of B and T, matBB is the outer product of B and B, matTT is the outer product of T and T
  iq1 = matBT[idx[2], idx[2]] - matBT[idx[2], idx[1]] - vec[idx[2]] + vec[idx[1]]
  iq2 = matBT[idx[1], idx[1]] - matBT[idx[1], idx[2]] - vec[idx[1]] + vec[idx[2]]
  iq3 = matBT[idx[1], idx[1]] + matBT[idx[2], idx[2]] - vec[idx[1]] - vec[idx[2]] - matBB[idx[1],idx[2]]
  iq4 = vec[idx[1]] + vec[idx[2]] - matTT[idx[1], idx[2]]
  ineq = c(iq1, iq2, iq3, iq4)
  return(ineq)
}
=#
function with_ineq_one_sided(matBT, matBB, matTT, vec, idx)
    iq1 = matBT[idx[2], idx[2]] - matBT[idx[2], idx[1]] - vec[idx[2]] + vec[idx[1]]
    iq2 = matBT[idx[1], idx[1]] - matBT[idx[1], idx[2]] - vec[idx[1]] + vec[idx[2]]
    iq3 = matBT[idx[1], idx[1]] + matBT[idx[2], idx[2]] - vec[idx[1]] - vec[idx[2]] - matBB[idx[1],idx[2]]
    iq4 = vec[idx[1]] + vec[idx[2]] - matTT[idx[1], idx[2]]
    ineq = vcat(iq1, iq2, iq3, iq4)
    return ineq
end

#=
score_b_with_one = function(beta, data){    ## Computes the maximum score function, given data and beta
  A   = data$Ab%o%data$At
  B   = data$Bb%o%data$Bt

  Ab   = data$Ab%o%data$Ab
  Bb   = data$Bb%o%data$Bb

  At   = data$At%o%data$At
  Bt   = data$Bt%o%data$Bt

  prc = data$tarprice
  index_list = combinations(nrow(A),2)
  ineqs = matrix(rep(-1000, 4*length(index_list[,1])), ncol = 4)

  comper = A + beta*B
  comperBB = Ab +beta*Bb
  comperTT = At +beta*Bt
  for(j in 1:length(index_list[,1])){
    ineqs[j,] <- with_ineq_onesided(comper, comperBB, comperTT, prc, index_list[j,])
  }
  #return(sum(ineqs[,1]>0 & ineqs[,2]>0))
  return(sum(ineqs>0))
}
=#
beta = Vector{Float64}([2.0, 1.0])
function score_b_with_one(beta::Vector{Float64},data::DataFrame)
    #I = num_agents
    beta = beta[1] # for Optim
    #A   = data$Ab%o%data$At
    #B   = data$Bb%o%data$Bt
    A = kron(data.Ab, data.At') #Take care of row and column
    B = kron(data.Bb, data.Bt') #Take care of row and column
    #Ab   = data$Ab%o%data$Ab
    #Bb   = data$Bb%o%data$Bb
    Ab = kron(data.Ab, data.Ab') #Take care of row and column
    Bb = kron(data.Bb, data.Bb') #Take care of row and column
    #At   = data$At%o%data$At
    #Bt   = data$Bt%o%data$Bt
    At = kron(data.At, data.At') #Take care of row and column
    Bt = kron(data.Bt, data.Bt') #Take care of row and column

    #prc = data$tarprice
    #index_list = combinations(nrow(A),2)
    #ineqs = matrix(rep(-1000, 4*length(index_list[,1])), ncol = 4)
    prc = convert(Vector{Float64}, data.tarprice)
    temp = [Combinatorics.combinations(1:I,2)...]
    index_list = Array{Int64,2}(undef, length(temp), 2)
    for i in 1:length(temp)
        index_list[i,1] = temp[i][1]
        index_list[i,2] = temp[i][2]
    end
    #ineqs = matrix(rep(-1000, 4*length(index_list[,1])), ncol = 4)
    ineqs = Array{Float64,2}(undef,length(index_list[:,1]),4)
    comper = A + beta*Bb
    comperBB = Ab +beta*Bb
    comperTT = At +beta*Bt
    for j in 1:length(index_list[:,1])
        ineqs[j,:] = with_ineq_one_sided(comper, comperBB, comperTT, prc, index_list[j,:])
    end
    # level of ineq is too big
    #return(sum(ineqs>0))
    #display(ineqs)
    res = sum(ineqs.>0)
    return res
end

#=
maxscore_mc = function(num_agents, sd_err, num_its = 100, quietly = FALSE, withtrans = TRUE, ml = FALSE){
  ptm <- proc.time()
  myests = matrix(rep(NA, num_its*1), ncol=1)
  for(i in 1:num_its){
    # obsdat = "not yet"
    # ct = 1
    # while(obsdat == "not yet"){
    #  cat("Generating Data.  Try ", ct, "\n")
      obsdat = givemedata(num_agents, sd_err)  # was <<- not =
    #  ct = ct +1
    # }
    cat("\n Data Generation Success!\n")
    if(ml){
      est = maxLik(function(beta){mylogiter(beta, obsdat, sd_err)}, start = 2)$estimate
      if(!quietly){
        cat( "with transfer (low var) rep ", i, ": ",  as.numeric(est), " Time: ", proc.time()-ptm, "\n" )
      }
      myests[i,] = est
    }
    if(ml==FALSE){
     if(withtrans){ score_bthis = function(x){ -1*score_b_with(x, obsdat) + 100000} }
     else { score_bthis = function(x){ -1*score_b_with_one(x, obsdat) + 100000 }}
     m = findmin(score_bthis)
     if(!quietly){
        cat( "with transfer (low var) rep ", i, ": ",  as.numeric(m$optim$bestmem), " Time: ", proc.time()-ptm, "\n" )
     }
     myests[i,] = as.numeric(m$optim$bestmem)
    }
  }
  runtime = proc.time() - ptm
  res = c(sumfun(myests), runtime[1])
  myests <<- myests
  return(res)
}
=#
obsdat.tarprice
function maxscore_mc_one_sided(num_agents::Int64, sd_err::Float64;
                     num_its::Int64 = 10,
                     quietly::Bool = false,
                     withtrans::Bool = true,
                     ml::Bool = false)
    #ptm <- proc.time()
    start = Dates.unix2datetime(time())
    #myests <<- matrix(rep(NA, num_its*1), ncol=1)
    myests = Array{Float64,2}(undef, num_its*1, 1)
    #for(i in 1:num_its){
    for i in 1:num_its
       #obsdat <<- givemedata(num_agents, sd_err)
       println("Create obsdat for iteration $i \n" )
       obsdat = givemedata_one_sided(num_agents, sd_err)
       #if(ml){
       if ml == true
           #est = tryCatch(maxLik(function(beta){mylogiter(beta, obsdat)}, start = 2)$estimate, error = function(err) NA)
           @time @show est_res = Optim.optimize(beta -> mylogiter(beta, obsdat), zeros(length(beta)), LBFGS());
           est = Optim.minimizer(est_res)
           #if(!quietly){
           if quietly!== true
               #cat( "with transfer (low var) rep ", i, ": ",  as.numeric(est), " Time: ", proc.time()-ptm, "\n" )
               finish = convert(Int, Dates.value(Dates.unix2datetime(time())-start))/1000
               println("###############################################\n" )
               println("with transfer (low var) rep ", i, ": ", est, " Time: ", finish, "\n" )
               println("###############################################\n" )
               myests[i,:] = est
            end
       elseif ml!==true
           #if(withtrans){ score_bthis = function(x){ -1*score_b_with(x, obsdat) + 100000} }
           if withtrans == true
               function score_bthis_with3(beta::Vector{Float64}, obsdat::DataFrame)
                   res = -1.0*score_b_with(beta, obsdat, num_agents) + 100000.0 # need to be Float64 for bboptimize
                   println("score: $(res) \n" )
                   return res
               end
               m_res = BlackBoxOptim.bboptimize(beta -> score_bthis_with3(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta), Method = :de_rand_1_bin, MaxSteps=200)
               println("score: ",score_bthis_with3(m_res.archive_output.best_candidate, obsdat))
           #else { score_bthis = function(x){ -1*score_b(x, obsdat) + 100000 }}
           else#if withtrans == false
               function score_bthis_with_one(beta::Vector{Float64}, obsdat::DataFrame)
                   res = -1.0*score_b_with_one(beta, obsdat) + 100000.0 # need to be Float64 for bboptimize
                   println("score:$(res) \n" )
                   return res
               end
               println("******** test **********")
               m_res = BlackBoxOptim.bboptimize(beta -> score_bthis_with_one(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta), Method = :de_rand_1_bin,  MaxSteps=100)
               println("score: ", score_bthis_with_one(m_res.archive_output.best_candidate, obsdat))
               println("score of correct: ", score_b(m_res.archive_output.best_candidate, obsdat, num_agents))
               println("TRUE score of correct: ", score_b([1.5], obsdat, num_agents))
               #score_b([2.3], obsdat, 20)
           end
           #m = findmin(score_bthis)
           #@time @show m_res = Optim.optimize(beta -> score_bthis(beta, obsdat), zeros(length(beta)), LBFGS())
           #m = Optim.minimizer(m_res)
           #m_res = BlackBoxOptim.bboptimize(beta -> score_bthis(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta))
           #m_res = BlackBoxOptim.bboptimize(beta -> score_bthis(beta, obsdat); SearchRange = (0.0, 50.0), NumDimensions = length(beta), Method = :de_rand_1_bin)
           m = m_res.archive_output.best_candidate
           #if(!quietly){
           if !quietly == true
              #cat( "with transfer (low var) rep ", i, ": ",  as.numeric(m$optim$bestmem), " Time: ", proc.time()-ptm, "\n" )
              finish = convert(Int, Dates.value(Dates.unix2datetime(time())-start))/1000
              println("###############################################\n" )
              println("with transfer (low var) rep ", i, ": ", m, " \nTime: ", finish, "\n" )
              println("###############################################\n" )
              myests[i,:] = m
           end
           #myests[i,] = as.numeric(m$optim$bestmem)
           myests[i,:] = m
       end
   end
   #runtime = proc.time() - ptm
   finish = convert(Int, Dates.value(Dates.unix2datetime(time())-start))/1000
   #mean = mean(vec-meanmat)
   #sqrt = sqrt(mean((vec-meanmat)^2))
   meanmat = 1.5 # true parameter
   res_mean = mean(myests.-meanmat)
   res_sqrt = sqrt(mean((myests.-meanmat).^2))
   #res = c(sumfun(myests), runtime[1])
   #res = sumfun(myests)
   #myests <<- myests
   global myests
   @show myests
   return res_mean, res_sqrt
end

#x10 =maxscore_mc(num_agents = 100, sd_err = 1, ml = TRUE)
@time @show res_mean, res_sqrt = maxscore_mc_one_sided(10, 1.0,quietly= false,　withtrans= true,　ml= false)

x10 =maxscore_mc_one_sided(100, 1, true)
x11 =maxscore_mc(num_agents = 100, sd_err = 5, ml = TRUE)
x12 =maxscore_mc(num_agents = 100, sd_err = 20, ml = TRUE)
