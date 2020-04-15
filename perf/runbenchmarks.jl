using Rotations
using BenchmarkTools
import Base.Iterators: product
import Random

const T = Float64

const suite = BenchmarkGroup()
Random.seed!(1)

noneuler = suite["Non-Euler conversions"] = BenchmarkGroup()
rotationtypes = [RotMatrix3{T}, Quat{T}, UnitQuaternion{T}, SPQuat{T}, MRP{T}, AngleAxis{T}, RotationVec{T}]
for (from, to) in product(rotationtypes, rotationtypes)
    if from != to
        name = "$(string(from)) -> $(string(to))"
        # use eval here because of https://github.com/JuliaCI/BenchmarkTools.jl/issues/50#issuecomment-318673288
        noneuler[name] = eval(:(@benchmarkable $to(rot) setup = rot = rand($from)))
    end
end

euler = suite["Euler conversions"] = BenchmarkGroup()
eulertypes = [
    RotX{T}, RotY{T}, RotZ{T},
    RotXY{T}, RotYX{T}, RotZX{T}, RotXZ{T}, RotYZ{T}, RotZY{T},
    RotXYX{T}, RotYXY{T}, RotZXZ{T}, RotXZX{T}, RotYZY{T}, RotZYZ{T},
    RotXYZ{T}, RotYXZ{T}, RotZXY{T}, RotXZY{T}, RotYZX{T}, RotZYX{T}]
for from in eulertypes
    to = RotMatrix3{T}
    name = "$(string(from)) -> $(string(to))"
    # use eval here because of https://github.com/JuliaCI/BenchmarkTools.jl/issues/50#issuecomment-318673288
    euler[name] = eval(:(@benchmarkable $to(rot) setup = rot = rand($from)))
end


composition = suite["Composition"] = BenchmarkGroup()
# use eval here because of https://github.com/JuliaCI/BenchmarkTools.jl/issues/50#issuecomment-318673288
composition["RotMatrix{3} * RotMatrix{3}"] = eval(:(@benchmarkable r1 * r2 setup = begin
    r1 = rand(RotMatrix3{T})
    r2 = rand(RotMatrix3{T})
end))

paramspath = joinpath(@__DIR__, "benchmarkparams.json")
if isfile(paramspath)
    loadparams!(suite, BenchmarkTools.load(paramspath)[1], :evals, :samples);
else
    tune!(suite, verbose = true)
    BenchmarkTools.save(paramspath, BenchmarkTools.params(suite))
end

results = run(suite, verbose=true)
println()
for (groupname, groupresults) in results
    println("Group: $groupname")
    max_name_length = maximum(length, keys(groupresults))
    for (name, trial) in groupresults
        if minimum(trial).allocs != 0
            @warn("$name  allocates!")
        end
        println(rpad(name, max_name_length), " ", BenchmarkTools.prettytime(minimum(trial).time))
    end
    println()
end
