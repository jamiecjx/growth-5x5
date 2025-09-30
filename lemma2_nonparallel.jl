using JLD2, IntervalArithmetic, ProgressBars, StaticArrays, Base.Threads, DelimitedFiles, ProgressMeter
println("Running with $(Threads.nthreads()) threads")
println("Begin Lemma 2")
include("growth_check.jl")

@inline function calc_range(box, p2, p3, M, r)
    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12 = box

    v11 = x1 * x4 + p2 * x7 * x10 - p3 * M[1,1]
    if v11.lo > r || v11.hi < -r; return true; end

    v12 = x1 * x5 + p2 * x7 * x11 - p3 * M[1,2]
    if v12.lo > r || v12.hi < -r; return true; end

    v13 = x1 * x6 + p2 * x7 * x12 - p3 * M[1,3]
    if v13.lo > r || v13.hi < -r; return true; end

    v21 = x2 * x4 + p2 * x8 * x10 - p3 * M[2,1]
    if v21.lo > r || v21.hi < -r; return true; end

    v22 = x2 * x5 + p2 * x8 * x11 - p3 * M[2,2]
    if v22.lo > r || v22.hi < -r; return true; end

    v23 = x2 * x6 + p2 * x8 * x12 - p3 * M[2,3]
    if v23.lo > r || v23.hi < -r; return true; end

    v31 = x3 * x4 + p2 * x9 * x10 - p3 * M[3,1]
    if v31.lo > r || v31.hi < -r; return true; end

    v32 = x3 * x5 + p2 * x9 * x11 - p3 * M[3,2]
    if v32.lo > r || v32.hi < -r; return true; end

    v33 = x3 * x6 + p2 * x9 * x12 - p3 * M[3,3]
    if v33.lo > r || v33.hi < -r; return true; end

    return false
end
g = Meta.parse(ARGS[1]) |> eval
p3b = inf((3 - sqrt(9 - 4*convert(Interval, g)))/2)

function subdivide_box(box, dim)
    I1, I2 = bisect(box[dim], 0.5)
    return (
        setindex(box, I1, dim),
        setindex(box, I2, dim),
    )
end

matlist = JLD2.load(ARGS[3])["matlist"]

A3 = [1 1 1/2;1 -1/2 -1;1/2 -1 1]
A31 = [1 1/2 1;1 -1 -1/2;1/2 1 -1]
A32 = [1 1 1/2;1/2 -1 1;1 -1/2 -1]

function test_if_m(mat, d, p2v, p3v; show=false, lim=400, timeout=200000000)
    boxlist = [(SVector{12}([-1..1 for _=1:12]), 1)]
    t = 0
    println("starting test_if_m round")
    display(show)
    while 0 < length(boxlist) < lim
        box, dim = pop!(boxlist)
        b1, b2 = subdivide_box(box, dim)
        if !calc_range(b1, p2v, p3v, mat, d)
            push!(boxlist, (b1, dim == 12 ? 1 : dim + 1))
        end
        if !calc_range(b2, p2v, p3v, mat, d)
            push!(boxlist, (b2, dim == 12 ? 1 : dim + 1))
        end
        t += 1
        if show && t % 1000000 == 0
            println("$(t/1000000), $(length(boxlist)), $(diam(boxlist[end][1][1]))")
        end
        if t == timeout
            return boxlist 
        end
    end
    boxlist
end


function lemma2_test(mat, d, p2v, p3v; show=false)
    test1 = length(test_if_m(mat, 1, p2v, p3v; show=show))
    if test1 == 0
        return 0
    else
        if show
            println("refinement")
        end
        # more refined test: essentially repeating lemma 1 at a smaller scale
        tempmatlist = [mat]
        for k=1:8
            tempmatlist = refine_boxes(tempmatlist, k)
        end
        for k in eachindex(tempmatlist)
            if show
                println("k = $k")
            end
            test1 = search_box(tempmatlist[k], g, p3b; lim=100000, stacklim=1000)
            if length(test1) == 0
                tempmatlist[k] *= 0
            end
        end
        filter!(x -> !iszero(x), tempmatlist)
        t = 0
        for M in tempmatlist
            t += length(test_if_m(M, d, p2v, p3v; show=show)) > 0
        end
        return t
    end  
end


p3v = Meta.parse(ARGS[2]) |> eval
p2v = sup(1.5 + sqrt(9/4 - convert(Interval, p3v)))




n = length(matlist)
validlist = fill(-1, n)
logfile = "lemma2_" * ARGS[3] * "_$(Float64(p3v))_1.csv"

if isfile(logfile)
    println("Reloading previous results from $logfile ...")
    data = readdlm(logfile, ',', Int)
    for row in eachrow(data)
        idx, val = row
        validlist[idx] = val
    end
end


iterator = [x for x in eachindex(matlist) if validlist[x] == -1]
println("Remaining to compute: ", length(iterator), " out of $n")


open(logfile, "a") do io
    for j in eachindex(iterator)
        i = iterator[j]
        println("i = $i")
        result = lemma2_test(matlist[i], 1, p2v, p3v; show=true)
        validlist[i] = result
        println(io, "$i,$result")
        flush(io)
    end
end

println("Finished")


# for k in 1:length(tempmatlist)
#     display(k)
#     test1 = search_box(tempmatlist[k], g, p3b; lim=100000, stacklim=1000, show=false)
#     if length(test1) == 0
#         tempmatlist[k] *= 0
#     end
# end
