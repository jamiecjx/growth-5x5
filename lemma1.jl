using IntervalArithmetic, ProgressBars, LinearAlgebra, IterTools, 
    JLD2, StaticArrays, Plots, ProgressMeter

g = Meta.parse(ARGS[1]) |> eval
p3b = inf((3 - sqrt(9 - 4*convert(Interval, g)))/2)
levels = parse(Int, ARGS[2])
function growth_factor_3(B; show=false)
    A = copy(B)
    # assume A_11 is 1
    A[2,:] -= A[1,:] * A[2,1]
    A[3,:] -= A[1,:] * A[3,1]
    # assume A_22 is maximum value
    if show
        display(A)
    end
    A[3,:] -= A[2,:] * A[3,2] / A[2,2]
    c1 = abs(A[3,3])
    return c1
end
function check_growth(A; show=false)
    # assume that A[1,1] = 1 and A[1,2], A[1,3], A[2,1], A[3,1] in [0,1] and everything else in [-1,1]
    A = copy(A)
    A[2,:] -= A[1,:] * A[2,1]
    A[3,:] -= A[1,:] * A[3,1]
    # now assume A[2,2] is the greatest magnitude entry because otherwise
    # you can pivot the original matrix to make it so
    p2 = sup(abs(A[2,2]))
    # entry must be greater than g/2
    if p2 < g/2
        return true
    end

    # from now on we know A[2,2] ≤ -g/2
    # assume now 2,2 is the greatest in magnitude
    # we discard boxes where this can never be true
    if maximum(inf.(abs.(A[2:3,2:3]))) > -inf(A[2,2])
        return true
    end
    # and now we restrict to the case where it is true
    A[2,2] = A[2,2] ∩ interval(-100, -maximum(inf.(abs.(A[2:3,2:3]))))
    A[2,3] = A[2,3] ∩ interval(-p2, p2)
    A[3,2] = A[3,2] ∩ interval(-p2, p2)
    A[3,3] = A[3,3] ∩ interval(-p2, p2)
    # # calculating last entry by finding its maximum
    # if maximise_p3(A) < g
    #     return true
    # end
    # Cohen paper tests 
    # 1: if they have different signs then p2 ≤ 2 
    if A[3,3] * A[3,2] * A[2,3] < 0
        return true
    end
    # 2: 3 separate cases left
    if A[3,2] > 0 && A[2,3] < 0 && A[3,3] < 0
        α1 = A[3,3] / -p2
        α2 = A[2,3] / -p2
        if α1 < 1/2
            return true
        end
    elseif A[3,2] < 0 && A[2,3] > 0 && A[3,3] < 0
        α1 = A[3,3] / -p2
        α2 = A[3,2] / -p2
        if α1 < 1/2
            return true
        end
    elseif A[3,2] < 0 && A[2,3] < 0 && A[3,3] > 0
        α1 = A[2,3] / -p2
        α2 = A[3,2] / -p2
        if α1 < 1/2 || α2 < 1/2
            return true
        end
    else
        α1 = NaN
        α2 = NaN
    end
    if !isnan(α1)
        α1 = α1 ∩ interval(-10,1)
        α2 = α2 ∩ interval(-10,1)
        if sup((α1 + α2) * p2 + α1 * α2 * p2 * (1-p2)) < g
            return true
        end
        if sup(α1* p2 + α2 * p2 + α1 * α2 * p2 - α1 * α2 * p2^2) < g
            return true
        end
        if α1 ⊂ interval(1/2, 1) && α2 ⊂ interval(1/4, 1)
            d1 = sup(1/4 * α1 * α2 * (1 + 1/α1 + 1/α2)^2)
            d2 = sup(1/4 * (α1 * α2 + α1 + α2) * (1 + 1/α1 + 1/α2))
            d3 = sup(1/4 * (α1 * α2 + α1 + α2)^2 / (α1 * α2))
            d4 = sup(1/4 * (α1 * α2 + α1 / α2 + α2 / α1 + 2α1 + 2α2 + 2))
            # (10) in cohen
            if min(d1,d2,d3,d4) < g
                return true
            end
            # (11) in cohen
            if sup(p2) < p3b
                return true
            end
        end
    end
    # finally, literally calculate the bounds on p3
    A[3,:] -= A[2,:] * A[3,2] / A[2,2]
    if sup(abs(A[3,3])) < g
        return true
    end
    # if you have reached here, then you haven't proven
    # that there are no growth matrices more than g
    false
end
println("Running with $(Threads.nthreads()) threads")

@inbounds @inline function check_growth_fast(A_in, local_g, local_p3b; show=false)
    # faster version of the above
    A = MMatrix{3,3,Interval{Float64}}(A_in)

    a21 = A[2,1]
    a31 = A[3,1]
    A[2,2] -= A[1,2] * a21
    A[3,2] -= A[1,2] * a31
    A[2,3] -= A[1,3] * a21
    A[3,3] -= A[1,3] * a31
    p2 = sup(abs(A[2,2]))
    if p2 < local_g/2
        return true
    end
    mx_inf_abs = max(
        inf(abs(A[2,2])), inf(abs(A[2,3])),
        inf(abs(A[3,2])), inf(abs(A[3,3]))
    )
    if mx_inf_abs > -inf(A[2,2])
        return true
    end
    A[2,2] = intersect(A[2,2], interval(-Inf, -mx_inf_abs))
    p2_interval = interval(-p2, p2)
    A[2,3] = intersect(A[2,3], p2_interval)
    A[3,2] = intersect(A[3,2], p2_interval)
    A[3,3] = intersect(A[3,3], p2_interval)
    prod_sup = sup(A[3,3] * A[3,2] * A[2,3])
    if prod_sup < 0
        return true
    end
    α1 = interval(0.0)
    α2 = interval(0.0)

    cant_tell = false

    if inf(A[3,2]) > 0 && sup(A[2,3]) < 0 && sup(A[3,3]) < 0
        α1 = A[3,3] / -p2
        α2 = A[2,3] / -p2
        if sup(α1) < 0.5
            return true
        end
    elseif sup(A[3,2]) < 0 && inf(A[2,3]) > 0 && sup(A[3,3]) < 0
        α1 = A[3,3] / -p2
        α2 = A[3,2] / -p2
        if sup(α1) < 0.5
            return true
        end
    elseif sup(A[3,2]) < 0 && sup(A[2,3]) < 0 && inf(A[3,3]) > 0
        α1 = A[2,3] / -p2
        α2 = A[3,2] / -p2
        if sup(α1) < 0.5 || sup(α2) < 0.5
            return true
        end
    else
        cant_tell = true
    end
    if !cant_tell
        α1 = intersect(α1, interval(-Inf, 1.0))
        α2 = intersect(α2, interval(-Inf, 1.0))

        expr1 = (α1 + α2) * p2 + α1 * α2 * p2 * (1.0 - p2)
        if sup(expr1) < local_g
            return true
        end

        expr2 = α1 * p2 + α2 * p2 + α1 * α2 * p2 - α1 * α2 * p2^2
        if sup(expr2) < local_g
            return true
        end

        if inf(α1) ≥ 0.5 && sup(α1) ≤ 1.0 && inf(α2) ≥ 0.25 && sup(α2) ≤ 1.0
            d1 = sup(1/4 * α1 * α2 * (1 + 1 / α1 + 1 / α2)^2)
            d2 = sup(1/4 * (α1 * α2 + α1 + α2) * (1 + 1 / α1 + 1 / α2))
            d3 = sup(1/4 * (α1 * α2 + α1 + α2)^2 / (α1 * α2))
            d4 = sup(1/4 * (α1 * α2 + α1 / α2 + α2 / α1 + 2 * α1 + 2 * α2 + 2))
            if min(d1, d2, d3, d4) < local_g
                return true
            end
            if sup(p2) < local_p3b
                return true
            end
        end
    end
    A[3,3] -= A[2,3] * A[3,2] / A[2,2]
    if sup(abs(A[3,3])) < local_g
        return true
    end
    return false
end

const VAR_POSITIONS = [(1,2),(1,3),(2,1),(2,2),(2,3),(3,1),(3,2),(3,3)]
function split_box(A, dim)
    i,j = VAR_POSITIONS[dim]
    I = A[i,j]
    m = mid(I)
    left  = copy(A); left[i,j]  = interval(inf(I), m)
    right = copy(A); right[i,j] = interval(m, sup(I))
    return (left, right)
end

function refine_boxes(boxes, dim)
    new_boxes = MMatrix{3,3}{Interval{Float64}}[]
    for A in ProgressBar(boxes)
        B1, B2 = split_box(A, dim)
        if !check_growth_fast(B1, g, p3b)
            push!(new_boxes, B1)
        end
        if !check_growth_fast(B2, g, p3b)
            push!(new_boxes, B2)
        end
    end
    return new_boxes
end

function search_box(box, g, p3b; lim=100000, stacklim=1000, show=false)
    stack = [(box,1)]
    t = 0
    while 0 < length(stack) < stacklim && t < lim
        (A,d) = pop!(stack)
        B1, B2 = split_box(A, d)
        if !check_growth_fast(B1, g, p3b)
            push!(stack, (B1, d == 8 ? 1 : d + 1))
        end
        if !check_growth_fast(B2, g, p3b)
            push!(stack, (B2, d == 8 ? 1 : d + 1))
        end
        t += 1
        if show && t % 100000 == 0
            display(length(stack))
        end
    end
    stack
end


A3 = [1 1 1/2;1 -1/2 -1;1/2 -1 1]
A31 = [1 1/2 1;1 -1 -1/2;1/2 1 -1]
A32 = [1 1 1/2;1/2 -1 1;1 -1/2 -1]


#### begin code here

println("Begin Lemma 1")
matlist = [MMatrix{3,3}([
    1..1 0..1 0..1
    0..1 -1..1 -1..1
    0..1 -1..1 -1..1
])]
matlist = refine_boxes(matlist, 4)
matlist = refine_boxes(matlist, 5)
matlist = refine_boxes(matlist, 7)
matlist = refine_boxes(matlist, 8)
println("Box width = 1")
display(length(matlist))

for k=1:levels
    for i=1:8
        global matlist = refine_boxes(matlist, i)
    end
    println("Box width = $(0.5^k)")
    display(length(matlist))
end

progress = Progress(length(matlist); showspeed=true)
Threads.@threads for i in 1:length(matlist)
    test1 = search_box(matlist[i], g, p3b; lim=10000, stacklim=1000, show=false)
    if length(test1) == 0
        matlist[i] *= 0
    end
    next!(progress)
end

filter!(x -> !iszero(x), matlist)

println("Finished")
display(length(matlist))
filename = "lemma1_$(Float64(g))_$(levels)_matlist.jld2"
JLD2.save(filename, "matlist", matlist)



# matlist2 = Vector{MMatrix{3,3,Interval{Float64}}}()
# for M in ProgressBar(matlist)
#     test1 = [M]
#     for i=1:8
#         test1 = refine_boxes(test1, i)
#     end
#     if length(test1) > 0
#         push!(matlist2, M)
#     end
# end
# length(matlist2)


# function measure_dist(matlist, distance)
#     dists = zeros(length(matlist))
#     for (i, M) in enumerate(matlist)
#         dists[i] = min(distance(M, A3), distance(M, A31), distance(M, A32))
#     end
#     dists
# end


# sort!(matlist2, by = box -> -min(
#     maximum(dist.(box, interval.(A3))),
#     maximum(dist.(box, interval.(A31))),
#     maximum(dist.(box, interval.(A32)))
# )
# )
# map(box -> min(
#     maximum(dist.(box, interval.(A3))),
#     maximum(dist.(box, interval.(A31))),
#     maximum(dist.(box, interval.(A32)))
# ), matlist2)


# checklist = fill(-1, length(matlist2));
# for i in ProgressBar(1:length(matlist2))
#     test1 = search_box(matlist2[i])
#     checklist[i] = length(test1)
#     if length(test1) > 0
#         println("fail")
#         display(i)
#         break
#     end
# end


# matlist3 = Vector{MMatrix{3,3,Interval{Float64}}}()
# for i in eachindex(matlist2)
#     if checklist[i] != 0
#         push!(matlist3, matlist2[i])
#     end
# end

# map(box -> min(
#     maximum(dist.(box, interval.(A3))),
#     maximum(dist.(box, interval.(A31))),
#     maximum(dist.(box, interval.(A32)))
# ), matlist3)













# # # ended here
# # JLD2.save("redo2_2_1_matlist3.jld2", "matlist3", matlist3)
# # matlist3 = JLD2.load("redo2_2_1_matlist3.jld2")["matlist3"]



# for i=1:8
#     display(i)
#     matlist3 = refine_boxes(matlist3, i)
# end
# sort!(matlist3, by = box -> -min(
#     maximum(dist.(box, interval.(A3))),
#     maximum(dist.(box, interval.(A31))),
#     maximum(dist.(box, interval.(A32)))
# )
# )

# matlist4 = Vector{MMatrix{3,3,Interval{Float64}}}()
# for M in ProgressBar(matlist3)
#     test1 = [M]
#     for i=1:8
#         test1 = refine_boxes(test1, i)
#     end
#     if length(test1) > 0
#         push!(matlist4, M)
#     end
# end
# length(matlist4)


# # function subdivide_box(box)
# #     halves = [bisect(I) for I in box[2:9]]
# #     subboxes = [
# #         reshape(vcat(1, [halves[d][c[d]] for d in 1:8]),3,3)
# #         for c in product((1:2 for _ in 1:8)...)
# #     ]
# #     return reshape(subboxes, 256)
# # end
# # function subdivide_box_growth(box)
# #     matlist = Vector{Matrix{Interval{Float64}}}()
# #     subboxes = subdivide_box(box)
# #     for A in subboxes
# #         if !check_growth_fast(A)
# #             push!(matlist, A)
# #         end
# #     end
# #     matlist
# # end

# # p = Progress(length(iterator); desc="Processing...", showspeed=true)

# # @Threads.threads for i in eachindex(matlist5)
# #     testlist = subdivide_box_growth(matlist5[i])
# #     # testlist = vcat(subdivide_box_growth.(testlist)...)
# #     # testlist = vcat(subdivide_box_growth.(testlist)...)
# #     # testlist = vcat(subdivide_box_growth.(testlist)...)
# #     if all(x -> length(x) == 0, testlist)
# #         matlist5[i] *= 0
# #     end
# #     next!(p)
# # end

# checklist = fill(-1, length(matlist4));
# for i in ProgressBar(1:length(matlist4))
#     test1 = search_box(matlist4[i])
#     checklist[i] = length(test1)
#     if length(test1) > 0
#         println("fail")
#         display(i)
#         break
#     end
# end

# matlist5 = Vector{MMatrix{3,3,Interval{Float64}}}()
# for i in eachindex(matlist4)
#     if checklist[i] != 0
#         push!(matlist5, matlist4[i])
#     end
# end


# sort!(matlist5, by=box -> -min(
#     maximum(dist.(box, interval.(A3))),
#     maximum(dist.(box, interval.(A31))),
#     maximum(dist.(box, interval.(A32)))
# ))


# # filter!(x -> !iszero(x), matlist5);

# filename = "v2_$(g)_matlist5.jld2"
# JLD2.save(filename, "matlist5", matlist5)




# growth_factor_3(mid.(search_box(matlist5[1])[end][1]);show=true)
# mid.(search_box(matlist5[1])[end][1]) - A32
# 1

# map(box -> min(
#     maximum(dist.(box, interval.(A3))),
#     maximum(dist.(box, interval.(A31))),
#     maximum(dist.(box, interval.(A32)))
# ), matlist5)


