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
    for A in boxes
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