using LightGraphs
using Random

function powerLaw(n; alp = 2.5, xmin = 1)
    Random.seed!(round(Int, time() * 10000))
    x = rand(n)
    for i = 1 : n
        x[i] = xmin * ((1 - x[i])^(-1.0/(alp - 1.0)))
    end
    xm = x[argmax(x)]
    x ./= xm
    return x
end

function Uniform(n)
    Random.seed!(round(Int, time() * 10000))
    x = rand(n)
    return x
end

function Exponential(n; lmd = 1, xmin = 1)
    Random.seed!(round(Int, time() * 10000))
    x = rand(n)
    for i = 1 : n
        x[i] = xmin - (1.0/lmd)*log(1-x[i])
    end
    xm = x[argmax(x)]
    x ./= xm
    return x
end

function generateEv(G, evn)
    Random.seed!(round(Int, time() * 10000))
    Ev = Array{Tuple{Int32, Int32}, 1}()
    n = nv(G)
    g2 = SimpleGraph(n)
    while size(Ev, 1) < evn
        ru = rand(1:n)
        rv = rand(1:n)
        if (ru != rv) && (!has_edge(G, ru, rv)) && (!has_edge(g2, ru, rv))
            add_edge!(g2, ru, rv)
            push!(Ev, (ru, rv))
        end
    end
    return Ev
end

function score(G, addE, s)
    gt = deepcopy(G)
    for (u, v) in addE
        add_edge!(gt, u, v)
    end
    W = getW(gt)
    return s' * W * s
end
