include("Util.jl")

using LightGraphs
using Random
using Laplacians

function randomSelect(Ev, k; randomSeed = round(Int, time() * 10000))
    Random.seed!(randomSeed)
    y = randperm(size(Ev, 1))
    S = Array{Tuple{Int32, Int32}, 1}()
    for i = 1 : k
        push!(S, Ev[y[i]])
    end
    return S
end

function exactOpinion(G, s, Ev, k)
    W = getW(G)
    S = Array{Tuple{Int32, Int32}, 1}()
    sev = size(Ev, 1)
    cho = zeros(Bool, sev)
    for i = 1 : k
        dz = zeros(sev)
        q = W * s
        for j = 1 : sev
            if cho[j]
                continue
            end
            (u, v) = Ev[j]
            dz[j] = (q[u] - q[v])^2 / (1.0 + W[u, u] + W[v, v] - 2*W[u, v])
        end
        xx = argmax(dz)
        cho[xx] = true
        (su, sv) = Ev[xx]
        W = updateW(W, su, sv)
        push!(S, Ev[xx])
    end
    return S
end

function approxOpinion(G, s, Ev, k; eps = 0.3)
    IpL = getSparseIpL(G)
    B = getSparseB(G)
    m = ne(G)
    n = nv(G)
    kkk = round(Int, 0.5 * log2(n) / (eps^2))
    S = Array{Tuple{Int32, Int32}, 1}()
    sev = size(Ev, 1)
    cho = zeros(Bool, sev)
    for rep = 1 : k
        dzfm = zeros(sev)
        f = approxchol_sddm(IpL, tol=1e-12)
        for i = 1 : kkk
            yy1 = B' * randn(m)
            yy2 = randn(n)
            zz1 = f(yy1)
            zz2 = f(yy2)
            for j = 1 : sev
                (uu, vv) = Ev[j]
                dzfm[j] += ((zz1[uu] - zz1[vv])^2 + (zz2[uu] - zz2[vv])^2)
            end
        end
        dzfm ./= kkk
        dzfm .+= 1.0
        q = f(s)

        dz = zeros(sev)
        for j = 1 : sev
            if cho[j]
                continue
            end
            (u, v) = Ev[j]
            dz[j] = (q[u] - q[v])^2 / dzfm[j]
        end
        xx = argmax(dz)
        cho[xx] = true
        (su, sv) = Ev[xx]
        IpL[su, su] += 1
        IpL[sv, sv] += 1
        IpL[su, sv] = -1
        IpL[sv, su] = -1
        push!(S, Ev[xx])
    end
    return S
end

function Optimum(G, s, Ev, k)
    W = getW(G)
    sev = size(Ev, 1)
    cho = zeros(Bool, sev)
    bst = 1e9
    optsn = zeros(Int, k)
    a = zeros(Int, k)
    stkW = []
    push!(stkW, W)
    top = 1

    dfs(dep) = begin
        if dep > k
            tmp = s' * stkW[top] * s
            if tmp < bst
                bst = tmp
                foreach(i -> optsn[i] = a[i], 1 : k)
            end
            return nothing
        end
        st = 1
        if dep > 1
            st = a[dep-1]
        end
        for i = st : sev
            if !cho[i]
                cho[i] = true
                (su, sv) = Ev[i]
                nW = copy(stkW[top])
                nW = updateW(nW, su, sv)
                top += 1
                push!(stkW, nW)
                a[dep] = i
                dfs(dep+1)
                pop!(stkW)
                top -= 1
                cho[i] = false
            end
        end
        return nothing
    end

    dfs(1)
    S = Array{Tuple{Int32, Int32}, 1}()
    foreach(i -> push!(S, Ev[optsn[i]]), 1 : k)
    return S
end

function DegSum(G, s, Ev, k)
    deg = degree(G)
    S = Array{Tuple{Int32, Int32}, 1}()
    sev = size(Ev, 1)
    cho = zeros(Bool, sev)
    for i = 1 : k
        dz = zeros(sev)
        for j = 1 : sev
            if cho[j]
                continue
            end
            (u, v) = Ev[j]
            dz[j] = deg[u] + deg[v]
        end
        xx = argmax(dz)
        cho[xx] = true
        (su, sv) = Ev[xx]
        deg[su] += 1
        deg[sv] += 1
        push!(S, Ev[xx])
    end
    return S
end

function DegProduct(G, s, Ev, k)
    deg = degree(G)
    S = Array{Tuple{Int32, Int32}, 1}()
    sev = size(Ev, 1)
    cho = zeros(Bool, sev)
    for i = 1 : k
        dz = zeros(sev)
        for j = 1 : sev
            if cho[j]
                continue
            end
            (u, v) = Ev[j]
            dz[j] = deg[u] * deg[v]
        end
        xx = argmax(dz)
        cho[xx] = true
        (su, sv) = Ev[xx]
        deg[su] += 1
        deg[sv] += 1
        push!(S, Ev[xx])
    end
    return S
end

function EdgeBetweenness(G0, Ev)
    G = deepcopy(G0)
    idx = Dict{Tuple{Int32, Int32}, Int32}()
    ID = 0
    for (u, v) in Ev
        ID += 1
        idx[(u, v)] = ID
        idx[(v, u)] = ID
        add_edge!(G, u, v)
    end
    Cb = zeros(size(Ev, 1))
    n = nv(G)
    m = ne(G)

    p = Array{Array{Int32, 1}, 1}(undef, n)
    d = zeros(Int32, n)
    S = zeros(Int32, n+10)
    sigma = zeros(n)
    d = zeros(Int32, n)
    Q = zeros(Int32, n+10)
    delta = zeros(n)

    for s = 1 : n
        foreach(i -> p[i] = [], 1 : n)
        top = 0
        sigma .= 0
        sigma[s] = 1.0
        d .= -1
        d[s] = 0
        front = 1
        rear = 1
        Q[1] = s

        while front <= rear
            v = Q[front]
            front += 1
            top += 1
            S[top] = v
            for w in neighbors(G, v)
                if d[w] < 0
                    rear += 1
                    Q[rear] = w
                    d[w] = d[v] + 1
                end
                if d[w] == (d[v] + 1)
                    sigma[w] += sigma[v]
                    push!(p[w], v)
                end
            end
        end

        delta .= 0

        while top > 0
            w = S[top]
            top -= 1
            for v in p[w]
                delta[v] += ((sigma[v] / sigma[w]) * (1 + delta[w]))
                if haskey(idx, (v, w))
                    Cb[idx[(v, w)]] += ((sigma[v] / sigma[w]) * (1 + delta[w]))
                end
            end
        end

    end

    return Cb
end

function topBetweenness(G, s, Ev, k)
    eb = EdgeBetweenness(G, Ev)
    S = Array{Tuple{Int32, Int32}, 1}()
    sev = size(Ev, 1)
    cho = zeros(Bool, sev)
    for i = 1 : k
        dz = zeros(sev)
        for j = 1 : sev
            if cho[j]
                continue
            end
            dz[j] = eb[j]
        end
        xx = argmax(dz)
        cho[xx] = true
        (su, sv) = Ev[xx]
        push!(S, Ev[xx])
    end
    return S
end