using LightGraphs
using SparseArrays
using LinearAlgebra

function readGraph(gn; dataPath = "networkData/") 
    return loadgraph(dataPath*gn*".lgz")
end

function getL(G)
    n = nv(G)
	L = zeros(n, n)
	for x in collect(edges(G))
        (u, v) = Tuple(x)
		L[u, u] += 1
		L[v, v] += 1
		L[u, v] -= 1
		L[v, u] -= 1
	end
	return L
end

function getW(G)
	return inv(getL(G)+I(nv(G)))
end

function getSparseIpL(G)
    n = nv(G)
    m = ne(G)
    d = ones(n)
    for x in collect(edges(G))
        (u, v) = Tuple(x)
        d[u] += 1
        d[v] += 1
    end
    Is = zeros(Int32, m*2+n)
    Js = zeros(Int32, m*2+n)
    Vs = zeros(m*2+n)
    ID = 0
    for x in collect(edges(G))
        (u, v) = Tuple(x)
        ID += 1
        Is[ID] = u
        Js[ID] = v
        Vs[ID] = -1
        Is[ID + m] = v
        Js[ID + m] = u
        Vs[ID + m] = -1
    end
    for i = 1 : n
        Is[m + m + i] = i
        Js[m + m + i] = i
        Vs[m + m + i] = d[i]
    end
    return sparse(Is, Js, Vs, n, n)
end

function getSparseL(G)
    n = nv(G)
    m = ne(G)
    d = zeros(n)
    for x in collect(edges(G))
        (u, v) = Tuple(x)
        d[u] += 1
        d[v] += 1
    end
    Is = zeros(Int32, m*2+n)
    Js = zeros(Int32, m*2+n)
    Vs = zeros(m*2+n)
    ID = 0
    for x in collect(edges(G))
        (u, v) = Tuple(x)
        ID += 1
        Is[ID] = u
        Js[ID] = v
        Vs[ID] = -1
        Is[ID + m] = v
        Js[ID + m] = u
        Vs[ID + m] = -1
    end
    for i = 1 : n
        Is[m + m + i] = i
        Js[m + m + i] = i
        Vs[m + m + i] = d[i]
    end
    return sparse(Is, Js, Vs, n, n)
end

function getSparseB(G)
    n = nv(G)
    m = ne(G)
    Is = zeros(Int32, m*2)
    Js = zeros(Int32, m*2)
    Vs = zeros(m*2)
    ID = 0
    for x in collect(edges(G))
        (u, v) = Tuple(x)
        ID += 1
        Is[ID] = ID
        Js[ID] = u
        Vs[ID] = 1
        Is[ID + m] = ID
        Js[ID + m] = v
        Vs[ID + m] = -1
    end
    return sparse(Is, Js, Vs, m, n)
end

function updateW(W, u, v)
    c = W[:,u] - W[:,v]
    fac = 1.0 / (1.0 + c[u] - c[v])
    return W - fac * c * c'
end