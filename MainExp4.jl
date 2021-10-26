include("Util.jl")
include("Algorithm.jl")
include("Gen.jl")

function doExp(G, s, Ev, k, lg; rsed = round(Int, time() * 10000))
    println(lg, "k = ", k)
    println(lg, "Random Score : ", score(G, randomSelect(Ev, k, randomSeed = rsed), s))
    println(lg, "Exact Score : ", score(G, exactOpinion(G, s, Ev, k), s))
    println(lg, "Approx Score : ", score(G, approxOpinion(G, s, Ev, k), s))
    println(lg, "Betweenness Score : ", score(G, topBetweenness(G, s, Ev, k), s))
    println(lg, "DegSum Score : ", score(G, DegSum(G, s, Ev, k), s))
    println(lg, "DegProduct Score : ", score(G, DegProduct(G, s, Ev, k), s))
end

function runMain(ags1, ags2)
    # Read Graph
    G = readGraph(ags1)

    # Generate initial opinion s AND candidate edge list Ev
    buf2 = split(ags2, ',')
    mk = parse(Int, buf2[1])
    evn = parse(Int, buf2[2])
    s = Uniform(nv(G))
    Ev = generateEv(G, evn)

    # do experiment
    lg = open("logExp4.txt", "a")
    println(lg, ags1, " ", nv(G), " ", ne(G))

    prsed = round(Int, time() * 10000)
    foreach(i -> doExp(G, s, Ev, i, lg, rsed = prsed), 1 : mk)
end

runMain(ARGS[1], ARGS[2])
