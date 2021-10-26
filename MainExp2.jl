include("Util.jl")
include("Algorithm.jl")
include("Gen.jl")

function runMain(ags1, ags2)
    # Read Graph
    G = readGraph(ags1)

    # Generate initial opinion s AND candidate edge list Ev
    buf2 = split(ags2, ',')
    k = parse(Int, buf2[1])
    evn = parse(Int, buf2[2])
    s = Uniform(nv(G))
    Ev = generateEv(G, evn)

    # do experiment
    lg = open("logExp2.txt", "a")
    println(lg, ags1, " ", nv(G), " ", ne(G))

    approxTime = time()
    ans2 = approxOpinion(G, s, Ev, k, eps=0.5)
    approxTime = time() - approxTime

    println(lg, "Approx Time : ", approxTime)
    println(lg)
end

runMain(ARGS[1], ARGS[2])
