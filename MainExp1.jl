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
    lg = open("logExp1.txt", "a")
    println(lg, ags1, " ", nv(G), " ", ne(G))

    exactTime = time()
    ans1 = exactOpinion(G, s, Ev, k)
    exactTime = time() - exactTime

    approxTime = time()
    ans2 = approxOpinion(G, s, Ev, k, eps=0.5)
    approxTime = time() - approxTime

    score0 = score(G, [], s)
    score1 = score(G, ans1, s)
    score2 = score(G, ans2, s)

    println(lg, "Exact Time : ", exactTime)
    println(lg, "Approx Time : ", approxTime)
    println(lg, "Exact score : ", score1)
    println(lg, "Approx score : ", score2)
    println(lg, "Exact delta : ", score0-score1)
    println(lg, "Approx delta : ", score0-score2)
    println(lg, "Ratio : ", (score0-score2) / (score0-score1))
    println(lg)
end

runMain(ARGS[1], ARGS[2])
