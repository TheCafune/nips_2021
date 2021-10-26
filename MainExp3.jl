include("Util.jl")
include("Algorithm.jl")
include("Gen.jl")

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
    lg = open("logExp3.txt", "a")
    println(lg, ags1, " ", nv(G), " ", ne(G))

    rsed = round(Int, time() * 10000)
    for k = 1 : mk
        println(lg, "k = ", k)
        ans1 = randomSelect(Ev, k, randomSeed = rsed)
        ans2 = Optimum(G, s, Ev, k)
        ans3 = exactOpinion(G, s, Ev, k)
        ans4 = approxOpinion(G, s, Ev, k)
        println(lg, "Random Score : ", score(G, ans1, s))
        println(lg, "Optimum Score : ", score(G, ans2, s))
        println(lg, "Exact Score : ", score(G, ans3, s))
        println(lg, "Approx Score : ", score(G, ans4, s))
    end
end

runMain(ARGS[1], ARGS[2])
