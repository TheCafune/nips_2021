#!/bin/sh
for var in GrQc USgrid Erdos992 Bcspwr10 Reality PagesGovernment WikiElec Dmela HepPh Anybeat PagesCompany CondMat Gplus GemsecRO ;do
    julia -O3 MainExp1.jl $var 50,10000
done
for var in Brightkite WikiTalk Douban Citeseer MathSciNet TwitterFollows Flickr FourSquare IMDB YoutubeSnap ;do
    julia -O3 MainExp2.jl $var 50,10000
done
for var in Karate Dolphins Netscience Diseasome ;do
    julia -O3 MainExp3.jl $var 8,30
done
for var in Yeast GridWorm Erdos992 Reality ;do
    julia -O3 MainExp4.jl $var 50,10000
done