Rabosky et al. 2007. Proc. Roy. Soc.

L(pendant) = sum(log(1-Bi)) + sum((ni - 1)*log(B))

where B = (exp(r*t.ip)-1)/(exp(r*t.ip) - epsilon) <- this is equivalent to Foote B in terms of b and d (I checked)
	t.ip is the length of the ith pendant edge

L(internal) = N.in * log(r) - r * sum(t.in) - sum(log(1 - epsilon * exp(-r*x.in)))
	N.in = 3 of internal edges; t.in = internal edge lengths;  x.in = birth time of ith lineage

L(total) = N.int * log(r) - r * sum(t.in) - sum(log(1 - epsilon * exp(-r*xi))) + sum(log(1-Bi)) + sum((ni - 1)*log(B))
	This seems to differ from Nee et al. 1994. Proc. Roy. Soc. and Paradis. 2003. Proc. Roy. Soc.

Consider the situation where all n tips represent a single species, and epsilon = zero:

L(pb) = N.int * log(r) - r * sum(t.in) - 0 + sum(log(1-Bi)) + 0
	log(1-Bi) simplifies to -r * sum(t.ip). so:

L(pb) = N.int * log(r) - r * (sum (t.in + t.ip)) <- the last term obviously representing the sum of all edge lengths
	I have confirmed that this identical to the equation in Sanderson & Wojciechowski. 1996. Am. J. Bot. (which Dan cites).

Look okay? I am not so sure. Likelihoods from MEDUSA, yule (APE), and pureBirth (LASER) are crazy different. For example, for the following tree, here is the output from turboMEDUSA:

tr <- read.tree(text = "((Blastocerus_dichotomus:4.573297053,((Mazama_temama:0.2777451736,(Mazama_pandora:0.2436685108,Mazama_americana:0.2436685108):0.03407666275):0.2116110583,(Mazama_nana:0.3329670708,((Mazama_rufina:0.05007015173,Mazama_gouazoubira:0.05007015173):0.1428284398,(Mazama_bricenii:0.0911496959,Mazama_chunyi:0.0911496959):0.1017488957):0.1400684793):0.1563891611):4.083940821):2.147809152,((Hippocamelus_antisensis:1.850706343,Hippocamelus_bisulcus:1.850706343):3.292789896,Ozotoceros_bezoarticus:5.143496239):1.577609966);")

  N.Models Shift.Node Cut.at       r epsilon LnLik.part
1        1         NA     NA 0.37471      NA    -19.816

Model fit summary for model #1:

	Log-likelihood = -19.81613
	aicc = 41.82273

and here is the output from yule in APE:

$lambda
[1] 0.3747063

$loglik
[1] -2.313819

pureBirth agrees with APE (which surprised me, as I am using Raboskys equations above). The difference get only greater with larger trees. For whales:

turboMEDUSA = -264.3723, yule = pureBirth = 22.52087

The paramter estimates are the same, but I am worried about the likelihoods. After re-reading some of the original papers, it now seems silly to even run optimize for Yule (at least when clades do not yet contain a break, and tips represent a single species), as MLEs are available (below; although the likelihood equations are often not given). This would make MEDUSA about a trillion times faster for the for the first few iterations. Regardless, they provide a good starting point for optimizing the Rabosky equations. However, I am confused over which MLE is correct.

lambda = (Num.nodes-2)/sum(time) i.e. the number of speciation events over the summed branch lengths. 

Or, is considering the equation Nt = No*exp(r*t):

lambda = log(n)/t.depth				[when including a stem age]
lambda = (log(n) - log(2))/t.depth	[when not; From Magallon and Sanderson. 2001. Evolution. They also offer a method-of-moments estimator for estimating r in a birth death model, but have no corresponding estimate of epsilon]

First, these do not work for 1 or two (the second equation) taxa. Second, it seems difficult to reconcile the first MLE with these last two, as the first incorporates edge length information (which I would think is integral). However, the third equation gives similar paramter estimates to pureBirth, yule, and MEDUSA. Also, I would be interested to know what assumptions, if any, these last two equation make about the distribution of edge lengths. It is apparently not equal edge lengths (I checked).

Although APE yule cites the first equation above as the MLE, it instead uses the following:

loglik <- -lambda * sum(all edges) + lfactorial(Num.nodes) + Num.nodes-1 * log(lambda)
	I am not sure exactly where this comes from; there are a number of very similar equations scattered among several papers. [Paradis. 2003. Proc. Roy. Soc.]
	This is identical to the Rabosky equation I am using above except for the middle term and one fewer nodes muliplied by log(lambda) in the third term.

I am not exactly sure what Dan is doing in pureBirth, as the code is very difficult to read and he does not provide a citation.

Fortunately, MEDUSA recovers this parameter estimate, but again with the much larger likelihood.

So, finally to the questions (I apologize for how long this is).

1. How do I know which MLE estimator is best? It seems reasonable to use Rabosky for everything, rather than mixing methods (especially when the likelihoods differ so much, and particularly when considering both yule and bd in the same analysis).

2. Now for the more serious question that I have been agonizing over all day. For pendant edges (without fossils) I have been using the Rabosky pendant equations above. These liklihoods agree with the FitzJohn code. For fossils I have been using the Foote equations. However, I happened to glance at some of your old code today and see that you were using Foote for pendant edges even with extant taxa alone. While I have been consistent, I am worried that this invalidates the method. Does it?