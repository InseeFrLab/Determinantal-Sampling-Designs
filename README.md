Introduction
=============

In the sequel, we provide some R-tools to help implement Determinantal Sampling Designs such as described in [Loonis and Mary (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0378375818300533) or in [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux). The tools consist of the translation in $R$ of tools that were written at first in _SAS_. The translation is due to Kim Antunez and Fabrice Nathan Tchazou Kamwa from _Ensae_, Jeanne Monnier from _Ecole Centrale de Lyon_, and Loik Acakpo Addra from _Gustave Eiffel University_. So far, all the programs are experimental. It means that they were strongly tested, but it may remains some bugs. Moreover, they were used first and foremost to help carry out the simulation studies described in [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux). Finaly, none of us is a professionnal coder, so the code itself might be improved. Should you have any questions or remarks, please feel free to contact Vincent Loonis at vincent.loonis (a)insee.fr . 

1. `Ppi` provides a real matrix $\Phi^N$, whose size is $(n\times N)$, and such that ${\Phi^N}^\intercal\Phi^N=P^\Pi$, where $P^\Pi$ is the projection matrix associated to $\Pi$, such as in described [Loonis and Mary (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0378375818300533).
2. `CaDsd` stands for _Contructing All Determintal Designs_. For a given vector $\Pi$,  `CaDsd` will provide __almost__ all the hermitian matrices whose diagonal is $\Pi$. `CaDsd` implements [Fickus and Al. (2013)]([https://www.sciencedirect.com/science/article/pii/S1063520312001297](https://link.springer.com/chapter/10.1007/978-0-8176-8373-3_2)) algorithm, according to the reparametrisation due to [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux).
3. `Drawing_Dsd` will draw one or several samples from a random variable $\mathbb{S}$ whose law is a determinantal sampling design and whose kernel is a __projection__ matrix $K$. `Drawing_Dsd` relies on an algorithm introduced by [ Lavancier et collab., (2015)](https://www.jstor.org/stable/24775312#metadata_info_tab_contents), and presented in both [Loonis and Mary (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0378375818300533) and [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux).



`Ppi`
=============

Let $\Pi$ be a vector in $\]0,1\[^N$, such that $\overset{N}{\underset{k=1}\sum}\Pi_k=n\in \mathbb{N}$. Ppi($\Pi$) will construct an orthomormal eigenbasis for the projection matrix $P^\Pi$ that was introduced by [Loonis and Mary (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0378375818300533). The construction relies on the _fast_ algorithm available in [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux) (Algorithm 7.1).  We give below an example for $\Pi=(0.5,0.75,0.75,0.2,0.4,0.6,0.8)^\intercal$. The first matrix is the eigen basis, the second is $K=P^\Pi$. 

```
> PhiNT<-Ppi(c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8))
> print(PhiNT)
          [,1]       [,2]      [,3]       [,4]
[1,] 0.7071068  0.0000000 0.0000000  0.0000000
[2,] 0.5000000 -0.7071068 0.0000000  0.0000000
[3,] 0.5000000  0.7071068 0.0000000  0.0000000
[4,] 0.0000000  0.0000000 0.4472136  0.0000000
[5,] 0.0000000  0.0000000 0.6324555  0.0000000
[6,] 0.0000000  0.0000000 0.5163978 -0.5773503
[7,] 0.0000000  0.0000000 0.3651484  0.8164966

> K<-PhiNT%*%t(PhiNT)
> print(K)
          [,1]       [,2]       [,3]      [,4]      [,5]       [,6]       [,7]
[1,] 0.5000000  0.3535534  0.3535534 0.0000000 0.0000000  0.0000000  0.0000000
[2,] 0.3535534  0.7500000 -0.2500000 0.0000000 0.0000000  0.0000000  0.0000000
[3,] 0.3535534 -0.2500000  0.7500000 0.0000000 0.0000000  0.0000000  0.0000000
[4,] 0.0000000  0.0000000  0.0000000 0.2000000 0.2828427  0.2309401  0.1632993
[5,] 0.0000000  0.0000000  0.0000000 0.2828427 0.4000000  0.3265986  0.2309401
[6,] 0.0000000  0.0000000  0.0000000 0.2309401 0.3265986  0.6000000 -0.2828427
[7,] 0.0000000  0.0000000  0.0000000 0.1632993 0.2309401 -0.2828427  0.8000000
> 
```


`CaDsd`
=============

Let $\Pi$ a vector in $\]0,1\[^N$, such that $\underset{k=1}{\overset{N}\sum}\Pi^k=\mu \in \mathbb{R}^+$. Let $M\leq \mu$ be an integer and $\Omega$ and $\rho$ two matrices whose size is respectively $M(\times N)$ and $M(\times (N-1))$ and whose entries lie in $\[0,1\]$. $CaDsd(\Pi,M,\Omega,\rho)$ provides a list that consists of : 

1. the hemitian matrix $K^{\Pi^\triangleright}(M,\Omega,\rho)$ such as described in  [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux), where the entries of $\Pi^\triangleright$ are those of $\Pi$ in descending order.
2. the sequence of spectra of all $K^{\Pi^\triangleright}(M,\Omega,\rho)$'s principal submatrices
3. an eigenbasis for $K^{\Pi^\triangleright}(M,\Omega,\rho)$.   

In the following example, we consider $\Pi=(0.5,0.75,0.75,0.2,0.4,0.6,0.8)^\intercal$ and observe that $K^{\Pi^\triangleright}(M,0^{(4 \times 7)},0^{(4\times 6})=H_N\odot P^{\Pi^\triangleright}$, where $\odot$ stands for the Hadamard product and $H_N$ is an hermitian matrix, whose 
diagonal entries are $1$ and off-diagonal entries's modulus is $1$.
```
> omega=matrix(0,4,7)
> rho=matrix(0,4,6)
> K<-CaDsd(omega,rho,M=4,pi=c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8))
> print(K$spectrum)
     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
[1,]  0.0 0.00  0.0  0.0  0.4  0.8    1
[2,]  0.0 0.00  0.3  0.9  1.0  1.0    1
[3,]  0.0 0.55  1.0  1.0  1.0  1.0    1
[4,]  0.8 1.00  1.0  1.0  1.0  1.0    1
> print(K$EigenBasis)
              [,1]             [,2]             [,3]             [,4]
phi -0.07968191+0i  2.817181e-01+0i -3.984095e-01+0i  7.453560e-01+0i
     0.08908709+0i -3.149704e-01+0i  4.454354e-01+0i  6.666667e-01+0i
    -0.08908708+0i  3.149704e-01+0i  8.017836e-01+0i  1.211648e-09+0i
    -0.21081850+0i  7.453560e-01+0i -2.199235e-08+0i -1.228836e-09+0i
     0.57735031+0i  4.082484e-01+0i -6.692061e-09+0i -3.739229e-10+0i
     0.63245558+0i -1.250641e-08+0i  4.887191e-09+0i  2.730747e-10+0i
     0.44721361+0i -8.843369e-09+0i  3.455766e-09+0i  1.930930e-10+0i
> print(K$K)
              phi                                                                                    
phi  0.8000000+0i  0.2236068+0i -0.2236068+0i  0.2267787+0i  0.0690066+0i -0.0503953+0i -0.0356348+0i
     0.2236068+0i  0.7500000+0i  0.2500000+0i -0.2535463+0i -0.0771517+0i  0.0563436+0i  0.0398410+0i
    -0.2236068+0i  0.2500000+0i  0.7499999+0i  0.2535462+0i  0.0771517+0i -0.0563436+0i -0.0398410+0i
     0.2267787+0i -0.2535463+0i  0.2535462+0i  0.5999999+0i  0.1825742+0i -0.1333333+0i -0.0942809+0i
     0.0690066+0i -0.0771517+0i  0.0771517+0i  0.1825742+0i  0.5000001+0i  0.3651484+0i  0.2581989+0i
    -0.0503953+0i  0.0563436+0i -0.0563436+0i -0.1333333+0i  0.3651484+0i  0.4000000+0i  0.2828427+0i
    -0.0356348+0i  0.0398410+0i -0.0398410+0i -0.0942809+0i  0.2581989+0i  0.2828427+0i  0.2000000+0i
> V<-Ppi(c(0.8, 0.75, 0.75, 0.6, 0.5, 0.4, 0.2))
> print(V%*%t(V))
           [,1]        [,2]        [,3]       [,4]        [,5]        [,6]        [,7]
[1,] 0.80000000  0.22360680  0.22360680  0.2267787  0.06900656  0.05039526  0.03563483
[2,] 0.22360680  0.75000000 -0.25000000 -0.2535463 -0.07715167 -0.05634362 -0.03984095
[3,] 0.22360680 -0.25000000  0.75000000 -0.2535463 -0.07715168 -0.05634362 -0.03984095
[4,] 0.22677869 -0.25354628 -0.25354628  0.6000000  0.18257418  0.13333333  0.09428090
[5,] 0.06900656 -0.07715167 -0.07715168  0.1825742  0.50000000 -0.36514837 -0.25819889
[6,] 0.05039526 -0.05634362 -0.05634362  0.1333333 -0.36514837  0.40000000  0.28284271
[7,] 0.03563483 -0.03984095 -0.03984095  0.0942809 -0.25819889  0.28284271  0.20000000
```

`Drawing_Dsd`
=============


```
> V<- Ppi(c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8))
> K<- V%*%t(V)
> echant1<-Drawing_Dsd(V,s=7,B=TRUE)
> print(echant1)
     Sample 1 Sample 2 Sample 3 Sample 4 Sample 5 Sample 6 Sample 7
[1,]        0        0        0        0        0        0        0
[2,]        1        1        1        1        1        1        1
[3,]        1        1        1        1        1        1        1
[4,]        0        1        0        0        1        0        0
[5,]        0        0        1        1        0        1        0
[6,]        1        1        1        1        1        0        1
[7,]        1        0        0        0        0        1        1
> echant2<-Drawing_Dsd(V,s=7,B=FALSE)
> print(echant2)
     Sample 1 Sample 2 Sample 3 Sample 4 Sample 5 Sample 6 Sample 7
[1,]        2        1        1        1        2        1        1
[2,]        3        2        2        2        3        2        2
[3,]        6        4        6        5        5        6        5
[4,]        7        6        7        6        7        7        7
> omega=matrix(0,4,7)
> rho=matrix(0,4,6)
> K<-CaDsd(omega,rho,M=4,pi=c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8))
> echant3<-Drawing_Dsd(K$EigenBasis,s=4,B=FALSE)
> print(echant3)
     Sample 1 Sample 2 Sample 3 Sample 4
[1,]        2        1        1        1
[2,]        3        2        2        2
[3,]        5        4        3        4
[4,]        6        7        5        5
```

Notation
------------
We introduce the following notation. We note U the set on which the sampling designs are defined, N the size of the set, n the size of the sampling and K the kernel matrices. For any square matrix K indexed
by U and s &sube; U, K<SUB>\|s</SUB> denotes the submatrix of K whose rows and columns are indexed by s. 

Determinantal Sampling Designs
-------------
Determinantal Sampling Designs are a family of sampling designs that adresses all theses issues. They are indexed by Hermitian contracting matrices, called kernel, and have known inclusion probabilities for any order, among many others properties. 

__Definition:__  A sampling design &#8473; on a finite set U is a __determinantal sampling design__ if there exists a Hermitian matrix K indexed by U, called kernel, such that for all s &isin;  2<SUP>U</SUP>, &#8721; <SUB>s\' &supe;  s</SUB> &#8473; ( s\' ) = det( K<SUB>\|s</SUB> ). This sampling design is denoted by __DSD(K)__.

This Package
-------------
Because they have such good properties but their construction is not so easy, it is important to develop tools to make it available. This package provides two algorithms of constructions of DSD's kernel matrices and a sampling algorithm. All the results used in this package are presented in the scientific paper [_Determinantal Sampling Designs_](https://www.sciencedirect.com/science/article/abs/pii/S0378375818300533) written by Vincent Loonis and Xavier Mary.


<br>  </br>


Functions
=============
loonismary
-------------
This function focuses on __fixed size *sampling designs*__ with __prescribed first order inclusion probabilities__ as it is common to work with them in practice.

It has been proven that constructing a determinantal sampling design with prescribed first order inclusion probabilities is equivalent to constructing a projection matrix with a prescribed diagonal. 
<br>  </br>

_loonismary_ constructs such a matrix with as sole input a vector of N probabilities, which are the first order inclusion probabilities of the sampling design. The exact knowledge of the coefficients K<SUB>kl</SUB> enables a precise characterization of the sampling designs so constructed: the DSD is then totally defined by the kernel matrix in output.

_loonismary_ takes in input any prescribed vector of inclusion probabilities __&Pi;__ of size N such that &#8721; <SUB>k &isin;  U</SUB> &pi;<SUB>k</SUB> is an integer. This integer n will be the number of elements sampled of the DSD.

_loonismary_ offers two types of output, that both provide whole necessary information. The user indicates his choice between the two with the boolean __K__ in input.
<ul>
<li> First type of output is the __kernel matrix K__ described in the definition above. This option corresponds to K = TRUE.</li>
<li> Second type of output is also __a matrix, of size N*x*n__, that will be noted __V__. V is an n orthonormal basis of the range of K such that VV&#772;<SUP>T</SUP> = K. It is easily possible to retrieve the kernel matrix with this output. This option corresponds to K = FALSE.
*This output is suited to the input needed by the sampling function* lavancier *presented below.*
</ul>
<br>  </br>
In the particular case of equal probability DSDs of size n (&Pi;<SUB>k</SUB> = nN<SUP>-1</SUP> for all k in U) and when n divides N, the matrix K is a block diagonal matrix with n blocks, whose entries are nN<SUP>-1</SUP>. The resulting DSD is thus the 1-per-stratum sampling design. This design is known to be more efficient than systematic sampling of the population in natural order.

In the general case, a possible drawback of this general construction is that some of the joint probabilities equal 0 leading to difficulties in estimating the variance. Also, in the case of equal first order inclusion probability nN<SUP>-1</SUP>, the algorithm always provides the same matrix (whatever the reordering), which properties (of the associated sampling design) may be difficult to interpret unless n divides N, as explained before. The next function, _periodicdsd_ provides an other construction that circumvent these drawbacks, but only in the case of equal first order inclusion probability.


```
> loonismary(c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8))
          [,1]       [,2]      [,3]       [,4]
[1,] 0.7071068  0.0000000 0.0000000  0.0000000
[2,] 0.5000000 -0.7071068 0.0000000  0.0000000
[3,] 0.5000000  0.7071068 0.0000000  0.0000000
[4,] 0.0000000  0.0000000 0.4472136  0.0000000
[5,] 0.0000000  0.0000000 0.6324555  0.0000000
[6,] 0.0000000  0.0000000 0.5163978 -0.5773503
[7,] 0.0000000  0.0000000 0.3651484  0.8164966


> loonismary(c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8), TRUE)
          [,1]       [,2]       [,3]      [,4]      [,5]       [,6]       [,7]
[1,] 0.5000000  0.3535534  0.3535534 0.0000000 0.0000000  0.0000000  0.0000000
[2,] 0.3535534  0.7500000 -0.2500000 0.0000000 0.0000000  0.0000000  0.0000000
[3,] 0.3535534 -0.2500000  0.7500000 0.0000000 0.0000000  0.0000000  0.0000000
[4,] 0.0000000  0.0000000  0.0000000 0.2000000 0.2828427  0.2309401  0.1632993
[5,] 0.0000000  0.0000000  0.0000000 0.2828427 0.4000000  0.3265986  0.2309401
[6,] 0.0000000  0.0000000  0.0000000 0.2309401 0.3265986  0.6000000 -0.2828427
[7,] 0.0000000  0.0000000  0.0000000 0.1632993 0.2309401 -0.2828427  0.8000000


> loonismary(rep(120/4000, 4000))
```

periodicdsd
-------------
This function constructs __fixed size, equal probabilities *sampling designs*__, that exhibit some __periodic behavior__. The kernels involved are special [Toeplitz matrices](https://en.wikipedia.org/wiki/Toeplitz_matrix) constructed upon primitive N<SUP>th</SUP> roots of the unity.

_periodicdsd_ constructs a projection Toeplitz matrix with as inputs:
<ul>
<li> an integer __n__, the size of the sampling.</li>
<li> an integer __r__ that influences the periodicity of the DSD. The bigger r is, the bigger the frequency is. </li>
<li> an integer __N__, such that n &le; N, r < N and r and N are two relatively prime integers. N is the size of the set U.
</ul>

As the sampling is of fixed size n, the first order inclusion probabilities are nN<SUP>-1</SUP>.

The exact knowledge of the coefficients K<SUB>kl</SUB> enables a precise characterization of the sampling designs so constructed: the DSD is then totally defined by the kernel matrix in output.

_periodicdsd_ offers two types of output, that both provide whole necessary information. The user indicates his choice between the two with the boolean __K__ in input.
<ul>
<li> First type of output is the __kernel matrix K__ with the properties described above. This option corresponds to K = TRUE.</li>
<li> Second type of output is also __a matrix, of size N*x*n__, that will be noted V. V is an n orthonormal basis of the range of K such that VV&#772;<SUP>T</SUP> = K. It is easily possible to retrieve the kernel matrix with this output. This option corresponds to K = FALSE.
*This output is suited to the input needed by the sampling function* lavancier *presented below.*
</ul>

```
> periodicdsd(4, 1, 7)
             [,1]                  [,2]                  [,3]                  [,4]
[1,] 0.3779645+0i  0.2356570+0.2955045i -0.0841050+0.3684881i -0.3405342+0.1639926i
[2,] 0.3779645+0i -0.0841050+0.3684881i -0.3405342-0.1639926i  0.2356570-0.2955045i
[3,] 0.3779645+0i -0.3405342+0.1639926i  0.2356570-0.2955045i -0.0841050+0.3684881i
[4,] 0.3779645+0i -0.3405342-0.1639926i  0.2356570+0.2955045i -0.0841050-0.3684881i
[5,] 0.3779645+0i -0.0841050-0.3684881i -0.3405342+0.1639926i  0.2356570+0.2955045i
[6,] 0.3779645+0i  0.2356570-0.2955045i -0.0841050-0.3684881i -0.3405342-0.1639926i
[7,] 0.3779645+0i  0.3779645-0.0000000i  0.3779645-0.0000000i  0.3779645-0.0000000i

> periodicdsd(4, 1, 7, K = TRUE)
                       [,1]                   [,2]                   [,3]                   [,4]
[1,] 0.57142857+0.00000000i 0.07142857-0.31294902i 0.07142857+0.03439819i 0.07142857-0.08956860i
[2,] 0.07142857+0.31294902i 0.57142857+0.00000000i 0.07142857-0.31294902i 0.07142857+0.03439819i
[3,] 0.07142857-0.03439819i 0.07142857+0.31294902i 0.57142857+0.00000000i 0.07142857-0.31294902i
[4,] 0.07142857+0.08956860i 0.07142857-0.03439819i 0.07142857+0.31294902i 0.57142857+0.00000000i
[5,] 0.07142857-0.08956860i 0.07142857+0.08956860i 0.07142857-0.03439819i 0.07142857+0.31294902i
[6,] 0.07142857+0.03439819i 0.07142857-0.08956860i 0.07142857+0.08956860i 0.07142857-0.03439819i
[7,] 0.07142857-0.31294902i 0.07142857+0.03439819i 0.07142857-0.08956860i 0.07142857+0.08956860i
                       [,5]                   [,6]                   [,7]
[1,] 0.07142857+0.08956860i 0.07142857-0.03439819i 0.07142857+0.31294902i
[2,] 0.07142857-0.08956860i 0.07142857+0.08956860i 0.07142857-0.03439819i
[3,] 0.07142857+0.03439819i 0.07142857-0.08956860i 0.07142857+0.08956860i
[4,] 0.07142857-0.31294902i 0.07142857+0.03439819i 0.07142857-0.08956860i
[5,] 0.57142857+0.00000000i 0.07142857-0.31294902i 0.07142857+0.03439819i
[6,] 0.07142857+0.31294902i 0.57142857+0.00000000i 0.07142857-0.31294902i
[7,] 0.07142857-0.03439819i 0.07142857+0.31294902i 0.57142857+0.00000000i

> periodicdsd(15, 1, 365)
> periodicdsd(15, 2, 365)
```
###  Influence of the parameter r on the periodicity characteristics:

*Inclusion probabilities &Pi;<SUB>ik</SUB> (grey curve) and selected observation (black vertical line) at step i for the first four steps of sampling for constant values of n and N and different values of r.*

r = 1

![Alt text](r1s1.png){ width=20% }  ![Alt text](r1s2.png){ width=20% }
![Alt text](r1s3.png){ width=20% } ![Alt text](r1s4.png){ width=20% }

r = 2

![Alt text](r2s1.png){ width=20% } ![Alt text](r2s2.png){ width=20% }
![Alt text](r2s3.png){ width=20% } ![Alt text](r2s4.png){ width=20% }

r = 3


![Alt text](r3s1.png){ width=20% } ![Alt text](r3s2.png){ width=20% } ![Alt text](r3s3.png){ width=20% } ![Alt text](r3s4.png){ width=20% }

r = 4


![Alt text](r4s1.png){ width=20% } ![Alt text](r4s2.png){ width=20% } ![Alt text](r4s3.png){ width=20% } ![Alt text](r4s4.png){ width=20% }


lavancier
-------------
Finally, the package provides a sampling function suited to the outputs of the previous functions and to the construction method of determinantal sampling designs presented in this package. 

*lavancier* samples from fixed size determinantal sampling designs.

This function takes 5 inputs:
<ul>
<li> __V__, a matrix,</li>
<li> __Pi__, a vector of first order inclusion probabilities weights,</li> 
<li> __s__, an integer, </li>
<li> __B__, a boolean, </li>
<li> __C__, a vector of boolean of size n. </li>
</ul>

The function let the choice of two models of input to the user:
<ul>
<li> It can take in input a matrix __V__ of an n orthonormal basis of the range of K, the kernel of the DSD. This matrix V is given by the previous two functions for fixed size sampling designs with prescribed first order inclusion probabilities (*loonismary*) or for fixed size, equal probabilities sampling designs (*periodicdsd*).</li>
<li> Or it can take directly a vector of first order inclusion probabilities weights __Pi__ if the user doesn't have any matrix describing his sampling design. Then, the function calls an other function of the package (*loonismary*) to build the matrix __V__.</li>
</ul>

The function proposes to take one or more samples from the same DSD and returns, depending of the number of samples, 
<ul>
<li> a vector with the sample, if only one sample was asked ( __s__ = 1).</li>
<li> a dataframe with one column for one sample, if several samples were asked ( __s__ >1).</li>
</ul>

There is also a possibility to plot the evolution of the probability weights at every step of the draw. Therefor, the user has to indicate for each step of the draw if the probability weights have to be plotted or not with a vector of booleans __C__ of size n.

_Example_: For a DSD of size 4 in a set of size 7, C is a vector of 4 booleans. 
C = (TRUE, FALSE, FALSE, FALSE) will only plot the input probabilities.
C = (FALSE, FALSE, TRUE, TRUE) will plot the probabilities for the draws of the second and the fourth elements of the sample.

*By default, C is a vector of n boolean = FALSE*.

The last input is a boolean __B__, that specifies the choice of the user about the expression of the output. 
<ul>
<li> *__B__ = TRUE* returns __a vector of size N, with only ones and zeros__. The ones indicate the elements sampled.</li>
<li> *__B__ = FALSE* returns __a vector of size n, containing the elements sampled__.</li>
</ul>

```
> V <- loonismary(c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8))

#lavancier used with a V matrix
> lavancier(V)
[1] 1 2 6 7

#lavancier used with a vector of first order inclusion probabilities weights
> lavancier(Pi = rep(12/36, 36))
[1]  1  5  8 12 15 17 20 24 26 29 32 36

#first type of output: vector of zeros and ones
> lavancier(V, B = TRUE)
[1] 1 0 1 0 1 1 0
#here, the elements sampled are 1, 3, 5 and 6

#second type of output (default output): vector of elements sampled
> lavancier(V, B = FALSE)
[1] 1 3 6 7

#several samples drawn
> lavancier(V, 5)
     Sample 1 Sample 2 Sample 3 Sample 4 Sample 5
[1,]        1        2        1        1        2
[2,]        3        3        3        2        3
[3,]        6        5        6        4        5
[4,]        7        6        7        7        7

#using the option of representation of the evolution of probability weights at 2nd, 3rd and 4th draws.
> lavancier(V, B = TRUE, C = c(FALSE, TRUE, TRUE, TRUE))
[1] 1 1 0 0 0 1 1
```
![Alt text](tirage2.jpg){ width=30% }
![Alt text](tirage3.jpg){ width=30% }
![Alt text](tirage4.jpg){ width=30% }
```
> V <- periodicdsd(15, 1, 365)

#lavancier can also be used with V matrices produced by periodicdsd
> lavancier(V, C = c(rep(TRUE, 5), rep(FALSE, 10)))
[1]   3  26  43  68  89 101 147 163 192 223 258 269 297 323 341
```
![Alt text](tiragec1.jpg){ width=30% }
![Alt text](tiragec2.jpg){ width=30% }
![Alt text](tiragec3.jpg){ width=30% }
![Alt text](tiragec4.jpg){ width=30% }
![Alt text](tiragec5.jpg){ width=30% }
