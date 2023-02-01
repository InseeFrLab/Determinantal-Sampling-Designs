Introduction
=============

In the sequel, we provide some R-tools to help implement Determinantal Sampling Designs such as described in [Loonis and Mary (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0378375818300533) or in [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux). The tools consist of the translation in $R$ of tools that were written at first in _SAS_. The translation is due to Kim Antunez and Fabrice Nathan Tchazou Kamwa from _Ensae_, Jeanne Monnier from _Ecole Centrale de Lyon_, and Loik Acakpo Addra from _Gustave Eiffel University_. So far, all the programs are experimental. It means that they were strongly tested, but it may remain some bugs. Moreover, they were used first and foremost to help carry out the simulation studies described in [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux). Should you have any questions or remarks, please feel free to contact Vincent Loonis at vincent.loonis (a)insee.fr . 

1. `Ppi` provides a real matrix ${\Phi^N}^\intercal$, whose size is $(N\times n)$, and such that ${\Phi^N}^\intercal\Phi^N=P^\Pi$, where $P^\Pi$ is the projection matrix associated to $\Pi$, such as in described [Loonis and Mary (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0378375818300533).
2. `CaDsd` stands for _Contructing All Determintal Designs_. For a given vector $\Pi$,  `CaDsd` will provide __almost__ all the hermitian contracting matrices whose diagonal is $\Pi$. `CaDsd` implements [Fickus and Al. (2013)](https://link.springer.com/chapter/10.1007/978-0-8176-8373-3_2) algorithm, according to the reparametrisation due to [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux).
3. `Drawing_Dsd` will draw one or several samples from a random variable $\mathbb{S}$ whose law is a determinantal sampling design and whose kernel is a __projection__ matrix $K$. `Drawing_Dsd` relies on an algorithm introduced by [ Lavancier et collab., (2015)](https://www.jstor.org/stable/24775312#metadata_info_tab_contents), and presented in both [Loonis and Mary (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0378375818300533) and [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux).
4. `Reciprocal_CaDsd` will provide the set of parameters that led to an hermitian contracting matrix $K^{\Pi^\triangleright}$ according to [Fickus and Al. (2013)](https://link.springer.com/chapter/10.1007/978-0-8176-8373-3_2) algorithm, using the implementation described in [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux).

`Ppi`
=============

Let $\Pi$ be a vector in $\]0,1\[^N$, such that $\overset{N}{\underset{k=1}\sum}\Pi_k=n\in \mathbb{N}$. 
`Ppi`($\Pi$) will construct an orthomormal eigenbasis for the projection matrix $P^\Pi$ that was introduced by [Loonis and Mary (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0378375818300533). The construction relies on the _fast_ algorithm available in [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux) (Algorithm 7.1).  We give below an example for $\Pi=(0.5,0.75,0.75,0.2,0.4,0.6,0.8)^\intercal$. The first matrix is the eigen basis, the second is $K=P^\Pi$. 

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

Let $\Pi$ a vector in $\]0,1\[^N$, such that $\underset{k=1}{\overset{N}\sum}\Pi_k=\mu \in \mathbb{R}^+$. Let $M\geq \mu$ be an integer and $\Omega$ and $\rho$ two matrices whose size is respectively $(M\times N)$ and $(M\times (N-1))$ and whose entries lie in $\[0,1\]$. `CaDsd`$(\Pi,M,\Omega,\rho)$ provides a list that consists of : 

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

We show how to obtain the matrices $Q^{\Pi^\triangleright}$ described in [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux)

```

omega=matrix(1,4,30)
rho=matrix(0,4,29)
# probas constantes et n=4 ne divise pas N=30
K_divise_pas<-CaDsd(omega,rho,M=4,pi=matrix(4/30,1,30))
# probas constantes et n=5 ne divise N=30
omega=matrix(1,5,30)
rho=matrix(0,5,29)
K_divise<-CaDsd(omega,rho,M=5,pi=matrix(5/30,1,30))
# proba quelconques et n=5 
pi<-runif(matrix(0,1,30))
pi<-5*pi/sum(pi)
K_general<-CaDsd(omega,rho,M=5,pi=pi)



```


`Drawing_Dsd`
=============

Let $K$ be a __projection__ matrix, and $\overline{\Phi}^\intercal$ be one of its orthonormal eigenbasis. `Drawing_DSD`$(\overline{\Phi}^\intercal,s,B)$ will provide a set of $s$ drawings from $DSD(K)$. $B$ is a boolean variable, that will shape the output as shwon below. 
`Drawing_Dsd` can be used directly after `Ppi` or `CaDsd`. 

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

`Reciprocal_CaDsd`
=============

Let $K^\Pi$ be a $(N \times N)$ hermitian contracting matrix whose diagonal is $\Pi$. Let $\Sigma^\intercal$ be a permutation matrix such that $\Pi^\triangleright=\Sigma^\intercal\Pi$, where the entries of $\Pi^\triangleright$ are those of $\Pi$ in descending order. According to [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux), relying on [Fickus and Al. (2013)](https://link.springer.com/chapter/10.1007/978-0-8176-8373-3_2), there is a set of multidimensionnal paramaters $(\Omega,\\{V^k\\}_{k=1}^{N-1})$, so that    
