Introduction
=============
In the following section, we present a set of R-tools designed to facilitate the implementation of Determinantal Sampling Designs, as detailed in [Loonis and Mary (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0378375818300533) or in [Loonis (2024)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2023002/article/00008-eng.htm). These tools represent the translation of functionality originally developed in SAS. The translation effort was undertaken by Kim Antunez and Fabrice Nathan Tchazou Kamwa from Ensae, Jeanne Monnier from Ecole Centrale de Lyon, and Loik Acakpo Addra from Gustave Eiffel University. It's important to note that all the programs provided here are experimental in nature. While they have undergone rigorous testing, there may still be some undiscovered bugs. Additionally, these tools were primarily utilized to facilitate the simulation studies outlined in  [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux). If you have any questions, feedback or improvment, please don't hesitate to reach out to Vincent Loonis at vincent.loonis (at) insee.fr.

1. `Ppi` provides a real matrix ${\Phi^N}^\intercal$, whose size is $(N\times n)$, and such that ${\Phi^N}^\intercal\Phi^N=P^\Pi$, where $P^\Pi$ is the __projection matrix__ associated to $\Pi$, such as in described [Loonis and Mary (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0378375818300533).
2. `CaDsd` stands for _Contructing All Determintal Sampling Designs_. For a given vector $\Pi$,  `CaDsd` provides __almost__ all the hermitian contracting matrices whose diagonal is $\Pi$. `CaDsd` implements [Fickus and Al. (2013)](https://link.springer.com/chapter/10.1007/978-0-8176-8373-3_2) algorithm, according to the reparametrisation due to [Loonis (2024)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2023002/article/00008-eng.htm).
3. `Periodic_Dsd` builds a complex matrix $\overline{\Phi^{p,N,n}}^\intercal$ such that $\overline{\Phi^{p,N,n}}^\intercal\Phi^{p,N,n}=K^{p,N,n}$, where $p$ and $N$ are prime and $K^{p,N,n}$ is a __projection__ matrix described in [Loonis and Mary (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0378375818300533). The determinantal sampling design $DSD(K^{p,N,n})$ has some interesting periodic features.
4. `Drawing_Dsd` draws one or several samples from a random variable $\mathbb{S}$ whose law is a __fixed-size__ determinantal sampling design, that is to say whose kernel is a __projection__ matrix $K$. `Drawing_Dsd` relies on an algorithm introduced by [ Lavancier et collab., (2015)](https://www.jstor.org/stable/24775312#metadata_info_tab_contents), and presented in both [Loonis and Mary (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0378375818300533) and [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux). `Drawing_Dsd` can be used, for example, after any of the 3 previous functions.
5. `Reciprocal_CaDsd` will provide the set of parameters that led to an hermitian contracting matrix $K^{\Pi^\triangleright}$ according to [Fickus and Al. (2013)](https://link.springer.com/chapter/10.1007/978-0-8176-8373-3_2) algorithm, using the implementation described in[Loonis (2024)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2023002/article/00008-eng.htm).

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

1. the hemitian matrix $K^{\Pi^\triangleright}(M,\Omega,\rho)$ such as described in  [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux), where the entries of $\Pi^\triangleright$ are those of $\Pi$ in descending order. $K^{\Pi^\triangleright}(M,\Omega,\rho)$ will be a __projection__ matrix iff $M=\mu=n \in \mathbb{N}^*$.
2. the sequence of spectra of all $K^{\Pi^\triangleright}(M,\Omega,\rho)$'s principal submatrices
3. an eigenbasis $\overline{\Phi^{p,N,n}}^\intercal$ for $K^{\Pi^\triangleright}(M,\Omega,\rho)$, that directly be used to draw a sample, when $K^{\Pi^\triangleright}(M,\Omega,\rho)$ is a __projection__ matrix.    

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

`Periodic_Dsd`
=============

Let $p,N,n$ be three integers such that $n \leq N$ and $p$ and $N$ are prime. There exists a __projection__ matrix $K^{p,N,n}$ such that $K_{kl}^{p,N,n}=\frac{1}{N}\frac{\sin(\frac{np(k-l)\pi}{N})}{\sin(\frac{p(k-l)\pi}{N})}\exp({\frac{i p(n-1)(k-l)\pi}{N}})$ and $K_{kk}^{p,N,n}=\frac{n}{N}$.
`Periodic_Dsd` provides both $K^{p,N,n}$ and a set of eigenvectors $\overline{\Phi^{p,N,n}}^\intercal$ that can be directly used in `Drawing_Dsd`. The $k^{th}$ column of $\Phi^{p,N,n}$ is given by $\varphi_k^{p,N,n}=((z^{k})^0,\cdots,(z^k)^{j},\cdots,(z^{k})^{n-1})^\intercal/\sqrt{N}$, where $z=\exp({\frac{2i\pi p}{N}})$.

```
> test<-periodic_dsd(n=3, p=3, N=7) 
> print(test$KpNn)
                       [,1]                   [,2]                   [,3]                   [,4]
[1,]  0.4285714+0.00000000i  0.1032173+0.04970682i  0.2001384+0.25096563i -0.0176414-0.07729202i
[2,]  0.1032173-0.04970682i  0.4285714+0.00000000i  0.1032173+0.04970682i  0.2001384+0.25096563i
[3,]  0.2001384-0.25096563i  0.1032173-0.04970682i  0.4285714+0.00000000i  0.1032173+0.04970682i
[4,] -0.0176414+0.07729202i  0.2001384-0.25096563i  0.1032173-0.04970682i  0.4285714+0.00000000i
[5,] -0.0176414-0.07729202i -0.0176414+0.07729202i  0.2001384-0.25096563i  0.1032173-0.04970682i
[6,]  0.2001384+0.25096563i -0.0176414-0.07729202i -0.0176414+0.07729202i  0.2001384-0.25096563i
[7,]  0.1032173+0.04970682i  0.2001384+0.25096563i -0.0176414-0.07729202i -0.0176414+0.07729202i
                       [,5]                   [,6]                   [,7]
[1,] -0.0176414+0.07729202i  0.2001384-0.25096563i  0.1032173-0.04970682i
[2,] -0.0176414-0.07729202i -0.0176414+0.07729202i  0.2001384-0.25096563i
[3,]  0.2001384+0.25096563i -0.0176414-0.07729202i -0.0176414+0.07729202i
[4,]  0.1032173+0.04970682i  0.2001384+0.25096563i -0.0176414-0.07729202i
[5,]  0.4285714+0.00000000i  0.1032173+0.04970682i  0.2001384+0.25096563i
[6,]  0.1032173-0.04970682i  0.4285714+0.00000000i  0.1032173+0.04970682i
[7,]  0.2001384-0.25096563i  0.1032173-0.04970682i  0.4285714+0.00000000i
> print(test$Phi)
             [,1]                  [,2]                  [,3]
[1,] 0.3779645+0i -0.3405342+0.1639926i  0.2356570-0.2955045i
[2,] 0.3779645+0i  0.2356570-0.2955045i -0.0841050-0.3684881i
[3,] 0.3779645+0i -0.0841050+0.3684881i -0.3405342-0.1639926i
[4,] 0.3779645+0i -0.0841050-0.3684881i -0.3405342+0.1639926i
[5,] 0.3779645+0i  0.2356570+0.2955045i -0.0841050+0.3684881i
[6,] 0.3779645+0i -0.3405342-0.1639926i  0.2356570+0.2955045i
[7,] 0.3779645+0i  0.3779645-0.0000000i  0.3779645-0.0000000i
```

`Drawing_Dsd`
=============

Let $K$ be a __projection__ matrix, and $\overline{\Phi}^\intercal$ be one of its orthonormal eigenbasis. `Drawing_DSD`$(\overline{\Phi}^\intercal,s,B)$ will provide a set of $s$ drawings from $DSD(K)$. $B$ is a boolean variable, that will shape the output as shwon below. 
`Drawing_Dsd` can be used directly after `Ppi`, `CaDsd` or `Periodic_Dsd`. 

```
> Phi<- Ppi(c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8))
> K<- Phi%*%t(Phi)
> echant1<-Drawing_Dsd(Phi,s=7,B=TRUE)
> print(echant1)
     Sample 1 Sample 2 Sample 3 Sample 4 Sample 5 Sample 6 Sample 7
[1,]        0        0        0        0        0        0        0
[2,]        1        1        1        1        1        1        1
[3,]        1        1        1        1        1        1        1
[4,]        0        1        0        0        1        0        0
[5,]        0        0        1        1        0        1        0
[6,]        1        1        1        1        1        0        1
[7,]        1        0        0        0        0        1        1
> echant2<-Drawing_Dsd(Phi,s=7,B=FALSE)
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
=================

Let $K^\Pi$ be a $(N \times N)$ hermitian contracting matrix whose diagonal is $\Pi$. Let $\Sigma^\intercal$ be a permutation matrix such that $\Pi^\triangleright=\Sigma^\intercal\Pi$, where the entries of $\Pi^\triangleright$ are those of $\Pi$ in descending order. According to [Loonis (2023)](https://www.researchgate.net/publication/359095103_Construire_tous_les_plans_de_sondage_determinantaux),there is a set of multidimensionnal paramaters $(\Omega,\\{V^k\\}_{k=1}^{N-1})$, that leads to $K^{\Pi^\triangleright}=\Sigma^\intercal K^\Pi \Sigma$ when using [Fickus and Al. (2013)](https://link.springer.com/chapter/10.1007/978-0-8176-8373-3_2) algorithm (see `CaDsd`). `Reciprocal_CaDsd` provides a list that consists of these parameters, along with $K^{\Pi^\triangleright}$ and $\rho$, where $rho$ is a real matrix, whose $k^{th}$ column is $\rho_k$ and such that $V^k$'s diagonal is given by $\exp(2i\pi\rho_k)$. 

```
M<-5
N<-10
n<-3
omega=matrix(runif(N*M),nrow=M)
rho=matrix(runif((N-1)*M),nrow=M)
 
 m_pi=c()
for (k in 1:10){
   m_pi[k]=2*n*k/(N*(N+1))
 }
 
 
 fic1<-CaDsd(omega,rho,M,pi=m_pi)
 Rfic1=Reciprocal_CaDsd(fic1$K)
 fic2<-CaDsd(Rfic1$omega,Rfic1$rho,M,pi=m_pi)
 
 print(omega)
          [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]       [,9]      [,10]
[1,] 0.3923160 0.2156829 0.2416166 0.2554178 0.4725795 0.3066209 0.9301133 0.1675726 0.67872843 0.90085664
[2,] 0.1519017 0.9063653 0.4604553 0.2622157 0.4536076 0.1326980 0.6871964 0.7136560 0.67827077 0.12909050
[3,] 0.5060835 0.8751697 0.9850278 0.7559215 0.6187970 0.8321663 0.2215734 0.3717336 0.72666554 0.04357961
[4,] 0.1107567 0.1113637 0.8075908 0.5961480 0.5437580 0.9520614 0.2817436 0.7840146 0.09198391 0.42699150
[5,] 0.9853772 0.5131908 0.4567270 0.6186477 0.3110066 0.5680391 0.4236700 0.7168322 0.86583961 0.64028894
> print(Rfic1$omega)
     [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]      [,9]     [,10]
[1,]    0 0.0000000 0.0000000 0.0000000 0.4725797 0.3066194 0.9301118 0.1675742 0.6787289 0.9008571
[2,]    0 0.0000000 0.0000000 0.2622155 0.4536051 0.1326989 0.6871965 0.7136504 0.6782756 0.1290919
[3,]    0 0.0000000 0.9850274 0.7559161 0.6187984 0.8321581 0.2215780 0.3717475 0.7266663 0.0435790
[4,]    0 0.1109931 0.8075823 0.5961483 0.5437617 0.9520535 0.2817430 0.7840032 0.0919844 0.4269923
[5,]    0 0.0000000 1.0000000 0.0000000 0.0000000 0.0000000 0.0000000 1.0000000 0.0000000 0.0000000
 
print(rho)
          [,1]       [,2]      [,3]      [,4]      [,5]       [,6]       [,7]      [,8]       [,9]
[1,] 0.2357493 0.15930012 0.3643466 0.7077350 0.3639960 0.05490594 0.19111256 0.6926908 0.36236432
[2,] 0.6635870 0.37681199 0.8386177 0.1067836 0.6525946 0.52289769 0.56952704 0.4784134 0.20116611
[3,] 0.7687070 0.02339084 0.9799671 0.4075552 0.1280213 0.44504825 0.33892575 0.8805832 0.09092541
[4,] 0.6115447 0.76314229 0.7932975 0.6978815 0.4963314 0.74728106 0.07423949 0.8453083 0.98745994
[5,] 0.8215298 0.06298310 0.3683020 0.6088772 0.8523733 0.23277038 0.92333140 0.4824277 0.02205846

print(Rfic1$rho)
          [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]      [,9]
[1,] 0.2357498 0.1592902 0.3643563 0.7077356 0.3639949 0.0549060 0.1911124 0.6926893 0.3623693
[2,] 0.9782203 0.3768119 0.8386194 0.1067830 0.6525948 0.5228979 0.5695266 0.4784130 0.2011667
[3,] 0.5671799 0.6762216 0.9799672 0.4075540 0.1280211 0.4450481 0.3389240 0.8805802 0.0909189
[4,] 0.9159252 0.0671987 0.4355104 0.6978815 0.4963315 0.7472810 0.0742399 0.8453092 0.9874636
[5,] 0.5702553 0.8995156 0.7945955 0.1848984 0.8523734 0.2327704 0.9233324 0.4824268 0.0220587
 
print(fic1$K)
                      phi                                                                                                                                    
phi  0.5454545+0.0000000i  0.0002706+0.0030142i  0.0003731+0.0020612i -0.0001052-0.0305094i -0.0500329+0.0294792i  0.0336857-0.0325909i  0.0403131+0.0201735i
     0.0002706-0.0030142i  0.4909091+0.0000000i -0.0244842+0.0240354i  0.0460412-0.0471265i  0.0702979+0.0555195i -0.0480811+0.0843178i  0.0032255-0.0525607i
     0.0003731-0.0020612i -0.0244842-0.0240354i  0.4363635+0.0000000i  0.1520803-0.0230265i  0.0314145+0.0217865i -0.0547662-0.0915697i -0.0831685-0.0044739i
    -0.0001052+0.0305094i  0.0460412+0.0471265i  0.1520803+0.0230265i  0.3818181+0.0000000i -0.0354837-0.1216975i  0.1470988+0.0799118i  0.0648282+0.0528967i
    -0.0500329-0.0294792i  0.0702979-0.0555195i  0.0314145-0.0217865i -0.0354837+0.1216975i  0.3272728+0.0000000i -0.0363378-0.0504562i  0.0674037+0.1681836i
     0.0336857+0.0325909i -0.0480811-0.0843178i -0.0547662+0.0915697i  0.1470988-0.0799118i -0.0363378+0.0504562i  0.2727273+0.0000000i  0.0226656+0.0979715i
     0.0403131-0.0201735i  0.0032255+0.0525607i -0.0831685+0.0044739i  0.0648282-0.0528967i  0.0674037-0.1681836i  0.0226656-0.0979715i  0.2181819+0.0000000i
     0.0560351-0.1420535i  0.0218787+0.0228373i -0.0818317+0.0467617i -0.0077228-0.1237909i -0.0228359-0.0404993i  0.0438964-0.1170160i  0.1337501-0.0599889i
     0.0883756-0.0457672i -0.0457612-0.0462458i  0.0249098-0.0774406i  0.0413521+0.0452907i -0.0545141-0.0875429i -0.0486826+0.0939439i -0.0191516-0.0256328i
     0.0805378-0.0633756i  0.0005886+0.0236914i -0.0535510-0.0254798i -0.0359029+0.0314095i -0.0574232+0.0137265i -0.0352750+0.0606964i -0.0315574-0.0282031i
                                                                     
phi  0.0560351+0.1420535i  0.0883756+0.0457672i  0.0805378+0.0633756i
     0.0218787-0.0228373i -0.0457612+0.0462458i  0.0005886-0.0236914i
    -0.0818317-0.0467617i  0.0249098+0.0774406i -0.0535510+0.0254798i
    -0.0077228+0.1237909i  0.0413521-0.0452907i -0.0359029-0.0314095i
    -0.0228359+0.0404993i -0.0545141+0.0875429i -0.0574232-0.0137265i
     0.0438964+0.1170160i -0.0486826-0.0939439i -0.0352750-0.0606964i
     0.1337501+0.0599889i -0.0191516+0.0256328i -0.0315574+0.0282031i
     0.1636356+0.0000000i -0.0318199-0.0230229i -0.0028158-0.0019392i
    -0.0318199+0.0230229i  0.1090908+0.0000000i  0.0498577+0.0218657i
    -0.0028158+0.0019392i  0.0498577-0.0218657i  0.0545449+0.0000000i

print(fic2$K)
                      phi                                                                                                                                    
phi  0.5454545+0.0000000i  0.0002706+0.0030142i  0.0003731+0.0020612i -0.0001053-0.0305089i -0.0500327+0.0294792i  0.0336858-0.0325909i  0.0403126+0.0201735i
     0.0002706-0.0030142i  0.4909091+0.0000000i -0.0244844+0.0240356i  0.0460414-0.0471268i  0.0702982+0.0555197i -0.0480817+0.0843178i  0.0032252-0.0525606i
     0.0003731-0.0020612i -0.0244844-0.0240356i  0.4363637+0.0000000i  0.1520799-0.0230263i  0.0314145+0.0217865i -0.0547660-0.0915700i -0.0831687-0.0044733i
    -0.0001053+0.0305089i  0.0460414+0.0471268i  0.1520799+0.0230263i  0.3818182+0.0000000i -0.0354836-0.1216974i  0.1470985+0.0799119i  0.0648281+0.0528967i
    -0.0500327-0.0294792i  0.0702982-0.0555197i  0.0314145-0.0217865i -0.0354836+0.1216974i  0.3272726+0.0000000i -0.0363377-0.0504563i  0.0674037+0.1681838i
     0.0336858+0.0325909i -0.0480817-0.0843178i -0.0547660+0.0915700i  0.1470985-0.0799119i -0.0363377+0.0504563i  0.2727273+0.0000000i  0.0226653+0.0979713i
     0.0403126-0.0201735i  0.0032252+0.0525606i -0.0831687+0.0044733i  0.0648281-0.0528967i  0.0674037-0.1681838i  0.0226653-0.0979713i  0.2181819+0.0000000i
     0.0560358-0.1420542i  0.0218784+0.0228375i -0.0818316+0.0467622i -0.0077226-0.1237912i -0.0228357-0.0404987i  0.0438971-0.1170165i  0.1337501-0.0599886i
     0.0883755-0.0457675i -0.0457614-0.0462458i  0.0249095-0.0774406i  0.0413518+0.0452906i -0.0545147-0.0875425i -0.0486826+0.0939440i -0.0191518-0.0256330i
     0.0805373-0.0633765i  0.0005894+0.0236918i -0.0535509-0.0254794i -0.0359024+0.0314093i -0.0574241+0.0137267i -0.0352752+0.0606967i -0.0315576-0.0282041i
                                                                     
phi  0.0560358+0.1420542i  0.0883755+0.0457675i  0.0805373+0.0633765i
     0.0218784-0.0228375i -0.0457614+0.0462458i  0.0005894-0.0236918i
    -0.0818316-0.0467622i  0.0249095+0.0774406i -0.0535509+0.0254794i
    -0.0077226+0.1237912i  0.0413518-0.0452906i -0.0359024-0.0314093i
    -0.0228357+0.0404987i -0.0545147+0.0875425i -0.0574241-0.0137267i
     0.0438971+0.1170165i -0.0486826-0.0939440i -0.0352752-0.0606967i
     0.1337501+0.0599886i -0.0191518+0.0256330i -0.0315576+0.0282041i
     0.1636364+0.0000000i -0.0318201-0.0230229i -0.0028156-0.0019389i
    -0.0318201+0.0230229i  0.1090908+0.0000000i  0.0498580+0.0218660i
    -0.0028156+0.0019389i  0.0498580-0.0218660i  0.0545453+0.0000000i
 
print(Rfic1$MatV$V5)
                            [,1]                        [,2]                        [,3]                        [,4]                        [,5]
[1,] -6.565439e-01+7.542522e-01i  2.118896e-05+6.484054e-06i  2.020036e-05+8.102525e-06i -1.103583e-05+1.326503e-06i  2.322065e-07-7.836604e-07i
[2,]  2.273816e-05-5.381057e-06i -5.745070e-01-8.184763e-01i  1.095546e-05+1.741139e-05i -8.069126e-06-4.376801e-06i  5.164651e-07-2.574969e-07i
[3,] -1.977464e-05+9.112115e-06i -1.467428e-05-1.362413e-05i  6.935471e-01+7.203904e-01i  8.590383e-06+2.763065e-06i -6.169201e-07+1.000569e-06i
[4,] -7.771306e-06+7.295390e-06i -7.963292e-06-3.405068e-06i -7.913839e-06-3.951831e-06i -9.997300e-01+2.304726e-02i -7.609696e-08+1.590219e-07i
[5,] -6.560823e-07+3.208232e-07i -4.429434e-07-1.805546e-07i -7.009369e-07-8.816060e-07i  2.124768e-07+2.959982e-08i  5.997837e-01-8.001621e-01i
```
