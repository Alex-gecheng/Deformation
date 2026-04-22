# Computing Locally Injective Mappings by Advanced MIPS

Xiao-Ming Fu∗ † 

∗University of Science and Technology of China 

Yang Liu† 

Baining Guo† ∗ 

† Microsoft Research 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-22/6e241e10-428b-4950-9823-6b5eb3548fa1/019bc8d07442a01271fc54c6c4acb54d2e6a78e493c2da794a60556276494eb1.jpg)



(a) AMIPS


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-22/6e241e10-428b-4950-9823-6b5eb3548fa1/2ea264a71955c3f0b475a427d9e6691a2a38c8cfb69967e8d10b2dcfc6c6a872.jpg)



(b) BDM [Lipman 2012]


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-22/6e241e10-428b-4950-9823-6b5eb3548fa1/03c0465ee39afb83f25b46baf981263a43d3b68b88fab51dbea3c277f1d0efc8.jpg)



(c) LIM [Schuller et al. 2013 ¨ ]


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-22/6e241e10-428b-4950-9823-6b5eb3548fa1/0677915d1bd855fecaed9bd98af9a84f31464d0f759b38297f4a0a8fe63c1d4c.jpg)



(d) STRICT [Levi and Zorin 2014]



Figure 1: 2D mesh deformation: comparison with the state-of-the-art methods. The input mesh $\tau$ is triangulated from a 2D regular quadrilateral mesh, and an interior handle is moved from the right side to the left side while fixing all the boundary vertices. Left to right: the strict minimizer [Levi and Zorin 2014]. Our method achieves the lowest maximal isometric distortion deformation results from our AMIPS method, bounded distortion mapping [Lipman 2012], locally injective mapping [Schuller et al ¨ . 2013] and $\delta _ { \mathrm { m a x } } ^ { i s \bar { o } } = \operatorname* { m a x } _ { { \bf t } \in \mathcal { T } } \operatorname* { m a x } \{ \sigma _ { 1 , { \bf t } } ^ { - 1 } , \sigma _ { 2 , { \bf t } } \}$ among all the methods with the least computation time, where $\sigma _ { 1 , \mathbf { t } }$ and $\sigma _ { 2 , \mathbf { t } }$ $( \sigma _ { 1 , \mathbf { t } } \le \sigma _ { 2 , \mathbf { t } } )$ ) are singular values of the Jacobian of the mapping associated with triangle t. Our method also achieves more evenly and smoothly distributed distortion (see the standard deviation of the isometric distortion $\delta _ { \mathrm { d e v } } ^ { i s o }$ ). The color on triangles encodes the isometric distortion, with white being optimal.


# Abstract

Computing locally injective mappings with low distortion in an efficient way is a fundamental task in computer graphics. By revisiting the well-known MIPS (Most-Isometric ParameterizationS) method, we introduce an advanced MIPS method that inherits the local injectivity of MIPS, achieves as low as possible distortions compared to the state-of-the-art locally injective mapping techniques, and performs one to two orders of magnitude faster in computing a mesh-based mapping. The success of our method relies on two key components. The first one is an enhanced MIPS energy function that penalizes the maximal distortion significantly and distributes the distortion evenly over the domain for both mesh-based and meshless mappings. The second is a use of the inexact block coordinate descent method in mesh-based mapping in a way that efficiently minimizes the distortion with the capability not to be trapped early by the local minimum. We demonstrate the capability and superiority of our method in various applications including mesh parameterization, mesh-based and meshless deformation, and mesh improvement. 

CR Categories: I.3.5 [Computer Graphics]: Computational Geometry and Object Modeling—Geometric algorithms, languages, and systems 

Keywords: locally injective mapping, inexact block coordinate descent methods, parameterization, deformation, mesh improvement 

# 1 Introduction

The task of computing a valid mapping is fundamental in many computer graphics and geometry processing applications, such as mesh parameterization, shape deformation and animation, image deformation and retargeting, remeshing and mesh improvement. A good mapping $f : \Omega \overline { { \subset } } \overline { { \mathbb { R } ^ { d } } }  \mathbb { R } ^ { d }$ usually satisfies the following properties: (1) $f$ is locally injective, i.e. det $J ( f ) > 0$ where $J ( f )$ is the Jacobian of $f$ ; (2) the mapping has as low as possible distortion, for instance, the singular values of $J ( f )$ should be close to 1 if an isometric mapping is desired; (3) $f$ is smooth, which is a property especially preferred in deformation; (4) the computation of $f$ is efficient enough to provide interactive feedback. It is also desired that the method of computing $f$ be applicable to both 2D and 3D mappings and easily generalized to higher dimensions. 

Mesh parameterization and deformation are typical applications in computer graphics where finding a low-distortion mapping is the central task. Numerous techniques have been developed in the past thirty years (cf. extensive surveys [Floater and Hormann 2005; Sheffer et al. 2006; Botsch and Sorkine 2008]). It is commonly known that linear based methods hardly satisfy all the above requirements although they are very efficient. For instance, the LSCM [Levy ´ et al. 2002] method cannot guarantee local injectivity and the mean value coordinate [Floater 2003a] may introduce high distortions. Instead, nonlinear approaches usually achieve lower distortion and some methods guarantee local injectivity, like MIPS [Hormann and Greiner 2000; Degener et al. 2003] and $\mathrm { \ A B F + + }$ [Sheffer et al. 2005]. However, inefficient computation and uncontrolled maximal distortions are the two main obstacles in using nonlinear approaches. Recently Schuller et al. propose a fast method to compute locally ¨ injective mappings [2013], but high distortion cannot be prevented (for instance, the region around the handle in Fig. 1c possesses high distortions). Lipman and his colleagues propose novel ways to bound the distortion explicitly [2012; 2013; 2014]. Levi and Zorin prioritize higher distortion for minimization [2014]. These methods significantly reduce distortion but with costly computation. 

The compromise of low distortion and computational efficiency is 

a major obstacle in computing locally injective mappings. In this paper, we tackle this problem by proposing a simple and efficient nonlinear approach that works for both 2D and 3D. Our method is based on the well-known MIPS method [Hormann and Greiner 2000] and the keys to our success are (1) the modification of the MIPS energy that penalizes the maximal isometric/conformal distortion significantly for both mesh-based mapping and meshless deformation; (2) the replacement of the vertex optimization scheme of MIPS by the inexact block coordinate descent method for efficient optimization that avoids early entrapment by local minimum for mesh-based applications. We call our method AMIPS (Advanced MIPS). AMIPS naturally inherits the local injectivity of MIPS and achieves comparable or lower distortion compared to state-of-the-art locally injective mapping methods like [Lipman 2012; Schuller et al¨ . 2013; Levi and Zorin 2014] with one or two orders of magnitude speedup. 

We successfully apply our method to various applications, including planar parameterization, 2D triangle and 3D tetrahedral handlebased mesh deformation, 2D and 3D handle-based meshless deformation, anisotropic tetrahedral mesh and all-hexahedral mesh improvement. We also achieve interactive response in handle-based deformation for moderately sized meshes both in 2D and 3D. 

# 2 Related Work

Locally injective mapping is desired in mesh parameterization and shape deformation. Given a disk-topology mesh, a convex combination mapping in 2D is bijective based on Tutte’s theorem if the target domain is convex [Tutte 1963; Floater 2003b]. The bijectivity is broken for concave or self-overlapping domains and has no guarantee even on 3D convex domains [Floater and Pham-Trong 2006]. The well-known MIPS method [Hormann and Greiner 2000] and its variations [Sander et al. 2001; Degener et al. 2003] provide boundary-free methods for isometric/conformal parameterization. The MIPS-type energy also finds its applications in mesh deformation [Eigensatz and Pauly 2009] and mesh improvement [Freitag and Knupp 2002; Jiao et al. 2011]. But MIPS is easily trapped by a local minimum and the computational speed is slow, as observed in [Sheffer et al. 2006]. 

Parameterization methods like $\mathrm { \ A B F + + }$ [Sheffer et al. 2005] and circle pattern [Kharevych et al. 2006] can achieve local injectivity but are hard to extend to 3D and are not applicable to mesh deformation. Schuller et al. [ ¨ 2013] design a barrier term to prevent inverted elements in their nonlinear optimization so that the mapping is always locally injective. Their method is fast enough to generate interactive results for moderate-size meshes but lacks control over extreme distortions. Lipman et al. [2012; 2013; 2014] construct a maximal convex subspace for bounding the maximal distortion and ensuring no inverted mesh elements. When their algorithms converge, locally injective mapping with bounded distortion is guaranteed but the computational cost is high, since a quadratic programming or semidefinite programming problem needs to solve at each iteration. Furthermore, a feasible solution space may not exist if the convexified subspace is empty or the upper bound of distortion is too small. 

Poranne and Lipman also introduce the bounded distortion mapping technique to 2D meshless deformation [2014] and achieve smooth and low-distortion mappings. But their method is hard to extend to 3D meshless deformation due the complexity of the bound constraint. For a fixed boundary 2D mapping, Weber et al. [2014] triangulate both the source and target domain to find a valid mapping. Jin et al. [2014] update mesh connectivity in 2D mesh deformation to make sure that a valid map exists always. These two methods are limited to only 2D and have no control on maximal distortion. Our method starts with a valid local injective mapping as in MIPS and suppresses the maximal distortion effectively. 

Distortion measurement The standard 2D MIPS energy measures the conformality of the mapping: $\sigma _ { 1 } \sigma _ { 2 } ^ { - 1 } + \sigma _ { 2 } \sigma _ { 1 } ^ { - \tilde { \Gamma } }$ where $\sigma _ { 1 } , \sigma _ { 2 }$ are the singular values of the Jacobian of the mapping associated with a triangle. For measuring isometric distortion, Degener et al. [2003] propose to minimize $\stackrel { \smile } { ( } \sigma _ { 1 } \sigma _ { 2 } ^ { - 1 } + \sigma _ { 2 } \sigma _ { 1 } ^ { - 1 } ) ( \sigma _ { 1 } \sigma _ { 2 } +$ $\sigma _ { 1 } ^ { - 1 } \sigma _ { 2 } ^ { - 1 } ) ^ { \theta }$ σ1 σ 2 , which penalizes both conformal distortion and area distortion. Other type energies like Dirichlet energy $\sigma _ { 1 } ^ { 2 } + \sigma _ { 2 } ^ { 2 }$ , stretch energy m $\operatorname { a x } \{ \sigma _ { 1 } , \sigma _ { 2 } \}$ , and Green-Lagrange energy $( \bar { \sigma _ { 1 } ^ { 2 } } - 1 ) ^ { \bar { 2 } } + ( \sigma _ { 2 } ^ { 2 } -$ $1 ) ^ { 2 }$ also build on the measurement of singular values. As rigid as Possible (ARAP), as similar as possible (ASAP) and as killing as possible (AKAP) are three alterative ways to minimize isometric/conformal distortions [Alexa et al. 2000; Igarashi et al. 2005; Sorkine and Alexa 2007; Solomon et al. 2011], but their commonly used solvers [Liu et al. 2008; Solomon et al. 2011] cannot guarantee that the mapping is locally injective. All of the above methods sum the distortions in a least squares sense and do not penalize the maximum distortion or control the distortion distribution. The strict minimizer [Levi and Zorin 2014] targets these problems and provides a minimal $L _ { \infty }$ -norm solution for distortion control. Weber et al. [2012] compute extremal quasiconformal maps to minimize the maximal conformal distortion while satisfying boundary constraints. Our method achieves an effect similar to the strict minimizer by penalizing the modified MIPS energy (see Sec. 3.3). 

Block coordinate descent methods (BCD) also called nonlinear block Gauss-Seidel methods are popular optimization tools for solving large-scale problems, especially highly structured ones. The basic idea is to update variables block by block while other variables are fixed. The standard BCD minimizes the sub-problems until they reach convergence and then switches to the next block. It is also called exact BCD. It has been widely used in many computer graphics applications, e.g., MIPS parametrization [Hormann and Greiner 2000] and mesh improvement [Freitag Diachin et al. 2006]. The disadvantages of exact BCD are mainly due to their slow convergence and the costly computation in finding exact solutions of sub-problems. Recently, inexact BCD has been proposed to overcome these problems [Bonettini 2011; Tappenden et al. 2013; Cassioli et al. 2013]. Xu and Yin show that under certain circumstances inexact BCD can achieve lower objectives and is not trapped by local minimum early [2013]. In this paper, we demonstrate that the inexact BCD solver is very suitable for computing low-distortion locally injective mappings. 

# 3 Advanced MIPS

We review the standard MIPS algorithm in Sec. 3.1 and present inexact block coordinate descent (BCD) solvers for minimizing the MIPS-type energy function efficiently in Sec. 3.2. In Sec. 3.3 we propose advanced MIPS (AMIPS) energy functions for suppressing the maximal distortion of mapping. 

# 3.1 A revisit of MIPS

The MIPS method constructs a piecewise linear interpolation function $f : \mathbb { R } ^ { 3 } \to \mathbb { R } ^ { 2 }$ that maps a triangle mesh $\tau$ with a boundary and no holes to a 2D parameterization region. The map $f$ is constituted by piecewise constant affine maps $g _ { \mathbf { t } }$ defined on triangles $\mathbf { t } \in \mathcal T$ By introducing a local coordinate system at each $\mathbf { t }$ , each $g _ { \mathbf t }$ has a linear form $g _ { \mathbf { t } } ( \mathbf { x } ) = J _ { \mathbf { t } } \mathbf { x } + b _ { \mathbf { t } }$ where $J _ { \mathbf { t } }$ is a $2 \times 2$ affine matrix that is also the Jacobian of $f$ . With the singular values of $J _ { \mathbf { t } }$ denoted as $\sigma _ { 1 }$ and $\sigma _ { 2 }$ , the MIPS energy $E _ { m i p s , \mathbf { t } } ^ { 2 D }$ defined on triangle t is 

$$
E _ {m i p s, \mathbf {t}} ^ {2 D} = \frac {\sigma_ {1}}{\sigma_ {2}} + \frac {\sigma_ {2}}{\sigma_ {1}} = \| J _ {\mathbf {t}} \| _ {F} \| J _ {\mathbf {t}} ^ {- 1} \| _ {F} = \frac {\mathrm {t r a c e} (J _ {\mathbf {t}} ^ {T} J _ {\mathbf {t}})}{\mathrm {d e t} J _ {\mathbf {t}}}.
$$

Here $\| \cdot \| _ { F }$ denotes the Frobenius norm. The MIPS energy is invariant under rigid transformation and uniform scaling. It attains 


mum of 2 when forms the overal $g _ { \mathbf t }$ is a conforIPS energy apping. The sum of. A good property of $E _ { m i p s , \mathbf { t } } ^ { 2 D }$ $E _ { m i p s } ^ { 2 D }$ goes to infinity when det $J _ { \mathbf { t } }$ approaches to zero. 

The standard MIPS algorithm for mesh parameterization is as follows: 

• The algorithm starts with a valid bijective mapping, for instance, using a barycentric mapping to map $\tau$ to a convex 2D domain. 

• Pick a vertex Update the po $\mathbf { p } \in \mathcal { T }$ randomly and  by minimizing $E _ { m i p s , \mathbf { t } } ^ { 2 D }$ emaini. Since $\breve { E } _ { m i p s , \mathbf { t } } ^ { 2 D }$ es.is a locally convex function with respect to p [Hormann 2001, Theorem 1.12], a Newton-type method can be employed to find the unique minimum while preventing triangles from being inverted. 

• The algorithm converges when no vertex can be moved further. 

Fig. 2a shows an example of a camel model parameterized using the MIPS algorithm where each iteration optimizes all the vertices once. 

# 3.2 Inexact block coordinate descent method

The block coordinate descent (BCD) method (also called nonlinear Gauss-Seidel method) [Bonettini 2011; Xu and Yin 2013] is a popular optimization tool in solving large-scale optimization problems. We consider the optimization problem 

$$
\min  _ {\mathbf {x}} F \left(x _ {1}, x _ {2}, \dots , x _ {n}\right) \tag {1}
$$

where $\mathbf { x } = ( x _ { 1 } , x _ { 2 } , \ldots , x _ { n } ) \in \Omega = \Omega _ { 1 } \times \Omega _ { 2 } \times \ldots \times \Omega _ { m } \subseteq \mathbb { R } ^ { n }$ . Here $\Omega _ { i }$ is a convex subset of $\Omega$ . Variable set x can be partitioned into $m$ blocks $\{ B _ { i } , i = 1 , \dots , m \}$ . Suppose that $F$ is multi-convex, i.e. $F$ is strictly convex with respect to variables of $B _ { i }$ in $\Omega _ { i }$ . A standard BCD method for solving Eq. 1 is described as follows. 

1. Start with an initial guess $\{ x _ { 0 } ^ { 0 } , \ldots , x _ { n } ^ { 0 } \}$ . 

2. For every $l , l \in \{ 1 , \ldots , m \}$ , solve the convex subproblem: 

$$
\min _ {B _ {l} \in \Omega_ {i}} F \left(B _ {0} ^ {k}, \dots , B _ {l - 1} ^ {k}, B _ {l}, B _ {l + 1} ^ {k - 1}, \dots , B _ {m} ^ {k - 1}\right)
$$

and update $B _ { l } ^ { k }$ with the solution. Here $k$ denotes the number of iterations and $B _ { l }$ are free variables in the subproblem. 

3. If there is no change of variables, stop the algorithm. Otherwise, go to Step 2. 

Theoretically the limit points of the sequence $\{ \mathbf { x } ^ { k } \}$ are only stationary points of $F$ , but they are also mostly local minimum of $F$ in practice. It is worth mentioning that the multi-convex condition of $F$ guarantees that the limit point exists. 

In step 2 of BCD, solving the convex subproblem exactly is usually time-consuming. Finding an approximate solution is a common way to accelerate the computation, for instance, applying only one step of gradient descent. The BCD method is thus categorized into two types: exact BCD and inexact BCD. Recently, Xu and Yin pointed out that inexact BCD usually yields lower objective values because the local linear approximation from one-step gradient descent helps avoid small regions around certain local minimum [2013]. 

The standard MIPS algorithm actually uses the exact BCD where the coordinates of each vertex form a block of variables and the MIPS energy function is locally convex with respect to each vertex on its one ring neighborhood [Hormann 2001]. The subproblem is solved by the Newton method. The inefficiency of MIPS is mainly because the local Newton iteration requires that it converges at every step, thus more iterations and evaluations of the Hessian are required. By replacing the Newton iterations with one step of gradient descent, MIPS can achieve lower distortions efficiently, as we show in Fig. 2. 

We present the details of the inexact BCD method for computing mesh-based locally injective mappings as follows: 

Blocks of variables In the original MIPS method, the coordinates of each vertex form a block of variables. Notice that the energy function defined on the simplex only involves neighboring vertices and the disconnected vertices actually can be updated simultaneously. We partition the mesh vertices into several blocks where any two vertices in the same block have no connected edges in the mesh. This partition is computed by the standard graph coloring algorithm. We use the Boost Graph library in our implementation. 

Initialization For mesh parameterization, we map the mesh onto a 2D circular domain by a convex combinatorial mapping so that the initial mapping is bijective. For mesh deformation and improvement, our input meshes are inverted-element free. 

Vertex updating For each vertex p in a block, we update it by one step of gradient descent: $\mathbf { p } _ { n e w } : = \mathbf { p } - \lambda \nabla _ { \mathbf { p } } F$ . The initial $\lambda$ is chosen such that $\mathbf { p } - \lambda \nabla _ { \mathbf { p } } F$ is on the boundary of the one ring region of p. $\lambda$ is reduced by $8 5 \%$ if $F$ increases or one of the neighboring simplices of p inverts. 

Stopping condition Our algorithm terminates when the relative error of the energy value is less than $1 0 ^ { - 6 }$ or the iteration number exceeds the maximal iteration number. 

Remark 1. In our inexact BCD implementation we update each vertex by one step of gradient descent. We also experimented with another inexact updating strategy – one step of Newton method: $\mathbf { p } _ { n e w } : = \mathbf { p } - \lambda \bigl ( \mathbf { \hat { V } _ { p } ^ { 2 } } F \bigr ) ^ { - 1 } \nabla _ { \mathbf { p } } \bar { F }$ . We found that the latter usually yields higher objective function values similar to using the exact Newton method. For instance, its energy curve of the example in Fig. 2 is almost coincident to the exact BCD. This behavior is probably because the approximation solution from one step of the Newton method is very close to the exact solution of the subproblem. So it is easily trapped by local minimum early. 

# 3.3 Advanced MIPS energy

As mentioned in Sec. 3.1, the MIPS energy conformal distortion and goes to infinity whe $E _ { m i p s } ^ { 2 D }$ penalizes 2Driangle degenerates. In 3D, the conformal distortion can be defined as the ratio between the maximal and minimal singular values. To avoid determining which singular value is maximal and minimal, we define a 

symmetric conformal energy as 

$$
\begin{array}{l} E _ {m i p s} ^ {3 D} (J (f)) = \frac {1}{8} \left(\frac {\sigma_ {1}}{\sigma_ {2}} + \frac {\sigma_ {2}}{\sigma_ {1}}\right) \left(\frac {\sigma_ {3}}{\sigma_ {1}} + \frac {\sigma_ {1}}{\sigma_ {3}}\right) \left(\frac {\sigma_ {2}}{\sigma_ {3}} + \frac {\sigma_ {3}}{\sigma_ {2}}\right) \\ = \frac {1}{8} \left(\| J (f) \| _ {F} ^ {2} \cdot \| J (f) ^ {- 1} \| _ {F} ^ {2} - 1\right) \\ = \frac {1}{8} \left(\kappa_ {F} (J (f)) ^ {2} - 1\right). \\ \end{array}
$$

Here $\sigma _ { 1 }$ , $\sigma _ { 2 }$ and $\sigma _ { 3 }$ are singular values of the Jacobian $J ( f )$ of the mapping $f$ , and $\kappa _ { F } ( J ( f ) )$ is the Frobenius condition number of $J ( f )$ . When $\sigma _ { 1 } = \sigma _ { 2 } = \sigma _ { 3 }$ , $E _ { m i p s } ^ { 3 D }$ is at its lowest value of 1 and $f$ is a conformal mapping. 

The conformal energy does not require the determinant of $J ( f )$ being 1 as is necessary for isometric mapping. We add an energy term $\begin{array} { r } { \dot { E _ { \mathrm { d e t } } } ( J ( f ) ) : = \frac { 1 } { 2 } ( \dot { \mathrm { d e t } } ( J ( f ) ) + \mathrm { d e t } ( \hat { J } ( f ) ) ^ { - 1 } ) } \end{array}$ into the conformal energy to define the isometric energy: 

$$
E _ {i s o} ^ {d} (J (f)) = \alpha E _ {m i p s} ^ {d} (J (f)) + (1 - \alpha) E _ {\det } (J (f)),
$$

where ments. $d$ represents 2D or 3D and o scale the minimal value $\alpha$ in our experi- to 1, we scale $E _ { m i p s } ^ { 2 D } ( J ( f ) )$ the standard 2D MIPS energy by $\frac { 1 } { 2 }$ . For mesh-based applications, $E _ { \mathrm { d e t } } ( J ( f ) )$ can be rewritten as $\textstyle { \frac { 1 } { 2 } }$ $\begin{array} { r } { \frac { 1 } { 2 } \bar { \Big ( } | T _ { i } | | t _ { i } | ^ { - 1 } + | t _ { i } | | T _ { i } | ^ { - 1 } \Big ) } \end{array}$ , where $| t _ { i } |$ is the volume of the transformed mesh element and $| \dot { T } _ { i } |$ is the volume of the original mesh element. $E _ { \mathrm { d e t } }$ is at the minimal value 1 when the mapping is equiareal/equivolume. Since a mapping is isometric if and only if it is both conformal and equiareal, $\dot { E } _ { i s o } ^ { d } ( J ( f ) )$ is a good measurement of isometric distortion. 

Maximal distortion control It is known that a single mesh element with high distortion ruins the parametrization result or yields unpleasing shape deformation even when the average distortion is very low. Therefore suppressing the maximum distortion is important in many computer graphics applications. Bounded distortion methods [Lipman 2012; Aigerman and Lipman 2013; Kovalsky et al. 2014; Poranne and Lipman 2014] set the upper bound of the distortion to explicitly control its maximum. The strict minimizer [Levi and Zorin 2014] tries to minimize the maximal distortion and control the distribution of distortion by optimizing the bound of the Frobenius norm of the difference between $J ( \bar { f } )$ and the rotation matrix formed by local frames. These approaches are costly due to the heavy use of (semi-infinite) quadratic programming. We propose to use an exponential function combined with $E _ { m i p s }$ or $E _ { i s o }$ to penalize the maximal conformal or isometric distortion directly, called the Advanced MIPS (AMIPS) energy. The conformal and isometric AMIPS energy functions are defined as follows: 

$$
\begin{array}{l} E _ {m i p s} ^ {\star} := \sum \exp (s \cdot E _ {m i p s}), E _ {d e t} ^ {\star} := \sum \exp (s \cdot E _ {\det }) \\ E _ {i s o} ^ {\star} := \sum \exp (s \cdot E _ {i s o}) = \sum \exp (s \cdot E _ {m i p s} + s \cdot E _ {\det}) \\ \end{array}
$$

where $\scriptstyle \sum$ sums over all mesh elements for mesh-based applications or the sample points for meshless deformations, and $s$ is a parameter controlling the level of penalty. A small $s$ has little effect on penalizing the maximal distortion, and a large $s$ causes numerical instability. We choose $s = 5$ for 2D mapping and $s = 2$ for 3D mapthat $E _ { \mathrm { d e t } }$ Forand $E _ { m i p s } ^ { 2 D }$ based locally injective mappings, it is known are locally convex with respect to any mesh vertex [Hormann we conjecture that $E _ { m i p s } ^ { 3 D }$ iao et al. 2011] in its one-ring region andis multi-convex too. The same conclusion holds for $E _ { i s o } ^ { \star }$ and $E _ { m i p s } ^ { \star }$ . So the inexact BCD can be effective minimizers (the gradient computation can be found in the appendix). By optimizing AMIPS energy functions, the maximum conformal and volume distortion can be suppressed effectively and the distribution 

of distortion can be tightly controlled. Fig. 1 compares our AMIPS with other state-of-the-art methods. 

Remark 2. Using the $L _ { \infty }$ norm to penalize the maximal distortion is an alternative way of optimization. But as observed in [Levi and Zorin 2014], $L _ { \infty }$ lacks the control of the distortion below the maximum. 

# 4 Experiments and Comparisons

We compare our AMIPS method with the state-of-the-art algorithms. For mesh parametrization and deformation, we compare with bounded distortion method (BDM) [Lipman 2012; Aigerman and Lipman 2013], locally injective mapping (LIM) [Schuller et al ¨ . 2013] and the strict minimizer (STRICT) [Levi and Zorin 2014]. For meshless deformation, we select the latest bounded distortion method (PGPM) [Poranne and Lipman 2014] as the competitor. For mesh quality improvement, we compare with the approaches presented in [Brewer et al. 2003] and [Fu et al. 2014]. We use the authors’ implementations in our experiments. For our inexact block coordinate descent method, we parallelized the update of vertices in each block of variables using OpenMP in $\mathrm { C } { + + }$ , which is about $3 { \sim } 5$ times faster than the serialized version. The order of blocks has little effect on our results according to our experiences. Our experiments were conducted on a desktop PC with a 3.4 GHz Intel Core i7 and 16 GB of RAM. 

Distortion metrics We measure the isometric distortion associated to a point p by $\delta _ { \mathbf { p } } ^ { i s o } : = \operatorname* { m a x } \{ \sigma _ { \operatorname* { m a x } , p } , \sigma _ { \operatorname* { m i n } , p } ^ { - 1 } \}$ [Sorkine et al. 2002]. $\delta ^ { i s o }$ reaches to the minimum 1 when all the singular values are 1. Note that $\delta _ { \mathbf { p } } ^ { i s o }$ is better than the common ARAP energy $\textstyle \sum _ { i } ( \sigma _ { i , p } - 1 ) ^ { 2 }$ which cannot measure the shrinking distortion well, especially when $\sigma _ { \mathrm { m i n } , p }$ is near to zero. 

The conformal distortion at $\mathbf { p }$ is measured by $\delta _ { \mathbf { p } } ^ { c o n } : = \sigma _ { \mathrm { m a x } , p } .$ σ min,p. $\sigma _ { \mathrm { m i n } , p } ^ { - 1 }$ Similarly the conformal distortion is denoted by $\delta _ { \mathbf { t } } ^ { c o n }$ for mesh-based mapping. When all the singular values are the same, $\delta ^ { c o n }$ attains its minimum of 1. 

The determinant of the Jacobian of the mapping also reflects the isometric distortion and we measure it by $\operatorname* { m a x } \{ \operatorname* { d e t } ( J ( f ) ) , \operatorname* { d e t } ( J ( f ) ) ^ { - 1 } \}$ . For mesh-based mapping, it is $\begin{array} { r }  \delta _ { \mathbf { t } } ^ { v o l } : = \operatorname* { m a x } \\{ \frac { \mathrm { v o l } \mathbf { t } } { \mathrm { v o l } \mathbf { t } ^ { 0 } } , \frac { \mathrm { v o l } \mathbf { t } ^ { 0 } } { \mathrm { v o l } \mathbf { t } } \} } \end{array}$ where $\mathbf { t } ^ { 0 }$ is the source image of t before mapping and vol computes the volume of the simplex. When $\delta ^ { v o l } = 1$ , the mapping is called equiareal/equivolume mapping. 

We report the above distortion metrics including their maximum (worse case), average and standard deviation over all the mesh elements or the sample points, denoted as $\delta _ { \mathrm { m a x } } ^ { i s o } , \delta _ { \mathrm { a v g } } ^ { i s o } , \delta _ { \mathrm { d e v } } ^ { i s o }$ $\delta _ { \mathrm { m a x } } ^ { c o n }$ , $\delta _ { \mathrm { a v g } } ^ { c o n }$ , $\delta _ { \mathrm { d e v } } ^ { c o n }$ , and δvolmax, δvolavg, δvoldev. We also visualize $\delta ^ { i s o }$ and $\delta ^ { c o n }$ on meshes using the color-coding patterns in Fig. 1-left and Fig. 4-upperleft. Best values are shown in bold. 

# 4.1 Planar parameterization

Isometric parameterization We map a single-patch mesh to a unit circular domain by convex combinatorial mapping first, then we minimize $E _ { i s o } ^ { \star }$ by inexact BCD. Considering that $E _ { \mathrm { d e t } } ^ { \star }$ could be considerably larger or smaller than the conformal energy $E _ { m i p s } ^ { \star }$ at the beginning of iteration, we scale the parameterization mesh by the ratio between the average edge length of the source mesh and the parameterization mesh per 100 iterations in the first 1000 iterations. Our method usually converges within 2000 iterations. 




Figure 3: Isometric parameterizations of a Bunny head and Bimba surface using AMIPS, LIM and BDM. The color encodes $\delta ^ { i s o }$ . Blue edges in the parameterization meshes highlight triangles that have a higher isometric distortion value than the maximal distortion value of our AMIPS method. The Bunny head is from Stanford 3D Scanning Repository and the Bimba model is from Aim@Shape Repository.


via a harmonic function using the authors’ implementation. LIM uses the same circular mesh as initialization. Since we have no prior on the best upper bound of BDM, we test the conformal bound 5 and 9 in our experiments. Table 1 shows the statistics of distortion metrics and timings. The parameterization quality of our method (AMIPS) is comparable or better than other methods with the least computational time. The textured models show that our results better preserve isometry. Although LIM can get a lower volume distortion $\delta ^ { v o l }$ than AMIPS, its conformal distortion is much higher (see the stretched regions on the back of Bunny head and the neck of Bimba for comparison). BDM bounds the conformal distortion, but the maximal area distortion is the largest among all the me(see the bottom regions of the two models). The small value of $\delta _ { d e v } ^ { i s o }$ 

Conformal parameterization Similar to the isometric parameterization, the initial ctriangles. Optimizing $E _ { m i p s } ^ { \star }$ map may introduce highly distorteddirectly usually requires more iterations. We apply our isometric parameterization in the first 1000 iterations to quickly expand the parameterization mesh away from the circular shape and then optimize $E _ { m i p s } ^ { \star }$ only. 

Fig. 4 compares our result with the standard MIPS, angle based flatting $\left( \mathrm { A B F + + } \right)$ ) [Sheffer et al. 2005], bounded distortion method (BDM) and locally injective mapping (LIM). We use our block definition, the same stopping condition and block updating scheme for standard MIPS. The conformal bound of BDM is the same as in [Lipman 2012] for the Isis and Fandisk models. We select the smallest feasible integer conformal bound for the other four models, and the LSCM [Levy et al ´ . 2002] energy is optimized for BDM. The solver of LIM also optimizes the LSCM energy. 

From the results, we can see that the standard MIPS is trapped by local minimum with higher distortion. $\mathrm { \ A B F + + }$ uses the least time for all results, but it introduces high maximal distortion, even producing inverted triangles.BDM bounds distortions if the algorithm converges. The bound is not easy to choose, for example, the Male model with bound 42 and Camel model with bound 27. Again, BDM is the slowest method due to its expensive computation. LIM can generate small average distortion, but its maximal distortion is high, for example, see the Dino model where the maximum distortion 29.26 is much higher than the average distortion 1.03. Its efficiency is also affected by the size of meshes. For the Bunny model with 35190 vertices, LIM takes 744 seconds to converge and becomes 

<table><tr><td>Model</td><td>\( \delta_{\max}^{iso}/\delta_{\text{avg}}^{iso}/\delta_{\text{dev}}^{iso} \)</td><td>\( \delta_{\max}^{con}/\delta_{\text{avg}}^{con}/\delta_{\text{dev}}^{con} \)</td><td>\( \delta_{\max}^{vol}/\delta_{\text{avg}}^{vol}/\delta_{\text{dev}}^{vol} \)</td><td>time (s)</td></tr><tr><td>Bunny head (AMIPS)</td><td>3.89/2.57/0.56</td><td>5.92/3.23/0.97</td><td>5.30/2.25/1.03</td><td>2.21</td></tr><tr><td>Bunny head (BDM(5))</td><td>18.18/3.64/2.08</td><td>5.00/3.74/0.87</td><td>169.82/4.62/6.51</td><td>47.87</td></tr><tr><td>Bunny head (BDM(9))</td><td>8.97/2.84/1.14</td><td>9.02/4.62/2.01</td><td>48.21/1.96/1.66</td><td>43.15</td></tr><tr><td>Bunny head (LIM)</td><td>6.29/2.39/0.99</td><td>39.07/6.09/5.32</td><td>11.62/1.11/0.35</td><td>41.19</td></tr><tr><td>Bimba (AMIPS)</td><td>4.62/2.20/0.47</td><td>6.18/2.44/1.03</td><td>4.47/2.16/0.57</td><td>4.87</td></tr><tr><td>Bimba (BDM(5))</td><td>7.83/2.59/1.10</td><td>5.00/3.32/1.12</td><td>12.50/2.16/1.43</td><td>106.48</td></tr><tr><td>Bimba (BDM(9))</td><td>6.99/2.42/0.97</td><td>9.00/4.06/2.22</td><td>6.28/1.52/0.54</td><td>89.49</td></tr><tr><td>Bimba (LIM)</td><td>6.79/2.29/0.92</td><td>27.71/4.05/3.83</td><td>2.84/1.55/0.29</td><td>16.31</td></tr></table>


Table 1: Isometric parameterization comparisons to BDM and LIM. A number in boldface emphasizes the best result.


to the slowest method. Our AMIPS achieves the lowest maximal distortion and average distortion among all the methods, and the computation speed is only slower than $\mathrm { \ A B F + + }$ . 

We also use the LIM framework to optimize AMIPS energy for conformal mapping. The gradient-descent and the Newton method are tested and the barrier term of LIM is enabled. We found the gradient-descent and the Newton method require a good initialization, otherwise they are easily trapped due to the nearly-zero linesearch step. For instance, for the Camel model, these methods fail if they start with the circular parameterization mesh. This phenomenon appears in all the examples in Fig. 4. We provide our intermediate AMIPS results to these methods as input and find that our result (using 400 inexact BCD iterations) can serve as a reliable initialization for the gradient-descent and the Newton method (see Fig. 5). 

# 4.2 Mesh Deformation

We apply our AMIPS algorithm on 2D triangle and 3D tetrahedral handle-based mesh deformations. We use the substepping strategy [Schuller et al ¨ . 2013] for effective optimization. We first determine the temporary handle positions as the farthest points along the ray from the original positions to the desired positions that does not trigger any simplex inversion, then optimize $E _ { i s o } ^ { \star }$ using our inexact BCD solver within a few iterations ( $3 { \sim } 7$ in our experiments) to reduce the AMIPS energy. We repeat the above procedure until handles are in their desired positions. Then we fix the handles and optimize $E _ { i s o } ^ { \star }$ until the algorithm converges. Our algorithm achieves interactive feedback for moderate sized meshes: 14 FPS for Fig 1 



Figure 5: Experiments on the gradient descent and the Newton method in the LIM framework. (a) is the result after 400 iterations of our method. (b) and (c) are conformal parametrization results using (a) as initialization. Their maximal and average conformal distortion are (10.14, 2.50), (3.96, 1.81) and (3.96, 1.19) respectively.


(5000 triangles), 63 FPS for Fig 6a (1267 triangles), 4 FPS for Fig 7b (13052 tets) and 8 FPS for Fig 7c (7537 tets). 

2D deformation Fig. 1 shows a simple 2D deformation. An interior vertex in a grid mesh is moved from the right side to the left side. We compare our AMIPS with BDM, LIM and STRICT. ARAP energy is used in LIM, BDM and STRICT methods. For the BDM method, we tried different isometric bounds and choose the lowest one such that BDM converges. The results show that AMIPS suppresses the maximal and average isometric distortions significantly. Furthermore it distributes the distortion more evenly than STRICT, which is designed to control the distortion distribution. 

Fig. 6 shows comparisons of two extreme deformation cases. The setting of all the methods is similar to Fig. 1. We can see that although the results of BDM [Lipman 2012] bounded the conformal distortion, their volume distortion is very high on the Circle model (see statistics in Table 2). The distortion of the LIM results is very high around the handles. Without explicitly bounding the isometric distortion, the STRICT method generates inverted triangles: 37 for Fig. 6a and 1828 for Fig. 6b. AMIPS is the fastest method and does not introduce high conformal and volume distortion. 

We further test our AMIPS energy using the LIM framework. We found that LIM-AMIPS can generate low-distortion results similar to AMIPS although it is less efficient. But compared to LIM with ARAP energy, its performance and the isometric quality are significantly better. Our experiments indicate that optimizing the AMIPS 




Figure 8: Extreme 3D deformation. The boundary vertices are fixed and three interior vertices are pushed inside. We show a cut-view of the mesh, and color encodes $\delta ^ { i s o }$ . Blue means the isometric distortion is greater than 100 and yellow shows inverted tets.


energy is better than optimizing ARAP energy with barrier. 

3D deformation Fig. 7 and 8 show four examples of 3D tetrahedral mesh deformation. The settings of LIM, BDM3D [Aigerman and Lipman 2013] and STRICT are similar to 2D deformation. The initial meshes of BDM3D are the results of LIM and the isometric bounds are chosen as low as possible so that BDM3D can converge. Our AMIPS method again achieves better results more efficiently (see statistics in Table 2). The deformation of the Cube model is a challenging case where three handles are pushed inside the cube. LIM tends to concentrate high distortion around the handles. BDM3D generates results similar to AMIPS but with more time. STRICT fails to find a reasonable result with 145 inverted tets. 



# 4.3 Meshless Deformation

Handle-based meshless deformation is popular in computer animation and image processing. A mapping $f : \mathbb { R } ^ { d } \overset { \cdot } {  } \mathbb { R } ^ { d }$ can be defined by a set of basis functions $B = \{ B _ { i } \} _ { i = 1 } ^ { m } \colon f ( { \bf x } ) : =$ $\begin{array} { r } { \mathbf { x } + \sum _ { i = 1 } ^ { m } \dot { \mathbf { c } _ { i } } B _ { i } ( \mathbf { x } ) } \end{array}$ . The Jacobian of $f$ at $\mathbf { x }$ has a simple form: $\begin{array} { r } { J ( f ) \overline { { = } } \tilde { I _ { d } } + \sum _ { i = 1 } ^ { m ^ { \ast } } \mathbf { c } _ { i } \nabla _ { \mathbf { x } } B _ { i } ( \mathbf { x } ) } \end{array}$ where $I _ { d }$ is a $d \times d$ identity matrix and the distortion metric can be built on the singular values of $J ( f )$ . Inspired by [Poranne and Lipman 2014], we formulate our AMIPS energy on a set of collocation points $\mathcal { Z } = \{ \mathbf { z } _ { j } \} _ { j = 1 } ^ { n }$ : $\begin{array} { r } { E _ { i s o } ^ { \star } : = \sum _ { j = 1 } ^ { n } E _ { i s o , { \bf z } _ { j } } ^ { \star } } \end{array}$ j=1to control the overall distortion. The ideal deformation leads to the following optimization problem: 

$$
\min  _ {\mathbf {c}} \gamma E _ {\mathbf {p o s}} + E _ {i s o} ^ {\star} \tag {2}
$$

where $\gamma > 0$ is a parameter weighting soft position constraint $E _ { \mathbf { p o s } }$ (the sum of squared differences between handles and their desired position) and updates during the iteration (the updating strategy is the same as in [Schuller et al ¨ . 2013]), and $\mathbf c \equiv [ \mathbf c _ { 1 } , . . . , \mathbf c _ { m } ]$ are unknowns initialized with zeroes. 

Minimization We test uniform cubic tensor product B-splines as basis functions. Six control points are used on each dimension. So c has only 72 elements in 2D. Since the number of variables is small, we can minimize Eq. 2 using the Newton method directly. When the number of collocation points is big, 20k for instance, the evaluation of the objective function and its gradient and Hessian is expensive. We create a hierarchical structure on the collocation points and optimize the objective function from coarse to dense to achieve interactive response rates (three layers with $1 \mathrm { k } / 5 \mathrm { k } / 2 0 \mathrm { k }$ points in our implementation), 63 FPS in Fig. 9. In 3D, the evaluation of the Hessian is more expensive (200 ms for Fig. 10 on the coarsest layer), so we use the conjugate gradient method instead to balance the convergence speed and computation. In the line-search step of the conjugate gradient and the Newton method, we backtrack the solution if $\operatorname* { d e t } ( J ( f ) )$ on any collocation point becomes negative. 

The soft position constraint energy $\gamma E _ { \mathbf { p o s } }$ dominates Eq. 2 when handles approach their desired positions. In some extreme situations, especially in 3D, the AMIPS energy $E _ { i s o } ^ { \star }$ gets higher and higher although it is considerably smaller than $\gamma E _ { \mathbf { p o s } }$ . The solver may fail to reduce AMIPS energy any more and is trapped by a local minimum. To solve this issue, we fix every handle $\mathbf { h } _ { k }$ at its current position, i.e. we enforce the hard constraint $\begin{array} { r } { \mathbf { h } _ { k } = \mathbf { h } _ { k } ^ { 0 } + \sum _ { i = 1 } ^ { m } \mathbf { c } _ { i } B _ { i } \mathbf { \bar { ( } } \mathbf { h } _ { k } ^ { 0 } \mathbf { ) } } \end{array}$ where $\mathbf h _ { k } ^ { 0 }$ is the rest pose of $\mathbf { h } _ { k }$ and substitute it into $E _ { i s o } ^ { \star }$ through Gaussian elimination, and optimize $E _ { i s o } ^ { \star }$ only. We usually take 10 iterations for the latter optimization on the coarsest layer for achieving interactive feedback. Then we restart to optimize Eq. 2. This procedure is repeated until all the handles are in their desired positions. 

ACM Transactions on Graphics, Vol. 34, No. 4, Article 71, Publication Date: August 2015 




Figure 13: Anisotropic tetrahedral mesh improvement. Slivers (whose dihedral angle is smaller than $1 5 ^ { \circ }$ in inverted space) are shown in red. Our algorithm eliminates all slivers.


of all these tetrahedrons while keeping the boundary lying on the original shape. We use inexact BCD to minimize the energy. 

Fig. 12 shows the improved all-hex meshes by our method. The original meshes (Dragon, Fertility, Rockarm) are generated by [Huang et al. 2014], [Gregson et al. 2011], and [Li et al. 2012] respectively and are already improved by Mesquite [Brewer et al. 2003]. Our method significantly improves the quality of hexes (see the zoom-in views and the scaled Jacobian of hexes in Table 3, where the best value is 1). 

Anisotropic tetrahedral mesh improvement Anisotropic meshing is targeted at generating tetrahedrons that are regular in the inversely transformed space under an anisotropic metric [Fu et al. 2014]. Note that the transformation is also a locally injective mapping. We utilize our AMIPS method to improve anisotropic meshes. We iteratively reposition mesh vertices by optimizing $E _ { i s o } ^ { \bar { \star } }$ and perform tetrahedron flip operations to reduce $E _ { i s o } ^ { \star }$ . Here we do not use the conformal energy $E _ { m i p s } ^ { \star }$ because the size of the desired tet is associated to the anisotropy and a similar transformation is not desired. Another common energy in tetrahedral mesh improvement is the inverse mean ratio metric combined with $E _ { \mathrm { d e t } } ( \bar { J ( f ) } )$ [Jiao 

et al. 2011]: $\begin{array} { r } { E _ { i m r m } ( J ( f ) ) = \frac { 1 } { 3 } \frac { \| J ( f ) \| _ { F } ^ { 2 } } { \operatorname* { d e t } ( J ( f ) ) ^ { 2 / 3 } } + E _ { \mathrm { d e t } } ( J ( f ) ) } \end{array}$ 13 kJ(f)k Fdet(J(f))2/3 + Edet(J (f )). We refer to the sum of $E _ { i m r m } ( J ( f ) )$ over all tetrahedrons as the inverse mean ratio metric (IMRM) isometric energy $E _ { i m r m }$ . Combined with the exponential function, we can define a new energy: $\begin{array} { r } { E _ { e i m r m } = \sum \exp ( s \cdot E _ { i m r m } ( J ( f ) ) ) } \end{array}$ called EIMRM energy. We compare our method with LCT [Fu et al. 2014] that improves meshes by random perturbation and 5-4 flips, using IMRM and EIMRM energy functions in our AMIPS framework. 



Figure 15: Isometric parameterization results using the default parameters (a) and the tuned parameters (b).


# 5 Conclusion and discussion

By revisiting the well-known MIPS method, we propose advanced MIPS (AMIPS), which is able to suppress the conformal/isometric distortion maximally using the inexact BCD solver. We demonstrate the efficiency and effectiveness of our method on various applications, including mesh parameterization and deformation, meshless deformation and mesh improvement. We believe our method can stimulate more effective methods for many applications. 

Proper initialization Unlike bounded distortion methods, our method requires that the optimization starts with a valid locally injective mapping similar to the standard MIPS and LIM. Since we choose a circular mapping to guarantee the validity of the initialization in our experiments, our method may take more iterations when the average distortion of the initialization is extremely high. For instance, our method took 2 hours for parameterizing a Bunny model with 1.2M vertices conformally. On the other hand, some linear approaches, like LSCM [Levy et al ´ . 2002] or LABF [Zayer et al. 2007], may provide a good initialization. Our method only took 2 minutes on the same Bunny model and yielded a similar good result when the initialization came from LABF (see Fig. 14). But it is known that linear approaches cannot always guarantee the local injectivity. Using mesh untangling techniques like [Escobar et al. 2003] to eliminate inverted elements is a promising preprocess for our method; however, there is also no theoretical guarantee on the success of mesh untangling. We leave this problem for future research. 

Parameter choices of AMIPS Two parameters $s$ and $\alpha$ are used in our AMIPS energy function. Their default values as given in Section 3.3 usually work well for general inputs. We only observed that using these default values may yield high distortion in mesh isometric parameterization if there are extremely high or low volume distortion in the initial circular mapping that causes numerical instability in computing the AMIPS energy and its gradient. Fig. 15a shows an example on a hand model with an open boundary on the 




wrist. Our algorithm with default parameters fails to achieve a low distortion result. To resolve this issue, we chose a small $s = 0 . 1$ at the beginning 1000 iterations to relieve the instability and reduce the distortion, then used the default value $s = 5$ for further improvement. Fig. 15b shows the improved result. We would like to extend this strategy in an automatic way in the future. 

Distortion bound Since our method is a purely nonlinear optimization scheme for computing locally injective mappings, we have no explicit control on maximal distortion, unlike the boundeddistortion (BDM) work of [Lipman 2012; Aigerman and Lipman 2013; Aigerman et al. 2014; Poranne and Lipman 2014]. However, determining an as low as possible bound is not trivial and the bisection strategy proposed in [Lipman 2012] is very costly both in time and memory. Results from our AMIPS provide a good estimate on the lowest upper bound of distortion and can be used as initialization of other methods (e.g. BDM) for further optimization. 

Bijectivity By explicitly detecting the intersection of the boundary of the parameterization mesh and using a line-search backtracking to avoid self-intersection, we can keep bijectivity of the parameterization during the optimization. The intersection detection usually slows down the computation, especially when the number of boundary vertices is large. For 

example, to compute the bijective conformal parameterization of the Camel model (see the insert), our algorithm is 20 times slower. This performance issue can be more severe in 3D deformation. 

Trajectory of handles It is known that there is no locally injective mapping in mesh deformation when the trajectories of handles are not consistent with each other [Jin et al. 2014]. As an extreme and hard case, a fixed boundary map treats all the boundary vertices as handles. Our AMIPS may fail in these situations if handle positions are set as hard constraints. Combining mesh retriangulation [Weber and Zorin 2014] or adapting mesh connectivity [Jin et al. 2014] is a possible way to resolve this problem. 




# Gradient computation of AMIPS

Here we only give the essential part of the gradient computation, the final result can be obtained easily by the chain rule. The MIPS energy $E = \| J ^ { - 1 } \| _ { F } \| J \| _ { F } .$ . Its derivative with respect to variable $p$ is: 

$$
\partial_ {p} E = 2 \operatorname {t r a c e} \left(\left(\frac {E}{\operatorname {t r a c e} (J ^ {T} J)} J ^ {T} - \frac {\operatorname {t r a c e} (J ^ {T} J)}{E} (J J ^ {T} J) ^ {- 1}\right) \cdot \partial_ {p} J\right)
$$

In 2D, it can be rewritten as 

$$
\partial_ {p} E = \frac {2 \operatorname {t r a c e} (J ^ {T} \cdot \partial_ {p} J)}{\det  J} - E \cdot \operatorname {t r a c e} (J ^ {- 1} \cdot \partial_ {p} J).
$$

The derivative of the volume distortion $\begin{array} { r } { E _ { \mathrm { d e t } } = \frac { 1 } { 2 } ( \operatorname* { d e t } J + \operatorname* { d e t } J ^ { - 1 } ) } \end{array}$ is: 

$$
\partial_ {p} E _ {\mathrm {d e t}} = \frac {1}{2} (\det  J - \det  J ^ {- 1}) \operatorname {t r a c e} (J ^ {- 1} \cdot \partial_ {p} J).
$$