# 18.330-Problem-set-9
18.330 Problem set 9


**Download Link:https://programming.engineering/product/18-330-problem-set-9/**


Description
5/5 – (2 votes)
Exercise 1: Fourier series

In lectures we saw that we can write a periodic function in the form

̂

( ) = ∑

=−

where ̂ = 1 ∫2 ( ) − d if the period is 2 .

2 0

The FFTW library uses a different convention, taking the period to be 1. Sup-pose we sample a function at + 1 uniformly-spaced points between [0, 1], obtaining values 0, 1 , … , , where = ( ) with ℎ = 1 . By using the Fourier series we are implicitly assuming that is periodic, i.e. 0 = . Then

2 −1

̂2

( ) = ∑

=− 2

where the Fourier coefficients are

1

̂ = ∫ ( ) −2 d

0

This is the convention we will use in this problem.

Discretize the integral for the Fourier series coefficients using the trape-zoidal rule using + 1 points. Using the assumption that 0 = reduce this to a sum over for = 0, … , − 1.

Implement this as a function fourier_coefficients_trapezoidal(f::Vector) where f is the vector of samples of for = 0, … , − 1. Your function should return the coefficients for ̂for = − 2 , … , 2 − 1 as

a vector.

Consider the following three functions:

The sawtooth function

( ) = mod( , 1)

(The function mod( , 1) returns just the fractional part of a number.)


2. The “W” wave function

( ) =

mod

0≤ <0.5

{1 −

0.5≤ <1

{

( ,1))

(

otherwise

3. The function

ℎ( ) = exp(cos(2 ))

Plot them.

Given the properties of these functions, how would you expect their Fourier coefficients to decay?

3. Take

= 100

. Calculate

,

̂

ℎ

using your fourier_ceofficients_trapezoidal

and

≥ 0

only for

function. Plot the magnitude of the coefficients as a function of

. Do you see the

expected behavior?

̂

̂

How accurately does our function fourier_coefficients_trapezoidal(f::Vector) calculate the Fourier coefficients? Use the following analytic solutions to

calculate and plot the respective errors:

0̂

1

and

̂

for

≠ 0

0̂= 4

;

=

2

2

=

2

̂=− 2

;

̂= 0

1

1

for

odd

and

for

even

ℎ̂= 1( )

where 1( ) is the modified Bessel function of the first kind, calculated in Julia as besseli(n, 1) using the SpecialFunctions.jl package.

Calculate the error (using norm) as the number of points used changes between = 10 and = 1000.

5. Write

a

function

̂

reconstruct_fourier_series(fs::Vector,

xs::Vector) which reconstructs f(x) from the Fourier series coef-ficients.

̂

Make a plot of ||f(xt) – reconstruct_fourier_series(fs, xt)||/N

as the number of coefficients used

is varied from

4→200

for each

of the functions.

Sample the reconstructed Fourier series at the points xt = 0:0.001:1.

What do you see? Can you explain why?


Make an interactive visualiztion that plots the following on the same axes:

the points used to calculate the Fourier coefficients in fourier_coefficients_rectangle;

the reconstructed function from the Fourier coefficients, found using reconstruct_fourier_series; and

the true function as the number of points is varied from= 10 ∶

2 ∶ 250.

Does this help explain the results in [1.6]? In particular, what do you see for the sawtooth function? This is known as the Gibbs phenomenon, which occurs when a function is discontinuous.

What is the operation count for your naive fourier_coefficients_trapezoidal function? In general this will not behave well as grows very large.

The FFT, however, can calculate this in (log( )) steps. The FFT im-plemented in FFTW.jl calculates

−1

̂= ∑ −2 /

=0

for = 0 ∶ − 1, which should be related to your discretization in [1.1].

Note that the s are different here. But since −2 ( + ) / = −2 −2 / = −2 / we have the relationship /2+ ̂ = −̂ /2+ .

The FFT algorithm therefore outputs

f ̂= [ 0̂, 1̂,… , /2̂−1, −̂ /2, −̂ /2+1, … , −̂1]

The FFTWpackage defines a function fftshift that shifts this vector to the form

f ̂= [ −̂ /2, −̂ /2+1, … , −̂1, 0̂, 1̂,… , /2̂−1]

Implement a function fast_fourier_coefficients that outputs the same results as fourier_coefficients_trapezoidal but using the FFT from FFTW.jl.

Check that the output is the similar to before.

Time your two functions for= 210. Is one faster than the other? How large can you take N such that fast_fourier_coefficients runs for under 1 second? What about for fourier_coefficients_trapezoidal?


Exercise 2: Solving an ODE with a spectral method

In this problem we will solve the boundary-value problem

″ =

with boundary conditions (0) = (1) = 0 by using a spectral method, i.e. by expanding in suitable basis functions.

Due to the chosen boundary conditions we will consider the sine series

∞

( ) = ∑ sin( ),

=1

which is the same as the Fourier series for an odd function.

In practice we have to truncate the summation as

( ) = ∑ sin( )

=1

The coefficients are given by

1

= 2 ∫ ( ) sin( )d

0

Similarly to the DFT we discretize this using the rectangle rule to get the Discrete Sine Transform (DST):

−1

sin

=

2

(

)

∑ ( )

=1

This sum is implemented in julia using the FFTW.jl library (which you will need to install) as follows:

dst(x) = FFTW.r2r(x, FFTW.RODFT00) / length(x)

Its inverse is given by

idst(x) = FFTW.r2r(x, FFTW.RODFT00)

[r2r stands for “real to real”, meaning that the transform maps a real vector to a real vector. RODFT00 is a symbol that selects one particular type of transform.]

1. Assume that there is a solution of the form

( ) = ∑

sin

( )

.

̃

Substitute this into the ODE

=1

″

2

̂

2

( ) = ( )

to show that

̃= −

.


2. We can therefore solve the ODE for by first calculating

using the DST,

̂

using the iDST.

̂

then calculate , and finally invert

Write a function spectral_solver(b) that does this to solve the ODE, where b is the discretized version of g(x). Solve the ODE with ( ) = sin(2 ). Plot the result.

The right-hand side is given by

h = 1 / N

b = sin.(2π * (h:h:1-h))

Calculate the error as a function of and plot it. What rate of conver-gence do you see?

Generate the error plots again, now for the right-hand side ( ) = exp (sin (2 )) − 1. Use the solution from your spectral solver with = 213 as the true solution. Calculate the error for = 2 for = 1 → 12. Some care is needed to make sure you use the correct points from the “true” solution for comparison.

Exercise 3: Finding roots in a different way

In class we defined the Chebyshev expansion of a function as

( ) = ∑ ( )

=0

which is an th-degree polynomial. The Chebyshev polynomials are defined

as ( ) = cos( arccos( )).

In general for a smooth function the Chebyshev series converges rapidly. We therefore expect that the roots of ( ) should be close to the roots of ( ), provided that is indeed a good approximation to . We have already seen that we can find all the roots of a polynomial with various methods.

The Chebyshev polynomials satisfy the following recurrence relation:

+1( ) = 2 ( ) − −1( )

We will use this to find the companion matrix for ( ), from which we can find the roots of ( ).

Consider the polynomial ( ). This is a degree + 1 polynomial and hence can be reexpanded in Chebyshev polynomials.


Consider the vector of Chebyshev polynomials

T

0

( )

1

( )

Now we can write

T

( ) =

⋮

, where

is an

T

−1( )

( )

that

is a

( ) = ( ) +

̂

( −1)×( −1)

matrix. We need the

term to account for the fact

−1

degree-

̂

polynomial. Here,

is the standard basis

vector with zeros everywhere except in the

th component.

Use the recurrence relation to find the form of and .

Verify what you found in [3.1] numerically when = 10. Build the vector T( ) for a random number in [−1, 1]. Compute and check that it gives T( ).

This looks almost like an eigenvalue problem, except for theterm.

should you choose so that

−1 ( ) ̂

from the right-hand

To remove this we can add and subtract

side. Writing

( ) =( ) + ∑ =0( )

, what value of

T

T

( ) + ( )

̂ −

−1

(1)

( ) = T

̂ ∑( )

(2)

̃

( ) +

( )

=0

=

̂

̃

What is the new matrix ?

This becomes an eigenvalue problem when is a root of ( ). There-fore, the eigenvalues of ̃are the roots of .

Write a function buildM(c::Vector) that constructs the matrix ̃from the coefficients in the Chebyshev expansion. Use this to write a cheby-shev_roots(c) function that finds the roots of the polynomial defined us-ing the Chebyshev coefficients c. Finally write a function fN(x, c) that calculates the series expansion to find ( ) defined by the vector c.

We can calculate Chebyshev coefficients using the dct functions in FFTW. We will use the Chebyshev points = cos( / ) for = 0 ∶ . You can then calculate the Chebyshev coefficients using the following code:


chebyshev_points(N) = cos.(π*(0:1:N)/N)

function chebyshev_coefficients(x)

N = length(x)

c = FFTW.r2r(x, FFTW.REDFT00)/(N-1)

c[1] /= 2

c[N] /=2

return c

end

Consider the polynomial ( ) = ( − 1/2)2 ( 2 − 1/9) ( + 1/4).

Using 10 Chebyshev points calculate the Chebyshev coefficients and then calculate the roots using chebyshev_roots. What do you see. What about multiplicities?

Plot ( ) and scatter plot the roots you find on top.

Now consider solving the problem exp(cos(4 )) = 1. Using = 100 points, alculate the Chebyshev coefficients for ( ) = exp(cos(4 ))− 1. Do they decay quickly? Use these to calculate the roots of ( ). Plot ( ) and scatter the roots you find on top. Do you find all the roots?

(Hint: you will find 100 eigenvalues. Only plot those that are real and lie between -1 and 1.)

Make an interactive visualization as is varied between = 4 ∶ 150. Plot ( ), the Cheyshev approximation to ( ) using coefficients and the roots you find on the same axes. Comment on what you see.

At what value ofdo you find all the roots? Remember to plot only those roots that are real and between -1 and 1.

Exercise 4: Gram–Schmidt for polynomials

In lectures we discussed treating the set of polynomials {1, , 2, 3, …} as the basis of a vector space with the inner product

( , ) = ∫ ( ) ( ) ( )d

We can therefore carry out Gram–Schmidt orthogonalization on these polyno-mials to generate a family of orthogonal polynomials.

We will implemnent this using the Polynomials.jl package. Integrals can be performed using the polyint(f, a, b) function to integrate the polynomial


from to . Here, , and should all be Polynomials.

Write a function gram_schmidt(vs::Vector, ip) which accepts a vec-tor of “vectors” in the vector space and the inner product on the vector spacce. For standard vectors this would be dot(v1, v2). The function should implement the Gram–Schmidt algorithm and return a vector of the resulting orthonormal basis elements.

Test your function for standard vectors vs = [rand(10) for i = 1:10] using ip = dot. To check everything went according to plan, form the matrix = ⋅ . If everything worked this should be the identity matrix.

For polynomials we define the inner product

( , ) = ∫ ( ) ( ) ( )d

For Legendre polynomials we have = −1 , = 1, and ( ) = 1. Using the Polynomials.jl package, implement a function legen-dre_inner_product(f, g). Use this to orthogonalize the vector of mono-mials up to order 7.

Use the functions you found in [4.3] and your inner-product function, find the Legendre polynomial expansion coefficients for the function ( ) =

( − 1/2)2 ( 2 − 1/9)( + 1/4).

Plot the reconstructed polynomial using the first coefficients for = 1 ∶ 7 and the true function. What do you see? How good are the Legen-dre polynomials at approximating the function?

8

