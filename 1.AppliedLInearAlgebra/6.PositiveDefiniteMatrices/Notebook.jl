### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 1f74c496-6726-11eb-3894-119ff752210e
begin
	using SymPy;
	using Test;
	using LinearAlgebra;
	using Plots;
	using Calculus;
	using PlutoUI;
	
	include("../../lib/KTBC.jl");
	using .KTBC: CreateKTBC;
end

# ╔═╡ 0668a58e-6727-11eb-2b6c-41626284bd0d
md"""
## Init deps
"""

# ╔═╡ 1c3315a0-6774-11eb-0081-f328cca176ed
begin
	import Pkg
	Pkg.add("PlutoUI")
	Pkg.add("SymPy")
	Pkg.add("Test")
	Pkg.add("LinearAlgebra")
	Pkg.add("Plots")
	Pkg.add("Calculus")
end

# ╔═╡ 9020f922-6725-11eb-2bc5-c190c647a78a
md"""
# Task \#1
### Symbol calculations
"""

# ╔═╡ 4acbb302-6726-11eb-1af2-0157c1a07336
begin
	@vars u₁ u₂
	u⁺ = [u₁ u₂]
	u = [
	  u₁
	  u₂
	]
	
	T₂ = convert(Array{Int}, CreateKTBC(2)[2])
	
	u⁺T₂u = u⁺ * T₂ * u
end

# ╔═╡ e6a566ba-6726-11eb-0ef2-e148ca198a89
md"""
#### Find the result
"""

# ╔═╡ 1cac6756-6727-11eb-2d9a-3776545b32da
result₁ = [u₁ * (u₁ - u₂) + u₂ * (-u₁ + 2 * u₂)]

# ╔═╡ 263872ca-6727-11eb-1233-c37ef9371211
md"""
### The result is:

`` u^T \cdot T₂ \cdot u = u₁^2 - 2u₁u₂ + u₂^2 ``

### Or you can see it by diff

`` u^T \cdot T₂ \cdot u - (u₁^2 - 2u₁u₂ + u₂^2) = 0 ``

#### If result is equal, difference is equal zero

"""

# ╔═╡ 70f8a80e-6727-11eb-0e9c-833f93a4c9f5
@test u⁺T₂u[1] - result₁[1] == 0

# ╔═╡ 7ded0a64-6727-11eb-1dd7-a353f7a80ba7
md"""
# Task \#3
### Produce Circulant matrix C
"""

# ╔═╡ 965b12bc-6727-11eb-0099-ddd1176df887
A = [
  1 -1 0
  0 1 -1
  -1 0 1
]

# ╔═╡ 9d00019a-6727-11eb-2a6f-3fda51c74b01
C = (A') * A

# ╔═╡ b3db8876-6727-11eb-1e64-1fa9949feee3
z = zeros(Int, 3)

# ╔═╡ c2c35d50-6727-11eb-1f1e-1161a06fa745
md"""
#### A contains dependent columns

```math
\begin{bmatrix}
  a_{11} \\
  a_{21} \\
  a_{31}
 \end{bmatrix} + \begin{bmatrix}
  a_{12} \\
  a_{22} \\
  a_{32}
 \end{bmatrix} + \begin{bmatrix}
  a_{13} \\
  a_{23} \\
  a_{33}
 \end{bmatrix} = \begin{bmatrix}
  0 \\
  0 \\
  0
\end{bmatrix}

```

"""

# ╔═╡ 8b584ee2-6728-11eb-2015-17f4ce49fc83
@test A[1:3] + A[4:6] + A[7:9] == z

# ╔═╡ 93d52f1a-6728-11eb-2733-cd6a61aa2b89
md"""
#### The same for C

```math
 \begin{bmatrix}
  c_{11} \\
  c_{21} \\
  c_{31}
 \end{bmatrix} + \begin{bmatrix}
  c_{12} \\
  c_{22} \\
  c_{32}
 \end{bmatrix} + \begin{bmatrix}
  c_{13} \\
  c_{23} \\
  c_{33}
 \end{bmatrix} = \begin{bmatrix}
  0 \\
  0 \\
  0
 \end{bmatrix}
```
"""

# ╔═╡ cdc18b56-6728-11eb-1b95-c340bfb9d559
@test C[1:3] + C[4:6] + C[7:9] == z

# ╔═╡ e2ccc9c8-6728-11eb-11c0-7fe79dc8a96d
md"""
#### Find the vector u
"""

# ╔═╡ f2cb5eac-6728-11eb-1d1f-fbb306febaef
begin
	v = ones(Int, 3)
	@test A * v == z
end

# ╔═╡ 06eb6980-672a-11eb-0d64-5dc376fb893b
md"""
#### C -> is not Positive definite for Cholesky
"""

# ╔═╡ 1eff6daa-672a-11eb-2df1-b54942784659
@test_throws PosDefException cholesky(C)

# ╔═╡ 2f1b2972-672a-11eb-1df4-11ba022cdd62
md"""
# Task \#11
### Find symmetric matrix for function
`` f(x,y) = 2xy ``
"""

# ╔═╡ 6a32b714-672a-11eb-2e95-d99358a79a68
md"""
#### Partial derivatives table `` f(x, y) = 2xy ``

```math
	\begin{matrix}
		\frac{\mathrm \partial}{\mathrm \partial x}
			\left( f(x, y) \right) = 2y
		&
		\frac{\mathrm \partial}{\mathrm \partial y}
			\left( f(x, y) \right) = 2x
	\\
		\frac{\mathrm \partial}{\mathrm \partial x \partial y}
			\left( f(x, y) \right) = 2
	 	&
		\frac{\mathrm \partial}{\mathrm \partial y \partial x}
			\left( f(x, y) \right) = 2
	\\
		\frac{\mathrm \partial^2}{\mathrm \partial x^2}
			\left( f(x, y) \right) = 0
		&
		\frac{\mathrm \partial^2}{\mathrm \partial y^2}
			\left( f(x, y) \right) = 0
	\\
 \end{matrix}
```

#### Hessian matrix

```math

\begin{bmatrix}
	\frac{\mathrm \partial}{\mathrm \partial x}
		\left( f(x, y) \right)
	&
	\frac{\mathrm \partial}{\mathrm \partial x \partial y}
		\left( f(x, y) \right)
	\\
	\frac{\mathrm \partial}{\mathrm \partial y \partial x}
		\left( f(x, y) \right)
	&
	\frac{\mathrm \partial}{\mathrm \partial y}
		\left( f(x, y) \right)
\end{bmatrix} \cdot \begin{bmatrix}
	x
	\\
	y
\end{bmatrix}

= 

\begin{bmatrix}
	2 & 0
	\\
	0 & 2
\end{bmatrix} \cdot \begin{bmatrix}
	x
	\\
	y
\end{bmatrix}

= 

\begin{bmatrix}
	2x
	\\
	2y
\end{bmatrix}

```

"""

# ╔═╡ 60f181b4-672a-11eb-16d8-bd6b5b5725cd
f(x, y) = 2 * x * y

# ╔═╡ 561285de-672a-11eb-33f9-3d9265291867
@vars x y a b c

# ╔═╡ 476a3532-6779-11eb-1ff1-9d670dda615f
md"""
#### Let's init symbolic ``S^{symbolic}=S^s``
"""

# ╔═╡ 71a3c0f2-6779-11eb-2fd1-cb4e54b8cfa8
begin
	w⁺ = [x y]
	w = [
	  x
	  y
	]
	Sˢ = [
		a b
		b c
	]
	w⁺ * Sˢ * w
end

# ╔═╡ 01eec5e6-672b-11eb-1bda-b38590341dac
md"""

#### Let's find a, b, c


```math
	\begin{matrix}
		x \left(a x + b y\right) + y \left(b x + c y\right)
		&
		=
		&
		ax^2 + 2byx + cy^2
	\\
		ax^2 = cy^2 = 0
		&
		=>
		&
		a = c = 0
	\\
		2byx = 2yx
		&
		=>
		&
		b = 1
 \end{matrix}

```

"""

# ╔═╡ 4eee9e0a-672a-11eb-2ab1-39f67d4ff891
S = [
  0 1
  1 0
]

# ╔═╡ eacb5054-672a-11eb-0cc7-4b6b9661e804
@test (w⁺ * S * w) - [2*x*y] == [0]

# ╔═╡ 817df440-672c-11eb-16ef-eb398ab32a0c
begin
	vals, vects = eigen(S)
	@test vals == [
		-1.0
		1.0
	]
end

# ╔═╡ c2f65df8-6782-11eb-34eb-a914ea97d400
md"""
### Task \# 24
"""

# ╔═╡ d40dd6e8-6782-11eb-2460-6dec703edafc
begin
	K = CreateKTBC(5)[1]
end

# ╔═╡ 17e154ac-6789-11eb-1936-7176d5492f77
round.(inv(K), digits=2)

# ╔═╡ 17620080-678c-11eb-142c-0be670b9ec02
md"""
``K`` is symmetric as ``K^{-1}``

``K^{-1} = (K^{-1})^{T}``
"""

# ╔═╡ 28d4134e-6789-11eb-1561-dfad67ceb90c
@test round.(inv(K), digits=2) == round.(inv(K)', digits=2)

# ╔═╡ Cell order:
# ╟─0668a58e-6727-11eb-2b6c-41626284bd0d
# ╠═1c3315a0-6774-11eb-0081-f328cca176ed
# ╠═1f74c496-6726-11eb-3894-119ff752210e
# ╟─9020f922-6725-11eb-2bc5-c190c647a78a
# ╠═4acbb302-6726-11eb-1af2-0157c1a07336
# ╟─e6a566ba-6726-11eb-0ef2-e148ca198a89
# ╠═1cac6756-6727-11eb-2d9a-3776545b32da
# ╟─263872ca-6727-11eb-1233-c37ef9371211
# ╠═70f8a80e-6727-11eb-0e9c-833f93a4c9f5
# ╟─7ded0a64-6727-11eb-1dd7-a353f7a80ba7
# ╠═965b12bc-6727-11eb-0099-ddd1176df887
# ╠═9d00019a-6727-11eb-2a6f-3fda51c74b01
# ╠═b3db8876-6727-11eb-1e64-1fa9949feee3
# ╟─c2c35d50-6727-11eb-1f1e-1161a06fa745
# ╠═8b584ee2-6728-11eb-2015-17f4ce49fc83
# ╟─93d52f1a-6728-11eb-2733-cd6a61aa2b89
# ╠═cdc18b56-6728-11eb-1b95-c340bfb9d559
# ╟─e2ccc9c8-6728-11eb-11c0-7fe79dc8a96d
# ╠═f2cb5eac-6728-11eb-1d1f-fbb306febaef
# ╟─06eb6980-672a-11eb-0d64-5dc376fb893b
# ╠═1eff6daa-672a-11eb-2df1-b54942784659
# ╟─2f1b2972-672a-11eb-1df4-11ba022cdd62
# ╟─6a32b714-672a-11eb-2e95-d99358a79a68
# ╠═60f181b4-672a-11eb-16d8-bd6b5b5725cd
# ╠═561285de-672a-11eb-33f9-3d9265291867
# ╟─476a3532-6779-11eb-1ff1-9d670dda615f
# ╠═71a3c0f2-6779-11eb-2fd1-cb4e54b8cfa8
# ╟─01eec5e6-672b-11eb-1bda-b38590341dac
# ╠═4eee9e0a-672a-11eb-2ab1-39f67d4ff891
# ╠═eacb5054-672a-11eb-0cc7-4b6b9661e804
# ╠═817df440-672c-11eb-16ef-eb398ab32a0c
# ╟─c2f65df8-6782-11eb-34eb-a914ea97d400
# ╠═d40dd6e8-6782-11eb-2460-6dec703edafc
# ╠═17e154ac-6789-11eb-1936-7176d5492f77
# ╟─17620080-678c-11eb-142c-0be670b9ec02
# ╠═28d4134e-6789-11eb-1561-dfad67ceb90c
