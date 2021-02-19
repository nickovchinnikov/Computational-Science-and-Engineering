### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 9a052fc4-6c5e-11eb-24e7-3f1da02508bc
begin
	using Test
	using Pkg
	using LinearAlgebra
	
	Pkg.build("PyCall")
	Pkg.add("SymPy")
	Pkg.add("Images")
	Pkg.add("ImageIO")
	
	using SymPy
	using Images
	using ImageIO
	
	include("../../lib/KTBC.jl")
	import .KTBC.CreateKTBC
end

# ╔═╡ 437d7d1e-6c5e-11eb-0e0b-91d58c495478
md"""
## \#1
"""

# ╔═╡ 575d4bac-6c5e-11eb-060a-1dabfe5f167f
begin
	Aₜ = [
		-1  1 0
		 0 -1 1
		-1  0 1
	]
	
	Aₜᵀ = Aₜ'
	
	@test Aₜᵀ * Aₜ == [
		2 -1 -1
	   -1  2 -1
	   -1 -1  2
	]
	
	Aₛ = [
	   -1  1  0 0
	   -1  0  1 0
	   -1  0  0 1
		0 -1  0 1
		0  0 -1 1
	]
	
	Aₛᵀ = Aₛ'
	
	@test Aₛᵀ * Aₛ == [
		 3   -1  -1  -1
		-1    2   0  -1
		-1    0   2  -1
		-1   -1  -1   3
	]
end

# ╔═╡ 0119c440-6c5f-11eb-2c3d-f1319bccd42d
md"""
## \#2
"""

# ╔═╡ 03a129ce-6c5f-11eb-1883-ab1716ee30a3
begin
	C = sympy.Symbol("C")
	
	uₜ = ones(3, 1)
	
	@test convert(Array{Int}, Aₜ * (C * uₜ))[:] == [
		0
		0
		0
	]
	
	vₜ = [
		 1
		 1
		-1
	]
	
	@test convert(Array{Int}, Aₜᵀ * (C * vₜ)) == [
		0
		0
		0
	]
end

# ╔═╡ 77c8d8b0-6c5f-11eb-2623-61783b1d556b
md"""
## \#3
"""

# ╔═╡ 6a9bc3b2-6c61-11eb-11e8-9338d319794b
begin
	uₛ = ones(4, 1)
	
	@test convert(Array{Int}, Aₛ * (C * uₛ))[:] == [
		0
		0
		0
		0
		0
	]
	
	vₛ₁ = [
		1
		1
	   -2
		1
		1
	]
	
	@test convert(Array{Int}, Aₛᵀ * (C * vₛ₁)) == [
		0
		0
		0
		0
	]
	
	vₛ₂ = [
	   -2
		1
		1
	   -2
		1
	]
	
	@test convert(Array{Int}, Aₛᵀ * (C * vₛ₂)) == [
		0
		0
		0
		0
	]
end

# ╔═╡ 90d10c28-6c62-11eb-1f2e-b346125b3932
md"""
## \#8
"""

# ╔═╡ 9da4bf3a-6c6a-11eb-349a-153eec621c58
begin
	c₁, c₂, c₃, c₄, c₅, c₆ = symbols("c₁ c₂ c₃ c₄ c₅ c₆")
	
	K = [
		c₁+c₂+c₄       -c₁             -c₂           -c₄
	      -c₁        c₁+c₃+c₅          -c₃           -c₅
	      -c₂          -c₃           c₂+c₃+c₆        -c₆
	      -c₄          -c₅             -c₆         c₄+c₅+c₆
	]
	
	D = [
		c₁+c₂+c₄      0         0        0
	       0       c₁+c₃+c₅     0        0
	       0          0      c₂+c₃+c₆    0
	       0          0         0     c₄+c₅+c₆
	]
	
	W = [
		  0      c₁     c₂      c₄
	      c₁     0      c₃      c₅
	      c₂     c₃     0       c₆
	      c₄     c₅     c₆      0
	]
	
	@test K == D - W
end

# ╔═╡ 0a72fc88-6c63-11eb-0d6f-c790283ecd24
begin	
	a₁ = [
		1
	   -1
		0
		0
	]
	
	K₁ = [
		 c₁ -c₁ 0 0
		-c₁  c₁ 0 0
		 0   0  0 0
		 0   0  0 0
	]
	
	@test c₁ * a₁ * a₁' == K₁
	
	a₂ = [
		1
		0
	   -1
		0
	]
	
	K₂ = [
		c₂ 0 -c₂ 0
		0  0  0  0
	   -c₂ 0  c₂ 0
		0  0  0  0
	]
	
	@test c₂ * a₂ * a₂' == K₂
	
	K₃ = [
		0   0   0  0
		0   c₃ -c₃ 0
		0  -c₃  c₃ 0
		0   0   0  0
	]
	
	a₃ = [
		0
		1
	   -1
	    0
	]
	
	c₃ * a₃ * a₃'
	
	@test c₃ * a₃ * a₃' == K₃
	
	K₄ = [
	    c₄  0  0  -c₄
		0   0  0   0
		0   0  0   0
	   -c₄  0  0   c₄
	]
	
	a₄ = [
		1
		0
	    0
	   -1
	]
	
	@test K₄ == c₄ * a₄ * a₄'
	
	K₅ = [
	    0  0  0  0
		0  c₅ 0 -c₅
		0  0  0  0
		0 -c₅ 0  c₅	
	]
	
    a₅ = [
		0
		1
	    0
	   -1
	]
	
	@test K₅ == c₅ * a₅ * a₅'
	
	K₆ = [
         0  0  0   0
		 0  0  0   0
		 0  0  c₆ -c₆
		 0  0 -c₆  c₆
	]
	
	a₆ = [
		0
		0
	    1
	   -1
	]
	
	@test K₆ == c₆ * a₆ * a₆'
	
	@test K == K₁ + K₂ + K₃ + K₄ + K₅ + K₆
end

# ╔═╡ 23e4f7aa-6c6a-11eb-3fe1-8daa7a059e68
begin
	A₆ = [a₁ a₂ a₃ a₄ a₅ a₆]'
	
	@test A₆ == [
		 1  -1   0   0
 		 1   0  -1   0
 		 0   1  -1   0
 		 1   0   0  -1
 		 0   1   0  -1
 		 0   0   1  -1
	]
	
	C₆ = diagm(0 => [c₁, c₂, c₃, c₄, c₅, c₆])
	
	@test K == A₆' * C₆ * A₆
end

# ╔═╡ 03d87872-6c6c-11eb-1043-59a80e58642f
md"""
## \#9
"""

# ╔═╡ 13269750-6c6c-11eb-3262-2937ef0d3f5e
begin
	A₅ = [
		-1  1  0  0 0
		 0 -1  1  0 0
		 0  0 -1  1 0
		 0  0  0 -1 1
	]
	
	# 1) AᵀA == B₅
	B₅ = CreateKTBC(5)[3]
	
	A₅ᵀA₅ = A₅' * A₅
	
	@test B₅ == A₅ᵀA₅
	
	# 2) T₄ == (T₄ᵀT₄)_reduced
	
	A₄ = A₅[:,1:4]
	
	A₄ᵀA₄ = A₄' * A₄
	
	T₄ = convert(Array{Int}, CreateKTBC(4)[2])
	
	@test T₄ == A₄ᵀA₄
	
	# 3) inv(T₄)
	@test inv(T₄) == [
		4.0  3.0  2.0  1.0
 		3.0  3.0  2.0  1.0
 		2.0  2.0  2.0  1.0
 		1.0  1.0  1.0  1.0
	]
	
	# 4) eigenvalues of T₄
	values, vectors = eigen(T₄)
	
	round.(values, digits=2) == [
		0.12
		1.0
		2.35
		3.53
	]
	
	# 5) det(T₄)
	@test det(T₄) == 1.0
end

# ╔═╡ 59da7064-6c70-11eb-096b-4f532330b19c
md"""
## \#12
"""

# ╔═╡ 61271e4e-6c70-11eb-1dba-91041e4280fb
begin
	# Example
	dim = 5
	
	onesMatrix = ones(dim - 1, dim - 1)
	eyeMatrix = I(dim - 1)
	
	Kᵣ = diagm(0 => [dim for i=1:dim - 1]) - onesMatrix

	Kᵣinv = (1/dim) * (onesMatrix + eyeMatrix)
	
	round.(Int, Kᵣ * Kᵣinv) == eyeMatrix
	
	@test round(Int, det(Kᵣ)) == 125 # Positive definite
end

# ╔═╡ caa0ca62-6c77-11eb-1fc9-a79b32e5eeb3
md"""
```math
KK^{-1} = 
\begin{bmatrix}
   n-1 & -1 & \dots & -1 \\
   -1 & n-1 & \dots & -1 \\
   \vdots & \vdots & \ddots & \vdots \\
   -1 & -1 & \dots & n-1 \\
\end{bmatrix} \frac{1}{n}
\begin{bmatrix}
   2 & 1 & \dots & 1 \\
   1 & 2 & \dots & 1 \\
   \vdots & \vdots & \ddots & \vdots \\
   1 & 1 & \dots & 2 \\
\end{bmatrix} = 
```
##### K and inv K is (n - 1) * (n - 1)

```math
\frac{1}{n}
\begin{bmatrix}
   n & 0 & \dots & 0 \\
   0 & n & \dots & 0 \\
   \vdots & \vdots & \ddots & \vdots \\
   0 & 0 & \dots & n \\
\end{bmatrix} = I
```

##### Positive definite by det test

###### For n > 1
```math
det(n - 1) = n - 1 > 0
```

###### For 2 * 2, n > 2

```math
det\begin{pmatrix}
   n-1 & -1 \\
   -1 & n-1 \\
\end{pmatrix} = (n - 1)^2 - 1 = n^2 - 2n > 0
```

###### And etc

"""

# ╔═╡ d426ccea-6c7d-11eb-1d42-adde535f324c
md"""
## \#17
"""

# ╔═╡ 7e8ad1e6-6c7b-11eb-053b-01566d56c7f0
load("./17.png")

# ╔═╡ e8e78e24-6c7d-11eb-2a5f-7d9a9665ca06
begin
	A₉ = [
		1 1 0 0 0 0 0 0 0
		0 1 1 0 0 0 0 0 0
		0 0 1 0 0 1 0 0 0
		0 0 0 0 0 1 0 0 1
		0 0 0 0 0 0 0 1 1
		0 0 0 0 0 0 1 1 0
		0 0 0 1 0 0 1 0 0
		1 0 0 1 0 0 0 0 0
		0 0 0 1 1 0 0 0 0
		0 0 0 0 1 1 0 0 0
		0 1 0 0 1 0 0 0 0
		0 0 0 0 1 0 0 1 0
	]
	
	A₉ᵀA₉ = A₉' * A₉
	
	n, m = size(A₉ᵀA₉)
	
	# Non-zero values
	@test n * m - size(filter(iszero, A₉ᵀA₉))[1] == 33
	
	D₉ = diagm(0 => diag(A₉ᵀA₉))
	
	W₉ = A₉ᵀA₉ - D₉
	
	@test W₉ == [
		 0  1  0  1  0  0  0  0  0
		 1  0  1  0  1  0  0  0  0
		 0  1  0  0  0  1  0  0  0
		 1  0  0  0  1  0  1  0  0
		 0  1  0  1  0  1  0  1  0
		 0  0  1  0  1  0  0  0  1
		 0  0  0  1  0  0  0  1  0
		 0  0  0  0  1  0  1  0  1
		 0  0  0  0  0  1  0  1  0
	]
	
	@test D₉ == [
		 2  0  0  0  0  0  0  0  0
		 0  3  0  0  0  0  0  0  0
		 0  0  2  0  0  0  0  0  0
		 0  0  0  3  0  0  0  0  0
		 0  0  0  0  4  0  0  0  0
		 0  0  0  0  0  3  0  0  0
		 0  0  0  0  0  0  2  0  0
		 0  0  0  0  0  0  0  3  0
		 0  0  0  0  0  0  0  0  2
	]
	
	# node 5 has 4 edges connected to it
	@test D₉[5, 5] == 4
	
	# nodes 2, 4, 6, 8 connection in the row
	@test W₉[5:5, :] == [0  1  0  1  0  1  0  1  0]
end

# ╔═╡ Cell order:
# ╠═9a052fc4-6c5e-11eb-24e7-3f1da02508bc
# ╠═437d7d1e-6c5e-11eb-0e0b-91d58c495478
# ╠═575d4bac-6c5e-11eb-060a-1dabfe5f167f
# ╠═0119c440-6c5f-11eb-2c3d-f1319bccd42d
# ╠═03a129ce-6c5f-11eb-1883-ab1716ee30a3
# ╠═77c8d8b0-6c5f-11eb-2623-61783b1d556b
# ╠═6a9bc3b2-6c61-11eb-11e8-9338d319794b
# ╠═90d10c28-6c62-11eb-1f2e-b346125b3932
# ╠═9da4bf3a-6c6a-11eb-349a-153eec621c58
# ╠═0a72fc88-6c63-11eb-0d6f-c790283ecd24
# ╠═23e4f7aa-6c6a-11eb-3fe1-8daa7a059e68
# ╠═03d87872-6c6c-11eb-1043-59a80e58642f
# ╠═13269750-6c6c-11eb-3262-2937ef0d3f5e
# ╟─59da7064-6c70-11eb-096b-4f532330b19c
# ╠═61271e4e-6c70-11eb-1dba-91041e4280fb
# ╟─caa0ca62-6c77-11eb-1fc9-a79b32e5eeb3
# ╠═d426ccea-6c7d-11eb-1d42-adde535f324c
# ╠═7e8ad1e6-6c7b-11eb-053b-01566d56c7f0
# ╠═e8e78e24-6c7d-11eb-2a5f-7d9a9665ca06
