### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ ec794b0c-67b1-11eb-23e7-459641028308
begin
	using Markdown
	using InteractiveUtils
	using Test
	using LinearAlgebra
end

# ╔═╡ 47f1758c-67b1-11eb-0a89-27d2595e5912
md"""

## Task \# 8

Compute the norms and condition numbers

#### Norm

```math
||A||^2 = max \frac{||Ax||^2}{||x||^2} = max \frac{ x^T A^T A x }{ x^T x } = \lambda_{max}(A^T A) = \sigma^{2}_{max}
```

```math
||A|| = \sqrt{\sigma^{2}_{max}} = \sigma_{max}
```

```math
||A^{-1}|| = \frac{1}{\sigma_{min}}
```

"""

# ╔═╡ 9aa05ff6-67b5-11eb-0b36-95986372b105
function normCalculator(A::Array{Int})
	AᵀA = A' * A
	
	A_eigenvals, A_eigenvec = eigen(AᵀA)
	
	A_normFromEigenSquared = √max(A_eigenvals...)
	
	A_norm = opnorm(A)
	
	return (A_normFromEigenSquared, A_norm)
end

# ╔═╡ 4533fdcc-67b3-11eb-30f8-851ef7fc43c9
begin
	A₁ = [
		1 7
		1 1
	]
	
	A₁normFromEigenSquared, A₁norm = normCalculator(A₁)
	
	@test round.(A₁normFromEigenSquared, digits=5) == round.(A₁norm, digits=5)
	
	A₁normFromEigenSquared, A₁norm, cond(A₁)
end

# ╔═╡ bf65578a-67b8-11eb-2627-9139089fcf68
md"""
#### Cond is Inf, because of `` || A_2^{-1} || = \sigma_{min} = 0 ``
"""

# ╔═╡ 1c2ed1ac-67b5-11eb-1e77-67e96dbd2802
begin
	A₂ = [
		1 1
		0 0
	]
	
	A₂normFromEigenSquared, A₂norm = normCalculator(A₂)
	
	@test round.(A₂normFromEigenSquared, digits=5) == round.(A₂norm, digits=5)
	
	A₂normFromEigenSquared, A₂norm, cond(A₂)
end

# ╔═╡ 1568a204-67b9-11eb-3d55-d72928803e5e
md"""
#### Cond is 1, because of `` || A_2 || = || A_2^{-1} || ``
"""

# ╔═╡ ed7d86f6-67b3-11eb-273c-43770e6f299c
begin
	A₃ = [
		1 1
		-1 1
	]
	
	A₃normFromEigenSquared, A₃norm = normCalculator(A₃)
	
	@test round.(A₃normFromEigenSquared, digits=5) == round.(A₃norm, digits=5)
	
	A₃normFromEigenSquared, A₃norm, cond(A₃)
end

# ╔═╡ c3dc9efe-67bd-11eb-0224-af9d70653b73
md"""
### Task \# 13
Estimate the condition number
"""

# ╔═╡ e21e095c-67bd-11eb-12e4-e92b18b0b223
begin
	A₄ = [
		1 1
		1 1.0001
	]
	cond(A₄)
end

# ╔═╡ 67060602-67c2-11eb-046e-dfc437a3dc81
md"""
### Task \#17
Fibonacci matrix and the SVD
"""

# ╔═╡ 832167be-67c2-11eb-1693-ab8ea1a74995
begin
	F = [
		1 1
		1 0
	]
	svd(F)
end

# ╔═╡ 8d7f475a-67c3-11eb-243b-53c5d68004f1
md"""
### Task \#18
"""

# ╔═╡ 9a4b6df6-67c3-11eb-27d6-2546dd023fb4
begin
	K = [
		1 2
		2 5
	]
	lu(K)
end

# ╔═╡ Cell order:
# ╠═ec794b0c-67b1-11eb-23e7-459641028308
# ╟─47f1758c-67b1-11eb-0a89-27d2595e5912
# ╠═9aa05ff6-67b5-11eb-0b36-95986372b105
# ╠═4533fdcc-67b3-11eb-30f8-851ef7fc43c9
# ╟─bf65578a-67b8-11eb-2627-9139089fcf68
# ╠═1c2ed1ac-67b5-11eb-1e77-67e96dbd2802
# ╟─1568a204-67b9-11eb-3d55-d72928803e5e
# ╠═ed7d86f6-67b3-11eb-273c-43770e6f299c
# ╟─c3dc9efe-67bd-11eb-0224-af9d70653b73
# ╠═e21e095c-67bd-11eb-12e4-e92b18b0b223
# ╟─67060602-67c2-11eb-046e-dfc437a3dc81
# ╠═832167be-67c2-11eb-1693-ab8ea1a74995
# ╟─8d7f475a-67c3-11eb-243b-53c5d68004f1
# ╠═9a4b6df6-67c3-11eb-27d6-2546dd023fb4
