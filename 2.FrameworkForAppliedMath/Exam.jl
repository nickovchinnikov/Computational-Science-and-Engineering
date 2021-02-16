### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 5590424a-6d0b-11eb-13a5-ed80e8342a06
begin
	using Test
	using LinearAlgebra
end

# ╔═╡ 1d03f5fe-6d10-11eb-1918-c1478c31b083
md"""
[Exam pdf link](https://ocw.mit.edu/courses/mathematics/18-085-computational-science-and-engineering-i-fall-2008/exams/quiz1_18085f07.pdf)
"""

# ╔═╡ 0cac1572-6d10-11eb-2699-f1d5f5baa4e3
md"""
## \#1
"""

# ╔═╡ 67a6401a-6d0b-11eb-1803-756ee33e062d
begin
	A₀ = [
		-1  1  0 0
		 0 -1  1 0
		 0  0 -1 1
	]
	
	A₀ᵀ = A₀'
	
	A₀ᵀA₀ = A₀ᵀ * A₀
	
	@test A₀ᵀA₀ == [
	  1  -1   0   0
	 -1   2  -1   0
	  0  -1   2  -1
	  0   0  -1   1
	]
	
	@test A₀ᵀA₀ * [
		1
		1
		1
		1
	] == [
		0
		0
		0
		0
	]
	
	@test det(A₀ᵀA₀) == 0
	
	A₁ = [
		 1  0 0
		-1  1 0
		 0 -1 1
	]
	
	A₁ᵀA₁ = A₁' * A₁
	
	@test A₁ᵀA₁ == [
		 2  -1   0
 		-1   2  -1
  		 0  -1   1
	]
	
	@test det(A₁ᵀA₁) == 1
end

# ╔═╡ 3fb30142-6d10-11eb-2a9c-939b9a40a2da
md"""
## \#2
"""

# ╔═╡ 41e88e3c-6d10-11eb-27ea-254fd36ce6cb
begin
	B = [
		 1 -1  0
		-1  2 -1
		 0 -1  1
	]
	
	values, Q = eigen(B)
	
	Λ = diagm(0 => values)
	
	@test B == round.(Int, Q * Λ * Q')
end

# ╔═╡ Cell order:
# ╠═5590424a-6d0b-11eb-13a5-ed80e8342a06
# ╟─1d03f5fe-6d10-11eb-1918-c1478c31b083
# ╠═0cac1572-6d10-11eb-2699-f1d5f5baa4e3
# ╠═67a6401a-6d0b-11eb-1803-756ee33e062d
# ╠═3fb30142-6d10-11eb-2a9c-939b9a40a2da
# ╠═41e88e3c-6d10-11eb-27ea-254fd36ce6cb
