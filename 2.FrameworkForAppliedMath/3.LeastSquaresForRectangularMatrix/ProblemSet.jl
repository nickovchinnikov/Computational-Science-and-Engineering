### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ c4381b3c-6b8f-11eb-34c1-9760c2fc6d87
begin
	using Test
end

# ╔═╡ 968a3130-6b8e-11eb-1d6c-1317dbf89bc8
md"""
## \#7
"""

# ╔═╡ b4d8b1a2-6b8e-11eb-3273-bd06e5e16a9b
begin
	b = [
		4
		1
		0
		1
	]
	
	x = [
		0
		1
		2
		3
	]
	
	A = [
		ones(4, 1) x
	]
	
	u = round.(Int, inv(A' * A) * A' * b)
	
	@test u == [
		3
		-1
	]
end

# ╔═╡ d2278240-6b90-11eb-0bc3-df340513edd4
md"""
## \#8
"""

# ╔═╡ e6c51852-6b90-11eb-1e64-1b29b65808fd
begin
	p = A * u
	
	e = b - p
	
	@test A' * e == [
		0
		0
	]
end

# ╔═╡ 00970928-6b93-11eb-1707-f157198becdd
md"""
## \#12
"""

# ╔═╡ 09012d82-6b93-11eb-00ab-41d5b7a8615e
begin
	A₂ = [
		A x.^2
	]
	
	u₂ = round.(Int, inv(A₂' * A₂) * A₂' * b)
	
	@test u₂ == [
		4
		-4
		1
	]
	
	p₂ = A₂ * u₂
	
	e₂ = b - p₂
	
	# Exact solution
	e₂ == [
		0.0
		0.0
		0.0
		0.0
	]
	
	A₃ = [
		A₂ x.^3
	]
	
	u₃ = round.(Int, inv(A₃' * A₃) * A₃' * b)
	
	p₃ = A₃ * u₃
	
	e₃ = b - p₃
	
	# Exact solution for cubic
	e₃ == [
		0.0
		0.0
		0.0
		0.0
	]
end

# ╔═╡ a50a2640-6ba0-11eb-1247-2388f1b38906
md"""
## \#22
"""

# ╔═╡ aa83f89c-6ba0-11eb-32ae-8ba9f30015eb
begin
	A₄ = [
		1 -1
		1 1
		1 2
	]
	
	b₄ = [
		5
		13
		17
	]
	
	u₄ = round.(Int, inv(A₄' * A₄) * A₄' * b₄)
	
	@test u₄ == [
		9
		4
	]
	
	e₄ = A₄ * u₄ - b₄
	
	@test e₄ == [
		0
		0
		0
	]
end

# ╔═╡ Cell order:
# ╠═c4381b3c-6b8f-11eb-34c1-9760c2fc6d87
# ╟─968a3130-6b8e-11eb-1d6c-1317dbf89bc8
# ╠═b4d8b1a2-6b8e-11eb-3273-bd06e5e16a9b
# ╟─d2278240-6b90-11eb-0bc3-df340513edd4
# ╠═e6c51852-6b90-11eb-1e64-1b29b65808fd
# ╟─00970928-6b93-11eb-1707-f157198becdd
# ╠═09012d82-6b93-11eb-00ab-41d5b7a8615e
# ╠═a50a2640-6ba0-11eb-1247-2388f1b38906
# ╠═aa83f89c-6ba0-11eb-32ae-8ba9f30015eb
