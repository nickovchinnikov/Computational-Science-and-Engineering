### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 22f805e8-7a7b-11eb-05e8-3dd4205997e1
begin
	using Test
	using LinearAlgebra
end

# ╔═╡ edb0b6e8-7a7a-11eb-0a88-990e1165fc1c
md"""
## \#2
"""

# ╔═╡ 357c1f92-7a7b-11eb-22bf-cd420daded2c
begin
	A = [
		1 -1 1 -1
		0  1 2  3
		0  0 2  6
		0  0 0  6
	]
	inv(A) * [
		 0
		-1//2
		-1
		-1
	]
end

# ╔═╡ d1f578ee-7a82-11eb-1375-c314d8672dd4
md"""
## \#3
"""

# ╔═╡ d5880030-7a82-11eb-3b88-a90d14844ed0
begin
	A₂ = [
		 1  1  1 1
		-1  1 -1 1
		 3  2  1 0
		 3 -2  1 0
	]
	inv(A₂) * [
		-1//6
		 0
		-1//2
		 0
	]
end

# ╔═╡ Cell order:
# ╠═22f805e8-7a7b-11eb-05e8-3dd4205997e1
# ╠═edb0b6e8-7a7a-11eb-0a88-990e1165fc1c
# ╠═357c1f92-7a7b-11eb-22bf-cd420daded2c
# ╠═d1f578ee-7a82-11eb-1375-c314d8672dd4
# ╠═d5880030-7a82-11eb-3b88-a90d14844ed0
