### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ b9f84364-6b75-11eb-183f-9358dcc0e4d0
begin
	using Test
	using LinearAlgebra
end

# ╔═╡ 9c44adc6-6b7c-11eb-20b5-11ad83fee84b
ArrayFloatOrIntDim2 = Array{<:Union{Float64, Int64}, 2}

# ╔═╡ 95a22bc2-6b77-11eb-2cd5-c50658a942b9
md"""
## Gram - Schmidt
"""

# ╔═╡ 3d8a1668-6b70-11eb-293c-ab5ce07a3ce5
begin
	function GramSchmidt(A::ArrayFloatOrIntDim2)
		m, n = size(A)
		Q = zeros(m, n)
		R = zeros(n, n)
		for j = 1:n
			v = A[:,j]
			for i = 1:j-1
				R[i, j] = Q[:, i]' * v
				v = v - R[i, j] * Q[:, i]
			end
			R[j, j] = norm(v)
			Q[:,j] = v / R[j,j]
		end
		return (Q, R)
	end
end

# ╔═╡ 7356cb64-6b71-11eb-01f7-91f369a9fb74
begin
	A = [
		1 0
		1 1
		1 3
		1 4
	]
	Q, R = GramSchmidt(A)
	
	@test A == round.(Int, Q * R)
end

# ╔═╡ 78c7863c-6b7c-11eb-0be7-cbfb3c368065
md"""
## Householder (probably wrong)
"""

# ╔═╡ 83f357d4-6b7c-11eb-11b6-5ddf77041fe5
begin
	function householder(A::ArrayFloatOrIntDim2)
		A = convert(Array{Float64}, A)
		m, n = size(A)
		U = zeros(m, n)
		for k = 1:n
			w = A[k:m, k]
			wₙ = norm(w)
			w[1] = (w[1] - wₙ)
			u = w / wₙ
			U[k:m, k] = u
			A[k:m, k:n] = A[k:m, k:n] - 2 * u * (u' * A[k:m, k:n])
		end
		R = triu(A[:, 1:n])
		return U, R
	end
end

# ╔═╡ 43133fa8-6b7d-11eb-331e-73d6eabd7d83
begin
	U, Rₕ = householder(A)
end

# ╔═╡ Cell order:
# ╠═b9f84364-6b75-11eb-183f-9358dcc0e4d0
# ╠═9c44adc6-6b7c-11eb-20b5-11ad83fee84b
# ╟─95a22bc2-6b77-11eb-2cd5-c50658a942b9
# ╠═3d8a1668-6b70-11eb-293c-ab5ce07a3ce5
# ╠═7356cb64-6b71-11eb-01f7-91f369a9fb74
# ╟─78c7863c-6b7c-11eb-0be7-cbfb3c368065
# ╠═83f357d4-6b7c-11eb-11b6-5ddf77041fe5
# ╠═43133fa8-6b7d-11eb-331e-73d6eabd7d83
