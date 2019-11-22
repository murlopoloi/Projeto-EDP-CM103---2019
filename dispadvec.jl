function dispadvec(β,ϵ,L,T,σ)
  #=
     β é a velocidade da água no rio;
     ϵ é o coeficiente de dispersão longitudinal;
     L é o comprimento do rio;
     T é o tempo observado de decaimento;
     σ é o coeficiente de decaimento de poluentes.
  =#

  σN = 0.3 * sqrt(1/2)
  μ = L/3
  # σN e μ coeficientes da distribuição normal.

  Δx = min(ϵ / abs(β),L/100)
  Δt = 0.5 / (2*ϵ / Δx^2 + σ)
  # condições de estabilidade do modelo.

  n = ceil(Int, L / Δx) # número de partições do espaço.
  m = ceil(Int, T / Δt) # número de partições do tempo.

  println("Δx = $Δx, n = $n")
  if n > 1000 || m > 10000
    error("NAO")
  end

  u = zeros(n+1,m+1)
  u0, uL = 0, 0  # condições de contorno.
  u[1,:].= u0
  u[n+1,:].= uL
  x = Δx
  f(x) = 2 * exp(-(x - L)^2) # fonte externa.

  for i = 2:1:n
    u[i,1] =  10/(σN*(2*pi)^(1/2)) * exp(-0.5 * ((x - μ)/σN)^2) # distribuição normal.
    x += Δx
    #μ = ceil(Int, μ - L/2)  # iteração utilizada apenas no exemplo 2.
  end

  for i = 2:1:n
    for j = 1:1:m
      x = (i - 1) * Δx
      u[i,j+1] = u[i,j] + Δt * (f(x) - σ * u[i,j] - (β * (u[i+1,j] - u[i-1,j])/(2 * Δx)) + ϵ * ((u[i-1,j]-2 * u[i,j] + u[i+1,j])/(Δx^2)))
    end
  end
  return u
end
