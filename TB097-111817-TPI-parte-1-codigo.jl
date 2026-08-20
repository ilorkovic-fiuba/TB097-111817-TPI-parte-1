# --------- Variables termodinámicas de cada componente ---------
struct Componente
    nombre::String
    Tc::Float64 # K
    Pc::Float64 # Pa
    w::Float64  
end

struct MezclaTernaria
    nombre::String
    componentes::Vector{Componente}
    kij::Matrix{Float64}
end

metano = Componente("Metano", 190.564, 4599000.0, 0.0115)
etano = Componente("Etano", 305.32, 4872000.0, 0.0995)
propano = Componente("Propano", 369.83, 4248000.0, 0.1523)

kij = [ 0.0    -0.003  0.016;
     -0.003   0.0    0.001; 
      0.016   0.001  0.0  ]

gas_natural = MezclaTernaria("Gas Natural", [metano, etano, propano], kij)

# Despues agrego los componentes polares

# --------- Parámetros de Peng-Robinson ---------

function k(w::Float64)
    return 0.37464 + 1.54226*w - 0.26992* (w)^2
end

function alpha(Tc::Float64, T::Float64, k::Float64)
    return ( 1 + k*( 1 - ( T/Tc )^(1/2) ) )^2
end

function a(Tc::Float64, Pc::Float64, alpha::Float64)
    return 0.45724 * ( 8.314462618^2 * Tc^2 / Pc ) * alpha
    # Uso en todo el código Sistema Internacional (8.31 J/mol.k)
end

function b(Tc::Float64, Pc::Float64)
    return 0.0778 * 8.314462618 * Tc / Pc
end

# --------- Parámetros de mezcla de Peng-Robinson ---------

function b_mix(y::Vector{Float64}, componentes::Vector{Componente})
    bi = [b(c.Tc, c.Pc) for c in componentes]
    return sum(y .* bi)
end

function a_mix(y::Vector{Float64}, T::Float64, mezcla::MezclaTernaria)
    comps = mezcla.componentes
    k_ij = mezcla.kij
    n = length(comps)

    a = [a(c.Tc, c.Pc, alpha(c.Tc, T, k(c.w))) for c in comps]

    a_m = 0.0
    for i in 1:n
        for j in 1:n
            aij = sqrt(a[i]*a[j])*(1 - k_ij[i, j])
            a_m += y[i]*y[j]*aij 
        end
    end
end

