using Test

include("TB097-111817-TPI-parte-1-codigo.jl")

@testset "Parámetros de Peng-Robinson (fluido puro)" begin
    # Lo corre para el Amoníaco, Perry 8th ed 2-502
    k_amoniaco = k(0.252608)
    alpha_amoniaco = alpha(405.65, 353.15, k_amoniaco)
    @test isapprox(alpha_amoniaco, 1.103, rtol=1e-3)

    a_amoniaco = a(405.65, 112.8*100000, alpha_amoniaco)
    @test isapprox(a_amoniaco, 4.611e-1*1.103, rtol=1e-3)
    
    b_amoniaco = b(405.65, 112.8*100000)
    @test isapprox(b_amoniaco, 23.262e-6, rtol=1e-3)
end
     
