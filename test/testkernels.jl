# Testing kernel functions

@testset "Kernels" begin
    boxcar = OnlineMoments.Boxcar()
    @test boxcar(0.0) == 0.5
    @test boxcar(0.5) == 0.5
    @test boxcar(0.99) == 0.5
    @test boxcar(-0.99) == 0.5
    @test boxcar(-1.0) == 0.0
    @test boxcar(1.0) == 0.0

    epaneknikov = OnlineMoments.Epaneknikov()
    @test epaneknikov(0.0) ≈ 5*3*sqrt(5)/100
    @test epaneknikov(sqrt(5)) ≈ 0
    @test epaneknikov(sqrt(5)-0.01) > 0
    @test epaneknikov(sqrt(5)+0.01) == 0
    @test epaneknikov(sqrt(5)+1.0) == 0
    @test epaneknikov(-sqrt(5)+0.01) > 0
    @test epaneknikov(-sqrt(5)-0.01) == 0
    @test epaneknikov(-sqrt(5)-1.0) == 0
end