# Testing utility functions

@testset "Utility functions" begin
    X_stream = stream_data(X_small)
    @test X_small[1] == X_stream()
    @test X_small[2] == X_stream()
    @test X_small[3] == X_stream()
end