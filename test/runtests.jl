using Test
import LogBeta

@testset "incomplete logbeta" begin
    @test LogBeta.logbeta(1.0, 2.0, 0.2, 0.9) ≈ -1.15518264015650398955737181354
    @test LogBeta.logbeta(3.0, 2.0, 0.01, 0.99) ≈ -2.48550282746660144854598013619
    @test LogBeta.logbeta(3.0, 2.0, 0.001, 0.01) ≈ -14.9226584212253042301433768824
    @test LogBeta.logbeta(3.0, 2.0, 0.99, 0.999) ≈ -9.92703257898043949664632355807
    @test LogBeta.logbeta(3.0, 2.0, 0.2, 0.9) ≈ -2.56774492809700481550963695770
    @test LogBeta.logbeta(3.0, 2.0, 0.1, 0.2) ≈ -6.23566150762002408487996529558
    @test LogBeta.logbeta(3.0, 2.0, 0.9, 0.95) ≈ -5.74770170859103654995762029784
end
