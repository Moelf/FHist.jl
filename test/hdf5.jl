using Test, FHist, HDF5


@testset "HDF5 I/O" begin
    fname = tempname()
    h1 = Hist1D(rand(10), 0:0.1:1)
    h2 = Hist2D((rand(10), rand(10)), (0:0.1:1, 0:0.2:1))
    h3 = Hist3D((rand(10), rand(10), rand(10)), (0:0.1:1, 0:0.2:1, 0:0.3:1))
    h3_int = Hist3D((rand(1:10, 10), rand(1:5, 10), rand(1:9999, 10)), (0:10, 0:5, 0:3); overflow=true)
    h5writehist(fname, "h1", h1)
    h5writehist(fname, "a/h2", h2)
    h5writehist(fname, "a/b/h3", h3)
    h5writehist(fname, "a/b/c/h3_int", h3_int)

    h5open(fname) do f
        version = VersionNumber(read_attribute(f["h1"], "_h5hist_version"))
        @test v"1.0" == version
    end

    # Explicit type specifications
    h1_from_h5 = h5readhist(fname, "h1", Hist1D)
    h2_from_h5 = h5readhist(fname, "a/h2", Hist2D)
    h3_from_h5 = h5readhist(fname, "a/b/h3", Hist3D)
    h3_int_from_h5 = h5readhist(fname, "a/b/c/h3_int", Hist3D)
    @test h1 == h1_from_h5
    @test h2 == h2_from_h5
    @test h3 == h3_from_h5
    @test h3_int == h3_int_from_h5

    # Without type specification
    h1_from_h5 = h5readhist(fname, "h1")
    h2_from_h5 = h5readhist(fname, "a/h2")
    h3_from_h5 = h5readhist(fname, "a/b/h3")
    h3_int_from_h5 = h5readhist(fname, "a/b/c/h3_int")
    @test h1 == h1_from_h5
    @test h2 == h2_from_h5
    @test h3 == h3_from_h5
    @test h3_int == h3_int_from_h5
end
