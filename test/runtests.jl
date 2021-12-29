using CosmoParticles
using Rotations
using SortingAlgorithms
using StatsBase
using Test
using Unitful

const CP = CosmoParticles

@testset "CosmoParticles.jl" begin

    @testset "Particles" begin
        p = Particles(:dm)

        p.id = sample(1:1000, 100; replace=false)
        p.pos = rand(3, 100)

        # different id and pos for copied Particles
        pcid = rand(100)
        pcpos = rand(3, 100)

        # Particle specific
        pc = Particles(:dm, :id => p.id, :pos => p.pos)
        @test pc == p

        pc = copy(p)
        @test pc == p
        @test pc !== p
        @test pc.props == p.props
        @test pc.props !== p.props

        pc = empty(p)
        @test pc.type === :dm
        @test isempty(pc.props)
        @test isempty(pc)

        @test CP.particle_name(p) == "Particles"
        @test issetequal(propertynames(p), [:type, :id, :pos])
        @test issetequal(propertynames(p; private=true), [:type, :props, :id, :pos])

        io = IOBuffer()
        show(io, "text/plain", p)
        @test String(take!(io)) == "dm: 100 Particles\n id pos"

        # test AbstractParticles with basic implementation Particles
        @test CP.get_props(p) === p.props

        @test getproperty(p, :id) === p.props[:id]
        @test getproperty(p, :type) === p.type === :dm
        @test getproperty(p, :props) === p.props
        @test p.id === p.props[:id]

        pc = copy(p)
        setproperty!(p, :pos, pcpos)
        @test p.pos === pcpos
        pc.id = pcid
        @test pc.id === pcid
        @test p.id !== pcid

        @test_throws ErrorException("""
            setfield!: immutable struct of type Particles cannot be changed""") p.type = :gas

        @test p[:id] === p.props[:id]

        pc = copy(p)
        p[:id] = pcid
        @test p.id === pcid

        @test keys(p) == keys(p.props)
        @test haskey(p, :id)
        @test !haskey(p, :vel)
        @test !haskey(p, :type)

        @test values(p) == values(p.props)

        @test issetequal(propertynames(p), [:type, :id, :pos])

        # test default implementations for AbstractParticles
        @test issetequal(Base.@invoke(propertynames(p::AbstractParticles)), [:id, :pos])

        pc = deepcopy(p)
        empty!(pc)
        @test isempty(pc.props)
        @test isempty(pc)

        @test !isempty(p)

        @test CP.particle_number(p) == 100
        @test CP.particle_number(Particles(:gas)) == 0
        @test CP.particle_number(Particles(:gas, Dict{Symbol,Any}(:mass => 3))) == 0

        @test isnothing(Base.@invoke(CP.particle_name(p::AbstractParticles)))

        io = IOBuffer()
        CP.show_properties(io, "text/plain", p)
        @test String(take!(io)) == "100 Particles\n id pos"
    end


    @testset "Utils" begin
        @testset "Index particles" begin
            p = Particles(:dm)
            p.id = sample(1:1000, 100; replace=false)
            p.pos = rand(3, 100)
            p.mass = 4u"kg"

            mask = rand(100) .> 0.5
            ind = findall(mask)

            @test CP._applyind(3, mask) == 3
            @test CP._applyind(p.mass, mask) == p.mass
            @test CP._applyind(p.id, mask) == p.id[mask]
            @test CP._applyind(p.pos, mask) == p.pos[:, mask]
            @test CP._applyind(3, ind) == 3
            @test CP._applyind(p.mass, ind) == p.mass
            @test CP._applyind(p.id, ind) == p.id[ind]
            @test CP._applyind(p.pos, ind) == p.pos[:, ind]

            pc = CP.applyind(p, mask)
            @test pc.id == p.id[mask]
            @test pc.pos == p.pos[:, mask]
            @test pc.mass == p.mass

            @test CP.applyind(p, ind) == pc
            @test p[ind] == pc
            @test p[mask] == pc

            pcc = deepcopy(p)
            CP.applyind!(pcc, mask)
            @test pcc == pc

            pcc = deepcopy(p)
            CP.applyind!(pcc, ind)
            @test pcc == pc

            pc = CP.applyind(p, mask; affect=(:pos, :mass))
            @test !haskey(pc, :id)
            @test pc.pos == p.pos[:, mask]
            @test pc.mass == p.mass
        end

        @testset "Find all in" begin
            a = sample(1:10000, 6000; replace=false)
            set = sample(1:10000, 100; replace=false)
            inboth = intersect(a, set)

            ind = CP.findall_in(a, set)
            @test length(ind) == length(inboth)
            @test issetequal(a[ind], inboth)
            @test CP.findall_in(a, Set(set)) == ind

            sort!(a)
            sort!(set)
            ind = CP.findall_in(a, set)
            @test length(ind) == length(inboth)
            @test issetequal(a[ind], inboth)
            @test CP.findall_in(a, Set(set)) == ind
            @test CP.findall_in_sorted(a, set) == ind

            # test empty arrays
            @test CP.findall_in_sorted([], []) |> isempty
            @test CP.findall_in_sorted(a, []) |> isempty
            @test CP.findall_in_sorted([], set) |> isempty

            @test CP.findall_in([], Set([])) |> isempty
            @test CP.findall_in(a, Set([])) |> isempty
            @test CP.findall_in([], Set(set)) |> isempty
        end
    end


    @testset "Transformations" begin
        @testset "Rotation" begin
            rotmat = rand(RotMatrix{3})
            a = rand(3, 100)
            au = a * u"m"

            arot = rotmat * a
            aurot = rotmat * au

            @test CP.matrix_rotate(a, rotmat) == arot
            @test CP.matrix_rotate(au, rotmat) == aurot

            ac = copy(a)
            CP.matrix_rotate!(ac, rotmat)
            @test ac == arot

            auc = copy(au)
            CP.matrix_rotate!(auc, rotmat)
            @test auc == aurot


            p = Particles(:dm)
            p.pos = copy(au)
            p.vel = copy(a)

            pc = rotate(p, rotmat, :pos)
            @test p.pos == au
            @test pc.pos == aurot
            @test p.vel === pc.vel == a

            pc = rotate(p, rotmat)
            @test p.pos == au
            @test pc.pos == aurot
            @test p.vel == a
            @test pc.vel == arot


            pc = deepcopy(p)
            rotate!(pc, rotmat, :vel)
            @test pc.pos == au
            @test pc.vel == arot

            pc = deepcopy(p)
            rotate!(pc, rotmat)
            @test pc.pos == aurot
            @test pc.vel == arot
        end

        @testset "Translations" begin
            da = rand(3)
            dau = rand(3) * u"m"
            a = rand(3, 100)
            au = a * u"m"

            atrans = a .+ da
            autrans = au .+ dau

            p = Particles(:dm)
            p.pos = copy(au)
            p.vel = copy(a)

            pc = translate(p, dau)
            @test p.pos == au
            @test pc.pos == autrans
            @test p.vel === pc.vel == a

            pc = translate(p, da, :vel)
            @test p.pos === pc.pos == au
            @test p.vel == a
            @test pc.vel == atrans

            pc = deepcopy(p)
            translate!(pc, dau)
            @test pc.pos == autrans
            @test pc.vel == a

            pc = deepcopy(p)
            translate!(pc, da, :vel)
            @test pc.pos == au
            @test pc.vel == atrans
        end
    end


    @testset "Misc Operations" begin
        ids = sample(1:1000, 100; replace=false)
        a = rand(3, 100)
        bu = rand(100) * u"kg"

        p = Particles(:dm)
        p.id = ids
        p.pos = copy(a)
        p.mass = copy(bu)

        @testset "Sorting" begin
            pc = sort(p, :id)
            @test issorted(pc.id)
            @test p.id == ids
            ind = searchsortedfirst(pc.id, p.id[1])
            @test pc.pos[:, ind] == p.pos[:, 1]
            @test pc.mass[ind] == p.mass[1]

            pc = sort(p, :mass; affect=(:mass, :pos), alg=RadixSort)
            @test issorted(pc.mass)#
            @test haskey(pc, :pos)
            @test !haskey(pc, :id)

            pc = sort(p, :mass; affect=(:pos,), alg=RadixSort)
            @test !haskey(pc, :mass)

            pc = deepcopy(p)
            sort!(pc, :id; alg=RadixSort)
            @test issorted(pc.id)
            ind = searchsortedfirst(pc.id, p.id[1])
            @test pc.pos[:, ind] == p.pos[:, 1]
            @test pc.mass[ind] == p.mass[1]
        end

        @testset "Filtering" begin
            massmin = 0.5u"kg"
            pc = filter(p -> p.mass .> massmin, p)
            @test all(pc.mass .> massmin)
            @test p.mass == bu
            ind = findfirst(>(massmin), p.mass)
            @test pc.pos[:, 1] == p.pos[:, ind]
            @test pc.id[1] == p.id[ind]

            @test filter(p -> findall(p.mass .> massmin), p) == pc

            pc = filter(p -> p.mass .> massmin, p; affect=(:id, :mass))
            @test all(pc.mass .> massmin)
            @test haskey(pc, :id)
            @test !haskey(pc, :pos)

            pc = deepcopy(p)
            filter!(p -> p.mass .> massmin, pc)
            @test all(pc.mass .> massmin)
            @test p.mass == bu
            ind = findfirst(>(massmin), p.mass)
            @test pc.pos[:, 1] == p.pos[:, ind]
            @test pc.id[1] == p.id[ind]


            ids_wanted = sample(1:1000, 100; replace=false)

            pc = filter(p, ids=ids_wanted)
            @test all(in.(pc.id, (ids_wanted,)))
            @test !any(in.(setdiff(ids, pc.id), (ids_wanted,)))
            ind = findfirst(in(ids_wanted), p.id)
            @test pc.pos[:, 1] == p.pos[:, ind]
            @test pc.mass[1] == p.mass[ind]

            @test filter(p; ids=Set(ids_wanted)) == pc

            pc = filter(p, ids=ids_wanted; affect=(:id, :mass))
            @test all(in.(pc.id, (ids_wanted,)))
            @test !any(in.(setdiff(ids, pc.id), (ids_wanted,)))
            @test haskey(pc, :mass)
            @test !haskey(pc, :pos)

            pc = deepcopy(p)
            filter!(pc; ids=ids_wanted)
            @test all(in.(pc.id, (ids_wanted,)))
            @test !any(in.(setdiff(ids, pc.id), (ids_wanted,)))
            ind = findfirst(in(ids_wanted), p.id)
            @test pc.pos[:, 1] == p.pos[:, ind]
            @test pc.mass[1] == p.mass[ind]

            @test filter!(deepcopy(p); ids=Set(ids_wanted)) == pc
        end
    end

    @testset "Geometry" begin
        # create square/cube/hypercube of positions between 0 and 1
        pos2 = reshape(Iterators.product(0:0.1:1, 0:0.1:1) |> collect, :)
        pos2 = reduce(hcat, getindex.(pos2, i) for i in eachindex(pos2[1]))'

        pos3 = reshape(Iterators.product(0:0.1:1, 0:0.1:1, 0:0.1:1) |> collect, :)
        pos3 = reduce(hcat, getindex.(pos3, i) for i in eachindex(pos3[1]))'

        pos4 = reshape(Iterators.product(0:0.1:1, 0:0.1:1, 0:0.1:1, 0:0.1:1) |> collect, :)
        pos4 = reduce(hcat, getindex.(pos4, i) for i in eachindex(pos4[1]))'

        @testset "Hyperrectangle" begin
            hrect = CosmoHyperrectangle(0.5rand(4), 0.5(1 .+ rand(4)))
            @test hrect isa CosmoHyperrectangle{Float64,4}

            # abstract functions do nothing
            @test isnothing(Base.@invoke CP.geometry_enclosing_corners(hrect::AbstractCosmoGeometry))
            @test isnothing(Base.@invoke CP.geometry_enclosing_center(hrect::AbstractCosmoGeometry))
            @test isnothing(
                Base.@invoke CP.mask_in(pos2::AbstractMatrix{<:Number}, hrect::AbstractCosmoGeometry)
            )

            hrect = CosmoHyperrectangle([0, 0, 0] .// 100, [15, 25, 20] .// 100)
            @test hrect isa CosmoHyperrectangle{Rational{Int},3}
            @test hrect isa CosmoCuboid{Rational{Int}}

            # construct hyperrectangles from center and side lengths
            @test CosmoHyperrectangle([75, 125, 100] .// 1000, ([15, 25, 20] .// 100)...) == hrect
            @test CosmoCuboid([75, 125, 100] .// 1000, ([15, 25, 20] .// 100)...) == hrect
            @test CosmoRectangle([75, 125] .// 1000, ([15, 25] .// 100)...) ==
                  CosmoRectangle([0, 0] .// 100, [15, 25] .// 100)

            # test constructors
            @test CosmoHyperrectangle([0, 0, 0], [0.15, 0.25, 0.20]) isa CosmoCuboid{Float64}
            @test CosmoHyperrectangle{Float64,3}([0, 0, 0.0], [0.15, 0.25, 0.20]) isa CosmoCuboid{Float64}
            @test CosmoHyperrectangle{Float64,3}([0, 0, 0], [0.15, 0.25, 1 // 10]) isa CosmoCuboid{Float64}
            @test_throws AssertionError CosmoHyperrectangle{Float64,4}([0, 0, 0.0], [0.15, 0.25, 0.20])
            @test CosmoCuboid([0, 0, 0], [0.15, 0.25, 0.20]) isa CosmoCuboid{Float64}
            @test_throws AssertionError CosmoCuboid([0, 0], [0.15, 0.25])
            @test CosmoRectangle([0, 0], [0.25, 0.20]) isa CosmoRectangle{Float64}
            @test_throws AssertionError CosmoRectangle([0, 0, 0], [0.15, 0.25, 0.20])
            @test_throws AssertionError CosmoRectangle([0, 0], [0.15, 0.25, 0.20])
            @test CosmoCuboid{Float64}([0, 0, 0.0], [0.15, 0.25, 0.20]) isa CosmoCuboid{Float64}
            @test CosmoCuboid{Float32}([0, 0, 0], [0.15, 0.25, 1 // 10]) isa CosmoCuboid{Float32}
            @test CosmoRectangle{Float64}([0, 0.0], [0.25, 0.20]) isa CosmoRectangle{Float64}
            @test CosmoRectangle{Float64}([0, 0], [0.25, 0.20]) isa CosmoRectangle{Float64}

            @test CosmoHyperrectangle([0, 0], 0.15, 0.25) isa CosmoRectangle{Float64}
            @test CosmoHyperrectangle([0, 0, 0], 0.15, 0.25, 0.20) isa CosmoCuboid{Float64}
            @test CosmoHyperrectangle([0, 0, 0, 0], 0.15, 0.25, 0.20, 0.54) isa CosmoHyperrectangle{Float64,4}
            @test_throws Exception CosmoHyperrectangle([0, 0], 0.15, 0.25, 0.54)

            r = 5 // 2
            @test CosmoHypercube([0, 0, 0, 0], r) == CosmoHyperrectangle([-r, -r, -r, -r], [r, r, r, r])
            @test CosmoCube([0, 0, 0], r) == CosmoCuboid([-r, -r, -r], [r, r, r])
            @test CosmoSquare([0, 0], r) == CosmoRectangle([-r, -r], [r, r])


            cuboid = CosmoCuboid([0, 0, 0] .// 100, [15, 25, 19] .// 100)
            @test CP.geometry_enclosing_center(cuboid) == [15, 25, 19] .// 200
            @test CP.geometry_enclosing_center(CosmoCube(zeros(Int, 3), 4 // 1)) == zeros(3)
            @test CP.geometry_enclosing_corners(cuboid) == (zeros(3), [15, 25, 19] .// 100)

            pos3in = pos3[:, CP.mask_in(pos3, cuboid)]
            @test all(0 .≤ pos3in[1, :] .≤ 0.15) &&
                  all(0 .≤ pos3in[2, :] .≤ 0.25) &&
                  all(0 .≤ pos3in[3, :] .≤ 0.19)

            rect = CosmoRectangle([0, 0] .// 100, [25, 19] .// 100)
            pos2in = pos2[:, CP.mask_in(pos2, rect)]
            @test all(0 .≤ pos2in[1, :] .≤ 0.25) && all(0 .≤ pos2in[2, :] .≤ 0.19)

            hrect = CosmoHyperrectangle([0, 0, 0, 0] .// 100, [15, 25, 19, 35] .// 100)
            pos4in = pos4[:, CP.mask_in(pos4, hrect)]
            @test all(0 .≤ pos4in[1, :] .≤ 0.15) &&
                  all(0 .≤ pos4in[2, :] .≤ 0.25) &&
                  all(0 .≤ pos4in[3, :] .≤ 0.19) &&
                  all(0 .≤ pos4in[4, :] .≤ 0.35)

            @test_throws AssertionError CP.mask_in(pos3, hrect)
        end

        @testset "Hypersphere" begin
            hsphere = CosmoHypersphere(0.5rand(4), 0.5)
            @test hsphere isa CosmoHypersphere{Float64,4}

            hsphere = CosmoHypersphere([0, 0, 0] .// 100, 15 // 100)
            @test hsphere isa CosmoHypersphere{Rational{Int},3}
            @test hsphere isa CosmoSphere{Rational{Int}}

            @test hsphere == CosmoHypersphere([0, 0, 0], 15 // 100)

            # test constructors
            @test CosmoHypersphere([0.15, 0.25, 0.20], 2) isa CosmoSphere{Float64}
            @test CosmoHypersphere{Float64,3}([0.15, 0.25, 0.20], 0.2) isa CosmoSphere{Float64}
            @test CosmoHypersphere{Float32,3}([0.15, 0.25, 0.20], 0.2) isa CosmoSphere{Float32}
            @test_throws AssertionError CosmoHypersphere{Float64,4}([0.15, 0.25, 0.20], 0.2)
            @test CosmoSphere([0.15, 0.25, 0.20], 0.2) isa CosmoSphere{Float64}
            @test_throws AssertionError CosmoSphere([0.15, 0.25], 0.3)
            @test CosmoCircle([0.25, 0.20], 0.2) isa CosmoCircle{Float64}
            @test_throws AssertionError CosmoCircle([0.15, 0.25, 0.20], 0.2)
            @test CosmoSphere{Float64}([0.15, 0.25, 0.20], 0.2) isa CosmoSphere{Float64}
            @test CosmoSphere{Float32}([0.15, 0.25, 0.20], 1) isa CosmoSphere{Float32}
            @test CosmoCircle{Float64}([0.25, 0.20], 0.2) isa CosmoCircle{Float64}

            r = 12 // 100

            center = [10, 20, 30] .// 100
            sphere = CosmoSphere(center, r)
            @test CP.geometry_enclosing_center(sphere) == center
            @test CP.geometry_enclosing_corners(sphere) == ([-2, 8, 18] .// 100, [22, 32, 42] .// 100)

            pos3in = pos3[:, CP.mask_in(pos3, sphere)]
            @test all(sum(abs2, pos3in .- center; dims=1) .≤ r^2)
            pos3notin = pos3[:, .~CP.mask_in(pos3, sphere)]
            @test all(sum(abs2, pos3notin .- center; dims=1) .> r^2)

            center = [20, 30] .// 100
            circle = CosmoCircle(center, r)
            pos2in = pos2[:, CP.mask_in(pos2, circle)]
            @test all(sum(abs2, pos2in .- center; dims=1) .≤ r^2)

            center = [10, 20, 30, 40] .// 100
            hsphere = CosmoHypersphere(center, r)
            pos4in = pos4[:, CP.mask_in(pos4, hsphere)]
            @test all(sum(abs2, pos4in .- center; dims=1) .≤ r^2)

            @test_throws AssertionError CP.mask_in(pos3, hsphere)
        end

        @testset "Cylinder" begin
            scyl = CosmoStandingCylinder(0.5rand(3), 0.5, 0.5)
            @test scyl isa CosmoStandingCylinder{Float64}

            scyl = CosmoStandingCylinder([1, 2, 3] .// 10, 2 // 10, 1 // 10)
            @test scyl isa CosmoStandingCylinder{Rational{Int}}

            @test CosmoStandingCylinder(rand(3), 1, 1 // 10) isa CosmoStandingCylinder{Float64}
            @test CosmoStandingCylinder{Float32}(rand(3), 1, 1 // 10) isa CosmoStandingCylinder{Float32}
            @test_throws AssertionError CosmoStandingCylinder(rand(4), 1, 1 // 10)

            center = [1, 2, 3] .// 10
            h = 24 // 100
            r = 15 // 100
            scyl = CosmoStandingCylinder(center, h, r)

            @test CP.geometry_enclosing_center(scyl) == center
            @test CP.geometry_enclosing_corners(scyl) == ([-5, 5, 18] .// 100, [25, 35, 42] .// 100)

            pos3in = pos3[:, CP.mask_in(pos3, scyl)]
            @test @views all(sum(abs2, pos3in[1:2, :] .- center[1:2]; dims=1) .≤ r^2)
            @test @views all(18 // 100 .≤ pos3in[3, :] .≤ 42 // 100)
            pos3notin = pos3[:, .~CP.mask_in(pos3, scyl)]
            @test @views all(
                vec(sum(abs2, pos3notin[1:2, :] .- center[1:2]; dims=1) .> r^2) .|
                (18 // 100 .> pos3notin[3, :]) .|
                (pos3notin[3, :] .> 42 // 100),
            )
            @test_throws AssertionError CP.mask_in(pos2, scyl)


            cyl = CosmoCylinder(0.5rand(3), 0.5rand(3), 0.5)
            @test cyl isa CosmoCylinder{Float64}

            cyl = CosmoCylinder([1, 2, 3] .// 10, [2, 3, 4] .// 10, 1 // 10)
            @test cyl isa CosmoCylinder{Rational{Int}}

            @test CosmoCylinder(rand(3), rand(3), 1 // 10) isa CosmoCylinder{Float64}
            @test CosmoCylinder(rand(3), [0, 0, 0], 1 // 10) isa CosmoCylinder{Float64}
            @test CosmoCylinder{Float32}(rand(3), [0, 0, 0], 1 // 10) isa CosmoCylinder{Float32}
            @test_throws AssertionError CosmoCylinder(rand(4), rand(3), 0.1)
            @test_throws AssertionError CosmoCylinder(rand(3), rand(2), 0.1)

            center = center
            startpos = [10, 20, 18] .// 100
            endpos = [10, 20, 42] .// 100
            r = r
            cyl = CosmoCylinder(startpos, endpos, r)

            @test CP.geometry_enclosing_center(cyl) == center
            @test all(CP.geometry_enclosing_corners(cyl) .≈ ([-5, 5, 18] .// 100, [25, 35, 42] .// 100))

            # right mask for standing cylinder
            @test scyl == cyl
            @test CP.mask_in(pos3, scyl) == CP.mask_in(pos3, cyl)
            @test_throws AssertionError CP.mask_in(pos2, cyl)

            @test CosmoStandingCylinder(cyl) == scyl
            @test CosmoCylinder(scyl) == cyl

            # same result switching random start and end points of the cylinder
            p1 = rand(3)
            p2 = rand(3)
            r = 0.15
            c1 = CosmoCylinder(p1, p2, r)
            c2 = CosmoCylinder(p2, p1, r)
            @test CP.geometry_enclosing_center(c1) == CP.geometry_enclosing_center(c2)
            @test CP.geometry_enclosing_corners(c1) == CP.geometry_enclosing_corners(c2)
            @test CP.mask_in(pos3, c1) == CP.mask_in(pos3, c2)

            @test_throws AssertionError CosmoStandingCylinder(c1)


            scyl = CosmoStandingCylinder([1, 2, 3] .// 10, 2 // 10, 1 // 10)
            @test scyl == CosmoStandingCylinder([1, 2, 3] .// 10, 2 // 10, 1 // 10)
            @test scyl == CosmoCylinder([1, 2, 2] .// 10, [1, 2, 4] .// 10, 1 // 10)
            @test CosmoCylinder([1, 2, 2] .// 10, [1, 2, 4] .// 10, 1 // 10) == scyl
        end


        p = Particles(:dm)
        p.pos = pos3
        p.mass = rand(size(pos3, 2))

        pu = Particles(:dm)
        pu.pos = pos3 .* u"m"
        pu.mass = rand(size(pos3, 2))

        center = rand(3)
        sph = CosmoSphere(center, 0.2)
        sphu = CosmoSphere(center * u"m", 0.2u"m")
        mask = CP.mask_in(pos3, sph)

        @test filter(p, sph) == p[mask]
        @test filter(pu, sphu) == pu[mask]

        pc = filter(p, sph; affect=(:pos,))
        @test pc.pos == p.pos[:, mask]
        @test !haskey(pc, :mass)

        pc = deepcopy(p)
        filter!(pc, sph)
        @test pc == p[mask]
    end

    @testset "Particle Collection" begin
        dm = Particles(:dm)
        dm.id = sample(1:1000, 100; replace=false)
        dm.pos = rand(3, 100)

        gas = Particles(:gas)
        gas.id = sample(1001:2000, 100; replace=false)
        gas.pos = rand(3, 100)

        @testset "ParticleCollection" begin
            pc = ParticleCollection()
            @test pc isa ParticleCollection{Particles}
            @test pc.particles == Dict{Symbol,Particles}()

            @test ParticleCollection(Particles) == pc

            pc = ParticleCollection(:dm => dm, :gas => gas)
            @test pc isa ParticleCollection{Particles}
            @test pc.particles == Dict{Symbol,Particles}(:dm => dm, :gas => gas)

            @test ParticleCollection(dm, gas) == pc

            @test redshift(pc) == 0

            pc = ParticleCollection(dm)
            pcc = copy(pc)
            @test pcc.particles[:dm] === pc.particles[:dm]
            pcc.particles[:dm] = gas
            @test pcc.particles[:dm] !== pc.particles[:dm]

            @test empty(pc) == ParticleCollection()

            pc = ParticleCollection(dm, gas)
            io = IOBuffer()
            show(io, "text/plain", pc)
            @test String(take!(io)) ==
                  "ParticleCollection\ndm: 100 Particles\n id pos\ngas: 100 Particles\n id pos"
        end

        @testset "RedshiftParticleCollection" begin
            z = 0.5

            pc = RedshiftParticleCollection(0.5)
            @test pc isa RedshiftParticleCollection{Particles}
            @test pc.particles == Dict{Symbol,Particles}()
            @test pc.z == z

            @test RedshiftParticleCollection(Particles, z) == pc

            pc = RedshiftParticleCollection(z, :dm => dm, :gas => gas)
            @test pc isa RedshiftParticleCollection{Particles}
            @test pc.particles == Dict{Symbol,Particles}(:dm => dm, :gas => gas)

            @test RedshiftParticleCollection(z, dm, gas) == pc

            @test redshift(pc) == z

            pc = RedshiftParticleCollection(z, dm)
            pcc = copy(pc)
            @test pcc.particles[:dm] === pc.particles[:dm]
            pcc.particles[:dm] = gas
            @test pcc.particles[:dm] !== pc.particles[:dm]

            @test empty(pc) == RedshiftParticleCollection(z)

            pc = RedshiftParticleCollection(z, dm, gas)
            io = IOBuffer()
            show(io, "text/plain", pc)
            @test String(take!(io)) ==
                  "ParticleCollection at z = $z\ndm: 100 Particles\n id pos\ngas: 100 Particles\n id pos"

            @test issetequal(propertynames(pc), [:z, :dm, :gas])
            @test issetequal(propertynames(pc; private=true), [:z, :particles, :dm, :gas])


            rpc = RedshiftParticleCollection(0, dm, gas)
            pc = ParticleCollection(dm, gas)
            @test rpc == pc && pc == rpc
            rpc = RedshiftParticleCollection(z, dm, gas)
            @test rpc != pc && pc != rpc
        end

        @testset "AbstractParticleCollection" begin
            pc = ParticleCollection(dm, gas)

            @test CP.get_particles(pc) === pc.particles

            @test getproperty(pc, :dm) === pc.particles[:dm]
            @test getproperty(pc, :particles) === pc.particles
            @test pc.dm === pc.particles[:dm]
            @test pc[:dm] === pc.particles[:dm]

            @test_throws ErrorException("setfield!: immutable struct of type ParticleCollection cannot be changed") pc.particles =
                Dict{Symbol,Particles}()
            pc.dm = gas
            @test pc[:dm] === gas

            pc[:dm] = dm
            @test pc.dm === dm

            @test keys(pc) == keys(pc.particles)
            @test haskey(pc, :dm)
            @test !haskey(pc, :stars)
            @test !haskey(pc, :particles)

            @test values(pc) == values(pc.particles)

            @test issetequal(propertynames(pc), [:dm, :gas])

            pcc = deepcopy(pc)
            empty!(pcc)
            @test isempty(pcc.particles)
            @test isempty(pcc)

            @test !isempty(pc)

            io = IOBuffer()
            CP.show_properties(io, "text/plain", pc)
            @test String(take!(io)) == "dm: 100 Particles\n id pos\ngas: 100 Particles\n id pos"
        end
    end
end
