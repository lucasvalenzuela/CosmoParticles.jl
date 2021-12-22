using CosmoParticles
using Rotations
using StatsBase
using Test
using Unitful

@testset "CosmoParticles.jl" begin
    CP = CosmoParticles

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

        @test CP.particle_name(p) == "Particles"
        @test issetequal(propertynames(p), [:type, :id, :pos])
        @test issetequal(propertynames(p; private=true), [:type, :props, :id, :pos])

        io = IOBuffer()
        show(io, "text/plain", p)
        @test String(take!(io)) == "dm: 100 Particles\n id pos"

        # test AbstractParticles with basic implementation Particles
        @test CP.get_props(p) === p.props

        @test getproperty(p, :id) === p.props[:id]
        @test getproperty(p, :type) === p.type
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

        @test CP.particle_number(p) == 100
        @test CP.particle_number(Particles(:gas)) == 0
        @test CP.particle_number(Particles(:gas, Dict{Symbol,Any}(:mass => 3))) == 0

        io = IOBuffer()
        CP.show_properties(io, "text/plain", p)
        @test String(take!(io)) == "100 Particles\n id pos"
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

    end
end
