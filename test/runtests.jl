using CosmoParticles
using FillArrays
using LazyArrays
using LinearAlgebra
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

        copy!(pc, p)
        @test pc !== p && pc == p && !isempty(pc)

        @test CP.particle_number(p) == 100
        @test CP.particle_number(Particles(:gas)) == 0
        @test CP.particle_number(Particles(:gas, Dict{Symbol,Any}(:mass => 3))) == 0

        @test isnothing(Base.@invoke(CP.particle_name(p::AbstractParticles)))

        io = IOBuffer()
        CP.show_properties(io, "text/plain", p)
        @test String(take!(io)) == "100 Particles\n id pos"
    end

    @testset "AllParticles" begin
        dm = Particles(:dm)
        dm.id = sample(1:1000, 100; replace=false)
        dm.pos = rand(3, 100)
        dm.mass = rand()
        dm.test = rand()

        gas = Particles(:gas)
        gas.id = sample(1001:2000, 100; replace=false)
        gas.pos = rand(3, 100)
        gas.mass = rand(100)
        gas.temp = rand(100)
        gas.zs = rand(2, 100)
        gas.test = rand(2, 100)

        pcoll = ParticleCollection(dm, gas)
        p = AllParticles(pcoll)


        @test CP.particle_collection(p) === pcoll
        # @test get_props(p) === p.props # not part of the API

        if p.id[1] == dm.id[1]
            @test p.id == vcat(dm.id, gas.id)
            @test p.pos == hcat(dm.pos, gas.pos)
            @test p.mass == vcat(fill(dm.mass, length(dm.id)), gas.mass)
            @test isequal.(p.temp, vcat(fill(missing, length(dm.id)), gas.temp)) |> all
            @test isequal.(p.zs, hcat(fill(missing, 2, length(dm.id)), gas.zs)) |> all
            @test p.test == hcat(fill(dm.test, 2, length(dm.id)), gas.test)
        else
            @test p.id == vcat(gas.id, dm.id)
            @test p.pos == hcat(gas.pos, dm.pos)
            @test p.mass == vcat(gas.mass, fill(dm.mass, length(dm.id)))
            @test isequal.(p.temp, vcat(gas.temp, fill(missing, length(dm.id)))) |> all
            @test isequal.(p.zs, hcat(gas.zs, fill(missing, 2, length(dm.id)))) |> all
            @test p.test == hcat(gas.test, fill(dm.test, 2, length(dm.id)))
        end

        @test_throws ErrorException p.id = rand(100)
        @test_throws KeyError p.vel

        @test issetequal(keys(p), [:id, :mass, :pos, :temp, :zs, :test])
        @test keys(p) == keys(AllParticles(pcoll))

        @test haskey(p, :id) && haskey(AllParticles(pcoll), :id)
        @test !haskey(p, :vel) && !haskey(AllParticles(pcoll), :vel)

        @test all(val === p[key] for (key, val) in zip(keys(p), values(p)))
        @test all(val === p[key] for (key, val) in pairs(p))

        @test_throws ErrorException copy(p)
        @test_throws ErrorException copy!(p, AllParticles(pcoll))
        @test_throws ErrorException copy!(p, dm)

        pc = Particles(:all)
        copy!(pc, p)
        collect(values(p)) # materializes properties in Dict
        @test pc.id == p.id
        @test pc.mass == p.mass
        @test pc.pos == p.pos
        @test all(isequal.(pc.temp, p.temp))
        @test all(isequal.(pc.zs, p.zs))
        @test pc.test == p.test

        copy!(pc, p, (:id, :mass, :temp))
        @test haskey.((pc,), [:id, :mass, :temp]) |> all
        @test !haskey(pc, :pos)
        @test isa.([pc.id, pc.mass, pc.temp], Array) |> all

        @test Particles(p, (:id, :mass, :temp)) == pc

        @test_throws ErrorException empty(p)
        @test_throws ErrorException empty!(p)

        @test !isempty(p)
        @test isempty(ParticleCollection() |> AllParticles)
        @test isempty(ParticleCollection(Particles(:dm), Particles(:gas)) |> AllParticles)

        @test p == AllParticles(pcoll)

        @test CP.particle_name(p) == "Particles"
        @test CP.particle_number(p) == 200

        list = [:id, :mass, :pos, :temp, :zs, :test]
        @test issetequal(propertynames(p), list)
        @test issetequal(propertynames(p; private=true), [fieldnames(AllParticles) |> collect; list])

        io = IOBuffer()
        show(io, "text/plain", p)
        @test String(take!(io)) == "all: 200 Particles\n id mass pos temp test zs"


        mask = isodd.(p.id)
        ind = findall(mask)
        @test_throws ErrorException p[mask]
        @test_throws ErrorException p[ind]


        # cover case where the dimension is not found from a vector or matrix
        # and where only one particle subtype exists
        pos = rand(3, 100)
        mass = rand()
        pc = ParticleCollection(Particles(:dm, :pos => pos, :mass => mass))
        ap = pc.all
        @test ap.pos == pos
        @test ap.mass == fill(mass, 100)
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

            @test issetequal(propertynames(pc), [:z, :all, :dm, :gas])
            @test issetequal(propertynames(pc; private=true), [:z, :all, :particles, :dm, :gas])


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

            @test_throws ErrorException(
                "setfield!: immutable struct of type ParticleCollection cannot be changed",
            ) pc.particles = Dict{Symbol,Particles}()
            pc.dm = gas
            @test pc[:dm] === gas

            pc[:dm] = dm
            @test pc.dm === dm

            @test pc.all == AllParticles(pc)
            @test_throws ErrorException pc.all = dm
            @test_throws KeyError pc[:all]

            @test keys(pc) == keys(pc.particles)
            @test haskey(pc, :dm)
            @test !haskey(pc, :stars)
            @test !haskey(pc, :particles)

            @test values(pc) == values(pc.particles)

            @test issetequal(propertynames(pc), [:all, :dm, :gas])

            pcc = deepcopy(pc)
            empty!(pcc)
            @test isempty(pcc.particles)
            @test isempty(pcc)

            @test !isempty(pc)

            copy!(pcc, pc)
            @test pcc !== pc && pcc == pc && !isempty(pcc)

            io = IOBuffer()
            CP.show_properties(io, "text/plain", pc)
            @test String(take!(io)) == "dm: 100 Particles\n id pos\ngas: 100 Particles\n id pos"
        end
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

            a = rand(45)
            b = rand(55)
            c = ApplyArray(vcat, a, b)
            @test CP._applyind(c, mask) == Array(c)[mask]
            @test CP._applyind(c, ind) == Array(c)[ind]

            a = rand(3, 45)
            b = rand(3, 55)
            c = ApplyArray(hcat, a, b)
            @test CP._applyind(c, mask) == Array(c)[:, mask]
            @test CP._applyind(c, ind) == Array(c)[:, ind]

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

            pc = CP.applyind(p, mask; affect=(:pos, :mass, :foo))
            @test !haskey(pc, :id)
            @test !haskey(pc, :foo)
            @test pc.pos == p.pos[:, mask]
            @test pc.mass == p.mass


            pc = ParticleCollection(p)
            ap = pc.all
            @test_throws ErrorException CP.applyind(ap, ind)
            @test_throws ErrorException CP.applyind!(ap, ind)
            @test_throws ErrorException CP.applyind(ap, mask)
            @test_throws ErrorException CP.applyind!(ap, mask)
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

        @testset "Unitful operations" begin
            a = rand(100)
            au = a * u"m"
            af = convert(Vector{Float32}, a)
            afu = af * u"m"
            b = rand()

            @test CP.product_preserve_type(a, b) == a .* b
            @test CP.product_preserve_type(au, b) == au .* b
            @test CP.product_preserve_type(af, b) == af .* convert(Float32, b)
            @test CP.product_preserve_type(afu, b) == afu .* convert(Float32, b)

            ac = copy(a)
            CP.product_preserve_type!(ac, b)
            @test ac == a .* b

            auc = copy(au)
            CP.product_preserve_type!(auc, b)
            @test auc == au .* b

            afc = copy(af)
            CP.product_preserve_type!(afc, b)
            @test afc == af .* convert(Float32, b)

            afuc = copy(afu)
            CP.product_preserve_type!(afuc, b)
            @test afuc == afu .* convert(Float32, b)
        end

        @testset "Lazy ustrip" begin
            @test CP.ustrip_lazy(3) == 3

            a = rand(4)
            b = rand(3)
            au = a * u"m"
            bu = b * u"m"
            @test CP.ustrip_lazy(a) === a
            @test CP.ustrip_lazy(au) == a
            @test CP.ustrip_lazy(Diagonal(au)) == Diagonal(a)
            @test CP.ustrip_lazy(Bidiagonal(au, bu, :U)) == Bidiagonal(a, b, :U)
            @test CP.ustrip_lazy(Tridiagonal(bu, au, bu)) == Tridiagonal(b, a, b)
            @test CP.ustrip_lazy(SymTridiagonal(au, bu)) == SymTridiagonal(a, b)

            auview = @view au[1:3]
            @test CP.ustrip_lazy(auview) == a[1:3]

            c = vcat(a, b)
            cu = ApplyArray(vcat, au, bu)
            @test CP.ustrip_lazy(cu) == c

            c = fill(2.5, 4)
            cu = Fill(2.5u"m", 4)
            @test CP.ustrip_lazy(cu) == c
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


            @testset "Particle Rotation" begin
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


            @testset "Particle Collection Rotation" begin
                dm = Particles(:dm)
                dm.pos = copy(au)
                dm.vel = copy(a)
                gas = Particles(:gas)
                gas.pos = copy(au)
                gas.vel = copy(a)

                pc = ParticleCollection(dm, gas)

                pcc = rotate(pc, rotmat, :pos)
                @test pc.dm.pos == pc.gas.pos == au
                @test pcc.dm.pos == pcc.gas.pos == aurot
                @test pc.dm.vel === pcc.dm.vel == a

                pcc = rotate(pc, rotmat)
                @test pc.dm.pos == pc.gas.pos == au
                @test pcc.dm.pos == pcc.gas.pos == aurot
                @test pc.dm.vel == pc.gas.vel == a
                @test pcc.dm.vel == pcc.gas.vel == arot


                pcc = deepcopy(pc)
                rotate!(pcc, rotmat, :vel)
                @test pcc.dm.pos == pcc.gas.pos == au
                @test pcc.dm.vel == pcc.gas.vel == arot

                pcc = deepcopy(pc)
                rotate!(pcc, rotmat)
                @test pcc.dm.pos == pcc.gas.pos == aurot
                @test pcc.dm.vel == pcc.gas.vel == arot


                pccc = deepcopy(pc)
                ap = pccc.all
                @test_throws ErrorException rotate(ap, rotmat)

                rotate!(ap, rotmat)
                @test pccc == pcc
            end
        end

        @testset "Translations" begin
            da = rand(3)
            dau = rand(3) * u"m"
            a = rand(3, 100)
            au = a * u"m"

            atrans = a .+ da
            autrans = au .+ dau

            @testset "Particle Translations" begin
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

            @testset "Particle Collection Translations" begin
                dm = Particles(:dm)
                dm.pos = copy(au)
                dm.vel = copy(a)
                gas = Particles(:gas)
                gas.pos = copy(au)
                gas.vel = copy(a)

                pc = ParticleCollection(dm, gas)

                pcc = translate(pc, dau)
                @test pc.dm.pos == pc.gas.pos == au
                @test pcc.dm.pos == pcc.gas.pos == autrans
                @test pc.dm.vel === pcc.dm.vel == a

                pcc = translate(pc, da, :vel)
                @test pc.dm.pos === pcc.dm.pos == au
                @test pc.dm.vel == pc.gas.vel == a
                @test pcc.dm.vel == pcc.gas.vel == atrans

                pcc = deepcopy(pc)
                translate!(pcc, dau)
                @test pcc.dm.pos == pcc.gas.pos == autrans
                @test pcc.dm.vel == pcc.gas.vel == a

                pcc = deepcopy(pc)
                translate!(pcc, da, :vel)
                @test pcc.dm.pos == pcc.gas.pos == au
                @test pcc.dm.vel == pcc.dm.vel == atrans


                pccc = deepcopy(pc)
                ap = pccc.all
                @test_throws ErrorException translate(ap, da)

                translate!(ap, dau)
                @test pccc == translate(pc, dau)

                pccc = deepcopy(pc)
                ap = pccc.all
                translate!(ap, da, :vel)
                @test pccc == pcc
            end
        end


        @testset "Comoving" begin
            z = 1.5
            a = rand(Float32, 3, 100)
            au = a * u"m"
            n = rand(Float32, 100)
            nu = n * u"m^-3"
            m = copy(n)

            propexp = ((:pos, 1), (:vel, 1), (:n, -3))

            acom = a .* (1 + z)
            aucom = au .* (1 + z)
            ncom = n .* (1 + z)^-3
            nucom = nu .* (1 + z)^-3

            aphys = a ./ (1 + z)
            auphys = au ./ (1 + z)
            nphys = n ./ (1 + z)^-3
            nuphys = nu ./ (1 + z)^-3

            @testset "Particle Comoving" begin
                p = Particles(:dm)
                p.pos = a
                p.n = n
                p.mass = m

                pc = to_comoving(p, z; propexp)
                @test pc.pos ≈ acom
                @test pc.n ≈ ncom
                @test pc.mass ≈ m

                pcc = deepcopy(p)
                to_comoving!(pcc, z; propexp)
                @test pcc == pc

                pc = to_physical(p, z; propexp)
                @test pc.pos ≈ aphys
                @test pc.n ≈ nphys
                @test pc.mass ≈ m

                pcc = deepcopy(p)
                to_physical!(pcc, z; propexp)
                @test pcc == pc

                p = Particles(:dm)
                p.pos = au
                p.n = nu
                p.mass = m

                pc = to_comoving(p, z; propexp)
                @test pc.pos ≈ aucom
                @test pc.n ≈ nucom
                @test pc.mass ≈ m

                pcc = deepcopy(p)
                to_comoving!(pcc, z; propexp)
                @test pcc == pc

                pc = to_physical(p, z; propexp)
                @test pc.pos ≈ auphys
                @test pc.n ≈ nuphys
                @test pc.mass ≈ m

                pcc = deepcopy(p)
                to_physical!(pcc, z; propexp)
                @test pcc == pc
            end

            @testset "Particle Collection Comoving" begin
                dm = Particles(:dm)
                dm.pos = copy(au)
                dm.n = copy(n)
                dm.mass = copy(m)
                gas = Particles(:gas)
                gas.pos = copy(a)
                gas.n = copy(nu)
                gas.mass = copy(m)

                pc = RedshiftParticleCollection(z, dm, gas)

                pcc = to_comoving(pc; propexp)
                @test pcc.dm == to_comoving(dm, z; propexp)
                @test pcc.gas == to_comoving(gas, z; propexp)

                pccc = deepcopy(pc)
                to_comoving!(pccc; propexp)
                @test pccc == pcc

                pcc = to_physical(pc; propexp)
                @test pcc.dm == to_physical(dm, z; propexp)
                @test pcc.gas == to_physical(gas, z; propexp)

                pccc = deepcopy(pc)
                to_physical!(pccc; propexp)
                @test pccc == pcc


                pccc = deepcopy(pc)
                ap = pccc.all
                @test_throws ErrorException to_comoving(ap, z; propexp)
                @test_throws ErrorException to_physical(ap, z; propexp)

                to_comoving!(ap, z; propexp)
                @test pccc == to_comoving(pc; propexp)

                pccc = deepcopy(pc)
                ap = pccc.all
                to_physical!(ap, z; propexp)
                @test pccc == to_physical(pc; propexp)
            end
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


        @testset "Particle Sorting" begin
            pc = sort(p, :id)
            @test issorted(pc.id)
            @test p.id == ids
            ind = searchsortedfirst(pc.id, p.id[1])
            @test pc.pos[:, ind] == p.pos[:, 1]
            @test pc.mass[ind] == p.mass[1]

            pc = sort(p, :mass; affect=(:mass, :pos), alg=RadixSort)
            @test issorted(pc.mass)
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

        @testset "Particle Collection Sorting" begin
            dm = p
            gas = Particles(:gas, p.props)
            pc = ParticleCollection(dm, gas)

            pcc = sort(pc, :id)
            @test pc != pcc
            @test pcc.dm == sort(dm, :id)
            @test pcc.gas == sort(gas, :id)

            affect = (:mass, :pos)
            pcc = sort(pc, :mass; affect, alg=RadixSort)
            @test pcc.dm == sort(dm, :mass; affect)
            @test pcc.gas == sort(gas, :mass; affect)

            pcc = deepcopy(pc)
            sort!(pcc, :id; alg=RadixSort)
            @test pcc == sort(pc, :id)


            ap = pc.all
            @test_throws ErrorException sort(ap, :id; affect)
            @test_throws ErrorException sort!(ap, :id; alg=RadixSort)
        end

        @testset "Particle Filtering" begin
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

            pc = filter(p; ids=ids_wanted)
            @test all(in.(pc.id, (ids_wanted,)))
            @test !any(in.(setdiff(ids, pc.id), (ids_wanted,)))
            ind = findfirst(in(ids_wanted), p.id)
            @test pc.pos[:, 1] == p.pos[:, ind]
            @test pc.mass[1] == p.mass[ind]

            @test filter(p; ids=Set(ids_wanted)) == pc

            pc = filter(p; ids=ids_wanted, affect=(:id, :mass))
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

        @testset "Particle Collection Sorting" begin
            dm = p
            gas = Particles(:gas, p.props)
            pc = ParticleCollection(dm, gas)

            massmin = 0.5u"kg"
            f = p -> p.mass .> massmin
            pcc = filter(f, pc)
            @test pcc.dm == filter(f, dm)

            @test filter(p -> findall(p.mass .> massmin), pc) == pcc

            affect = (:id, :mass)
            pcc = filter(f, pc; affect)
            @test pcc.dm == filter(f, dm; affect)

            pcc = deepcopy(pc)
            filter!(f, pcc)
            @test pcc == filter(f, pc)

            ids_wanted = sample(1:1000, 100; replace=false)

            pcc = filter(pc; ids=ids_wanted)
            @test pcc.dm == filter(dm; ids=ids_wanted)
            @test filter(pc; ids=Set(ids_wanted)) == pcc

            affect = (:id, :mass)
            pcc = filter(pc; ids=ids_wanted, affect)
            @test pcc.dm == filter(dm; ids=ids_wanted, affect)

            pcc = deepcopy(pc)
            filter!(pcc; ids=ids_wanted)
            @test pcc == filter(pcc; ids=ids_wanted)

            @test filter!(deepcopy(pc); ids=Set(ids_wanted)) == pcc


            ap = pc.all
            @test_throws ErrorException filter(f, ap; affect)
            @test_throws ErrorException filter!(f, ap)
            @test_throws ErrorException filter(ap; ids=ids_wanted)
            @test_throws ErrorException filter!(ap; ids=ids_wanted)
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


        @testset "Filtering" begin
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


            # Particle collection
            dm = p
            gas = Particles(:gas, p.props)
            pc = ParticleCollection(dm, gas)

            pcc = filter(pc, sph)
            @test pcc.dm == dm[mask]

            pcc = filter(pc, sph; affect=(:pos,))
            @test pcc.dm == filter(dm, sph; affect=(:pos,))

            pcc = deepcopy(pc)
            filter!(pcc, sph)
            @test pcc == filter(pc, sph)


            ap = pc.all
            @test_throws ErrorException filter(ap, sph)
            @test_throws ErrorException filter!(ap, sph)
        end
    end


    @testset "Properties" begin
        @testset "Norm" begin
            a2 = rand(2, 100)
            a3 = rand(3, 100)
            a4 = rand(4, 100)

            n2 = norm.(eachcol(a2))
            n3 = norm.(eachcol(a3))
            n4 = norm.(eachcol(a4))

            @test colnorm(a2) ≈ n2
            @test colnorm(a3) ≈ n3
            @test colnorm(a4) ≈ n4

            @test colnorm(a2 * u"m") ≈ n2 * u"m"
            @test colnorm(a3 * u"m") ≈ n3 * u"m"
            @test colnorm(a4 * u"m") ≈ n4 * u"m"

            @test colnorm2(a2) ≈ n2 .^ 2
            @test colnorm2(a3) ≈ n3 .^ 2
            @test colnorm2(a4) ≈ n4 .^ 2

            @test colnorm2(a2 * u"m") ≈ (n2 * u"m") .^ 2
            @test colnorm2(a3 * u"m") ≈ (n3 * u"m") .^ 2
            @test colnorm2(a4 * u"m") ≈ (n4 * u"m") .^ 2

            am2 = similar(a2, Union{Float64,Missing})
            am3 = similar(a3, Union{Float64,Missing})
            am4 = similar(a4, Union{Float64,Missing})
            am2 .= a2
            am3 .= a3
            am4 .= a4
            am2[:, 1] .= missing
            am3[:, 1] .= missing
            am4[:, 1] .= missing

            @test all(colnorm(am2) .≈ n2) |> ismissing
            @test all(colnorm(am3) .≈ n3) |> ismissing
            @test all(colnorm(am4) .≈ n4) |> ismissing

            @test all(colnorm(am2 * u"m") .≈ n2 * u"m") |> ismissing
            @test all(colnorm(am3 * u"m") .≈ n3 * u"m") |> ismissing
            @test all(colnorm(am4 * u"m") .≈ n4 * u"m") |> ismissing

            @test all(colnorm2(am2) .≈ n2 .^ 2) |> ismissing
            @test all(colnorm2(am3) .≈ n3 .^ 2) |> ismissing
            @test all(colnorm2(am4) .≈ n4 .^ 2) |> ismissing

            @test all(colnorm2(am2 * u"m") .≈ (n2 * u"m") .^ 2) |> ismissing
            @test all(colnorm2(am3 * u"m") .≈ (n3 * u"m") .^ 2) |> ismissing
            @test all(colnorm2(am4 * u"m") .≈ (n4 * u"m") .^ 2) |> ismissing
        end

        @testset "Norm around Origin" begin
            a2 = rand(2, 100)
            a3 = rand(3, 100)
            a4 = rand(4, 100)

            o2 = rand(2)
            o3 = rand(3)
            o4 = rand(4)

            n2 = norm.(eachcol(a2 .- o2))
            n3 = norm.(eachcol(a3 .- o3))
            n4 = norm.(eachcol(a4 .- o4))

            @test colnorm(a2, o2) ≈ n2
            @test colnorm(a3, o3) ≈ n3
            @test colnorm(a4, o4) ≈ n4

            @test colnorm(a2 * u"m", o2 * u"m") ≈ n2 * u"m"
            @test colnorm(a3 * u"m", o3 * u"m") ≈ n3 * u"m"
            @test colnorm(a4 * u"m", o4 * u"m") ≈ n4 * u"m"

            @test colnorm2(a2, o2) ≈ n2 .^ 2
            @test colnorm2(a3, o3) ≈ n3 .^ 2
            @test colnorm2(a4, o4) ≈ n4 .^ 2

            @test colnorm2(a2 * u"m", o2 * u"m") ≈ (n2 * u"m") .^ 2
            @test colnorm2(a3 * u"m", o3 * u"m") ≈ (n3 * u"m") .^ 2
            @test colnorm2(a4 * u"m", o4 * u"m") ≈ (n4 * u"m") .^ 2

            am2 = similar(a2, Union{Float64,Missing})
            am3 = similar(a3, Union{Float64,Missing})
            am4 = similar(a4, Union{Float64,Missing})
            am2 .= a2
            am3 .= a3
            am4 .= a4
            am2[:, 1] .= missing
            am3[:, 1] .= missing
            am4[:, 1] .= missing

            @test all(colnorm(am2, o2) .≈ n2) |> ismissing
            @test all(colnorm(am3, o3) .≈ n3) |> ismissing
            @test all(colnorm(am4, o4) .≈ n4) |> ismissing

            @test all(colnorm(am2 * u"m", o2 * u"m") .≈ n2 * u"m") |> ismissing
            @test all(colnorm(am3 * u"m", o3 * u"m") .≈ n3 * u"m") |> ismissing
            @test all(colnorm(am4 * u"m", o4 * u"m") .≈ n4 * u"m") |> ismissing

            @test all(colnorm2(am2, o2) .≈ n2 .^ 2) |> ismissing
            @test all(colnorm2(am3, o3) .≈ n3 .^ 2) |> ismissing
            @test all(colnorm2(am4, o4) .≈ n4 .^ 2) |> ismissing

            @test all(colnorm2(am2 * u"m", o2 * u"m") .≈ (n2 * u"m") .^ 2) |> ismissing
            @test all(colnorm2(am3 * u"m", o3 * u"m") .≈ (n3 * u"m") .^ 2) |> ismissing
            @test all(colnorm2(am4 * u"m", o4 * u"m") .≈ (n4 * u"m") .^ 2) |> ismissing
        end

        @testset "Dot Product" begin
            a2 = rand(2, 100)
            a3 = rand(3, 100)
            a4 = rand(4, 100)
            b2 = rand(2, 100)
            b3 = rand(3, 100)
            b4 = rand(4, 100)

            c2 = reduce(hcat, [dot(ai, bi) for (ai, bi) in zip(eachcol(a2), eachcol(b2))]) |> vec
            c3 = reduce(hcat, [dot(ai, bi) for (ai, bi) in zip(eachcol(a3), eachcol(b3))]) |> vec
            c4 = reduce(hcat, [dot(ai, bi) for (ai, bi) in zip(eachcol(a4), eachcol(b4))]) |> vec

            for (a, b, c) in zip([a2, a3, a4], [b2, b3, b4], [c2, c3, c4])
                @test coldot(a, b) ≈ c
                @test coldot(a * u"m", b) ≈ c * u"m"
                @test coldot(a, b * u"m") ≈ c * u"m"
                @test coldot(a * u"m", b * u"m") ≈ c * u"m^2"

                am = similar(a, Union{Float64,Missing})
                am .= a
                am[:, 1] .= missing
                bm = similar(b, Union{Float64,Missing})
                bm .= b
                bm[:, 1] .= missing

                @test all(coldot(am, b) .≈ c) |> ismissing
                @test all(coldot(a * u"m", bm) .≈ c * u"m") |> ismissing
                @test all(coldot(a, bm * u"m") .≈ c * u"m") |> ismissing
                @test all(coldot(am * u"m", bm) .≈ c * u"m") |> ismissing
                @test all(coldot(am * u"m", bm * u"m") .≈ c * u"m^2") |> ismissing
            end

            am = similar(a3, Union{Float64,Nothing})
            am .= a3
            @test coldot(am, b3) ≈ c3
        end

        @testset "Cross Product" begin
            a = rand(3, 100)
            b = rand(3, 100)
            au = a * u"m"
            bu = b * u"m"

            c = reduce(hcat, [cross(ai, bi) for (ai, bi) in zip(eachcol(a), eachcol(b))])

            @test colcross(a, b) ≈ c
            @test colcross(au, b) ≈ c * u"m"
            @test colcross(a, bu) ≈ c * u"m"
            @test colcross(au, bu) ≈ c * u"m^2"

            am = similar(a, Union{Float64,Missing})
            am .= a
            am[:, 1] .= missing
            bm = similar(b, Union{Float64,Missing})
            bm .= b
            bm[:, 1] .= missing

            @test all(colcross(am, b) .≈ c) |> ismissing
            @test all(colcross(au, bm) .≈ c * u"m") |> ismissing
            @test all(colcross(a, bm * u"m") .≈ c * u"m") |> ismissing
            @test all(colcross(am * u"m", bm) .≈ c * u"m") |> ismissing
            @test all(colcross(am * u"m", bm * u"m") .≈ c * u"m^2") |> ismissing

            am = similar(a, Union{Float64,Nothing})
            am .= a
            @test colcross(am, b) ≈ c
        end

        @testset "Angular momentum" begin
            a = rand(3, 100)
            b = rand(3, 100)
            au = a * u"m"
            bu = b * u"m"

            for m in [rand(), rand(100)]
                mu = m * u"kg"

                j = m' .* colcross(a, b)

                @test angmom(a, b, m) ≈ j
                @test angmom(au, b, m) ≈ j * u"m"
                @test angmom(a, bu, m) ≈ j * u"m"
                @test angmom(au, bu, m) ≈ j * u"m^2"
                @test angmom(a, b, mu) ≈ j * u"kg"
                @test angmom(au, b, mu) ≈ j * u"m*kg"
                @test angmom(a, bu, mu) ≈ j * u"m*kg"
                @test angmom(au, bu, mu) ≈ j * u"m^2*kg"

                for func in [angmomtot, angmomtot_stable]
                    @test func(a, b, m) ≈ sum(j; dims=2)
                    @test func(au, b, m) ≈ sum(j * u"m"; dims=2)
                    @test func(a, bu, m) ≈ sum(j * u"m"; dims=2)
                    @test func(au, bu, m) ≈ sum(j * u"m^2"; dims=2)
                    @test func(a, b, mu) ≈ sum(j * u"kg"; dims=2)
                    @test func(au, b, mu) ≈ sum(j * u"m*kg"; dims=2)
                    @test func(a, bu, mu) ≈ sum(j * u"m*kg"; dims=2)
                    @test func(au, bu, mu) ≈ sum(j * u"m^2*kg"; dims=2)
                end


                origin = rand(3)
                velorigin = rand(3)
                jorigin = m' .* colcross(a .- origin, b)
                @test angmom(au, bu, mu; origin=origin * u"m") ≈ jorigin * u"m^2*kg"
                @test angmomtot(au, bu, mu; origin=origin * u"m") ≈ sum(jorigin * u"m^2*kg"; dims=2)
                jorigin = m' .* colcross(a, b .- velorigin)
                @test angmom(au, bu, mu; velorigin=velorigin * u"m") ≈ jorigin * u"m^2*kg"
                @test angmomtot(au, bu, mu; velorigin=velorigin * u"m") ≈ sum(jorigin * u"m^2*kg"; dims=2)
                @test angmomtot_stable(au, bu, mu; velorigin=velorigin * u"m") ≈
                      sum(jorigin * u"m^2*kg"; dims=2)
                jorigin = m' .* colcross(a .- origin, b .- velorigin)
                @test angmom(au, bu, mu; origin=origin * u"m", velorigin=velorigin * u"m") ≈
                      jorigin * u"m^2*kg"
                @test angmomtot(au, bu, mu; origin=origin * u"m", velorigin=velorigin * u"m") ≈
                      sum(jorigin * u"m^2*kg"; dims=2)
                @test angmomtot_stable(au, bu, mu; origin=origin * u"m", velorigin=velorigin * u"m") ≈
                      sum(jorigin * u"m^2*kg"; dims=2)

                am = similar(a, Union{Float64,Missing})
                am .= a
                am[:, 1] .= missing
                bm = similar(b, Union{Float64,Missing})
                bm .= b
                bm[:, 1] .= missing

                @test all(angmom(am, b, m) .≈ j) |> ismissing
                @test all(angmom(au, bm, m) .≈ j * u"m") |> ismissing
                @test all(angmom(a, bm * u"m", m) .≈ j * u"m") |> ismissing
                @test all(angmom(am * u"m", bm, m) .≈ j * u"m") |> ismissing
                @test all(angmom(am * u"m", bm * u"m", m) .≈ j * u"m^2") |> ismissing
                @test all(angmom(am, b, mu) .≈ j * u"kg") |> ismissing
                @test all(angmom(au, bm, mu) .≈ j * u"m*kg") |> ismissing
                @test all(angmom(a, bm * u"m", mu) .≈ j * u"m*kg") |> ismissing
                @test all(angmom(am * u"m", bm, mu) .≈ j * u"m*kg") |> ismissing
                @test all(angmom(am * u"m", bm * u"m", mu) .≈ j * u"m^2*kg") |> ismissing

                am = similar(a, Union{Float64,Nothing})
                am .= a
                @test angmom(am, b, m) ≈ j
            end
        end


        @testset "Particles properties" begin
            dm = Particles(:dm, :pos => rand(3, 100), :vel => rand(Float32, 3, 100), :mass => 2)
            gas = Particles(
                :gas,
                :pos => rand(3, 100),
                :vel => rand(Float32, 3, 100),
                :mass => rand(100),
                :temp => rand(100),
                :test => 1,
                :mass2 => Fill(2, 100),
            )

            dmu = deepcopy(dm)
            gasu = deepcopy(gas)
            dmu.pos *= u"m"
            dmu.vel *= u"km/s"
            dmu.mass *= u"kg"
            gasu.pos *= u"m"
            gasu.vel *= u"km/s"
            gasu.mass *= u"kg"
            gasu.mass2 *= u"kg"
            gasu.temp *= u"K"

            for (dm, gas) in zip([dm, dmu], [gas, gasu])
                cp = ParticleCollection(dm, gas)
                ap = cp.all
                p = Particles(ap)

                @test meanpos(dm) == mean(dm.pos; dims=2)
                @test meanpos(dm; massweighted=false) == mean(dm.pos; dims=2)
                @test meanpos(gas) ≈ sum(gas.pos .* gas.mass'; dims=2) / sum(gas.mass)
                @test meanpos(gas; massweighted=false) == mean(gas.pos; dims=2)
                @test meanpos(ap) ≈ meanpos(p)

                @test meanprop(dm, :pos) == mean(dm.pos; dims=2)
                @test meanprop(gas, :pos; massprop=:mass2) == mean(gas.pos; dims=2)
                @test meanprop(gas, :temp) ≈ sum(gas.temp .* gas.mass) / sum(gas.mass)
                @test meanprop(gas, :temp; massweighted=false) == mean(gas.temp)
                @test meanprop(gas, :temp; massprop=:mass2) == mean(gas.temp)
                @test meanprop(dm, :mass) == dm.mass
                @test meanprop(dm, :mass; massweighted=false) == dm.mass
                @test meanprop(gas, :test) == gas.test
                @test meanprop(gas, :mass2) == gas.mass2[1]
                @test meanprop(gas, :mass2; massweighted=false) == gas.mass2[1]
                @test meanprop(gas, :mass2; massprop=:mass2) == gas.mass2[1]
                @test meanprop(ap, :mass) ≈ meanvel(p, :mass)

                @test meanvel(dm) == mean(dm.vel; dims=2)
                @test meanvel(dm; massweighted=false) == mean(dm.vel; dims=2)
                @test meanvel(gas) ≈ sum(gas.vel .* gas.mass'; dims=2) / sum(gas.mass)
                @test meanvel(gas; massweighted=false) == mean(gas.vel; dims=2)
                @test meanvel(ap) ≈ meanvel(p)

                @test sumprop(dm, :pos) == sum(dm.pos; dims=2)
                @test sumprop(gas, :pos) == sum(gas.pos; dims=2)
                @test sumprop(dm, :mass) == CP.particle_number(dm) * dm.mass
                @test sumprop(gas, :mass2) == sum(gas.mass2)
                @test sumprop(ap, :pos) == sumprop(p, :pos)

                @test angmom(dm; angmomprop=:j) == angmom(dm.pos, dm.vel, dm.mass)
                @test angmomtot(dm; angmomprop=:j) == angmomtot(dm.pos, dm.vel, dm.mass)
                @test angmomtot_stable(dm; angmomprop=:j) == angmomtot_stable(dm.pos, dm.vel, dm.mass)
                dm.j = angmom(dm)
                @test angmom(dm; angmomprop=:j) ≈ angmom(dm.pos, dm.vel, dm.mass)
                @test angmomtot(dm; angmomprop=:j) ≈ angmomtot(dm.pos, dm.vel, dm.mass)
                @test angmomtot_stable(dm; angmomprop=:j) ≈ angmomtot_stable(dm.pos, dm.vel, dm.mass)
                @test angmom(gas; angmomprop=:j) == angmom(gas.pos, gas.vel, gas.mass)
                @test angmomtot(gas; angmomprop=:j) == angmomtot(gas.pos, gas.vel, gas.mass)
                @test angmomtot_stable(gas; angmomprop=:j) == angmomtot_stable(gas.pos, gas.vel, gas.mass)

                origin = 2 * dm.pos[:, 1]
                velorigin = 2 * dm.vel[:, 1]

                for func in [angmom, angmomtot, angmomtot_stable]
                    @test func(gas; origin) == func(gas.pos, gas.vel, gas.mass; origin)
                    @test func(gas; velorigin) == func(gas.pos, gas.vel, gas.mass; velorigin)
                    @test func(gas; origin, velorigin) == func(gas.pos, gas.vel, gas.mass; origin, velorigin)
                end
            end
        end
    end
end
