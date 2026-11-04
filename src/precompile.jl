# Typical operations on particles, run during precompilation for the types of the arguments.
function _precompile_workload(pos, vel, mass)
    dm = Particles(:dm, :id => collect(1:4), :pos => pos, :vel => vel, :mass => mass)
    gas = Particles(:gas, :id => collect(5:8), :pos => copy(pos), :vel => copy(vel), :mass => copy(mass))
    pc = ParticleCollection(dm, gas)

    center = meanpos(dm)
    filter(pc, CosmoSphere(center, oneunit(eltype(pos))))
    sort(dm, :id)
    dm[[:id, :mass]]
    translate!(pc, -center)
    rotate(pc, [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0])
    angmomtot(dm)
    colnorm(dm.pos)
    Particles(pc.all)
    Tables.columntable(dm)
    return nothing
end

# The extension for UnitfulAstro runs the operations with units.
@setup_workload begin
    @compile_workload begin
        for T in (Float32, Float64)
            _precompile_workload(zeros(T, 3, 4), zeros(T, 3, 4), ones(T, 4))
        end
    end
end
