module CosmoParticlesUnitfulAstroExt

using CosmoParticles: CosmoParticles
using PrecompileTools: @compile_workload
using Unitful: @u_str
using UnitfulAstro: UnitfulAstro

# Positions in kpc, velocities in km/s, and masses in Msun, as in cosmological simulations.
@compile_workload begin
    CosmoParticles._precompile_workload(
        zeros(3, 4) * UnitfulAstro.kpc,
        zeros(3, 4) * u"km/s",
        ones(4) * UnitfulAstro.Msun,
    )
end

end
