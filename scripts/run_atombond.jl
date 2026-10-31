# Open a structure in an AtomBond window and keep it open until closed.
#
#   julia --project=. scripts/run_atombond.jl examples/structures/Si_conv.vasp

using AtomBond

path = isempty(ARGS) ? joinpath(@__DIR__, "..", "examples", "structures", "Si_conv.vasp") : ARGS[1]
v = visualize(path)
wait(v.screen)
