# Open AtomBond and keep the window open until it is closed: with a structure
# file, or without one to start from an empty cell (New structure tool).
#
#   julia --project=. scripts/run_atombond.jl examples/structures/Si_conv.vasp
#   julia --project=. scripts/run_atombond.jl

using AtomBond

v = isempty(ARGS) ? visualize() : visualize(ARGS[1])
wait(v.screen)
