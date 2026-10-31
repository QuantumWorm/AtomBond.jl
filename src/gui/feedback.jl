# Click feedback, so a click visibly (and audibly) registers:
#
# - buttons flash blue for a moment after a click (the colour Makie shows
#   while the mouse is held down is often never rendered on a quick click);
# - the chosen item of a popup menu flashes before the menu closes;
# - a click in the 3D view leaves a short expanding ring at the cursor:
#   deep orange on an atom, grey on empty space, blue for a right-click;
# - a short click sound plays through `paplay` when it is installed
#   (WSLg has a PulseAudio server; `sudo apt install pulseaudio-utils`).
#
# click_feedback!(flash=false, sound=false) turns them off.
#
# Effects are timed by drawn frames (the figure's render ticks), not by the
# clock alone: the first click on something can spend a while compiling, and
# a flash timed from the click would be over before it was ever drawn.

const CLICK_FLASH = Ref(true)
const CLICK_SOUND = Ref(true)
const FLASH_TIME = 0.15                          # s
const FLASH_COLOR = RGBf(0.26, 0.52, 0.96)
const RIPPLE_TIME = 0.35                         # s
const RIPPLE_COLORS = Dict(:atom => RGBf(0.85, 0.30, 0.05), :empty => RGBf(0.45, 0.45, 0.50),
                           :menu => RGBf(0.26, 0.52, 0.96))

"""
    click_feedback!(; flash=true, sound=true)

Turn the visual click effects (button flash, menu flash, ripple in the 3D
view) and the click sound on or off.
"""
function click_feedback!(; flash::Bool=CLICK_FLASH[], sound::Bool=CLICK_SOUND[])
    CLICK_FLASH[], CLICK_SOUND[] = flash, sound
    return (; flash, sound)
end

# --- sound -------------------------------------------------------------------

const CLICK_WAV = Ref{Union{Nothing,String}}(nothing)

"""A 25 ms click (decaying 1.8 kHz tone) as a 16-bit mono WAV file, written once."""
function click_wav()
    CLICK_WAV[] === nothing || isfile(CLICK_WAV[]) || (CLICK_WAV[] = nothing)
    CLICK_WAV[] === nothing || return CLICK_WAV[]
    rate = 44100
    n = round(Int, 0.025rate)
    s = [round(Int16, 9000 * sin(2π * 1800 * t / rate) * exp(-t / (0.005rate))) for t in 0:n-1]
    path = joinpath(tempdir(), "atombond_click.wav")
    open(path, "w") do io
        write(io, "RIFF", UInt32(36 + 2n), "WAVE", "fmt ", UInt32(16), UInt16(1), UInt16(1),
              UInt32(rate), UInt32(2rate), UInt16(2), UInt16(16), "data", UInt32(2n))
        write(io, s)
    end
    return CLICK_WAV[] = path
end

"""Play the click sound in the background (nothing happens without `paplay`)."""
function click_sound()
    CLICK_SOUND[] || return false
    player = Sys.which("paplay")
    player === nothing && return false
    try
        run(pipeline(`$player $(click_wav())`; stdout=devnull, stderr=devnull); wait=false)
    catch err
        err isa InterruptException && rethrow()
        return false
    end
    return true
end

# --- timing --------------------------------------------------------------------

"""
    on_frames(f, scene) -> stop()

Call `f(seconds)` on every drawn frame of `scene`'s window, with the time
since the first such frame, until `f` returns true. Returns a function that
stops it early.
"""
function on_frames(f, scene)
    t0, done = Ref(NaN), Ref(false)
    tick = Makie.events(scene).tick
    listener = Ref{Any}(nothing)
    stop() = done[] || (done[] = true; @async Makie.off(tick, listener[]))  # not while it is notifying
    listener[] = Makie.on(tick) do _
        done[] && return
        isnan(t0[]) && (t0[] = time())
        f(time() - t0[]) && stop()
    end
    return stop
end

"""Call `f()` once `seconds` have passed and at least `frames` frames were drawn since now."""
function after_frames(f, scene; seconds::Real, frames::Integer=2)
    n = Ref(0)
    return on_frames(scene) do t
        (n[] += 1) >= frames && t >= seconds || return false
        f()
        return true
    end
end

# --- buttons -------------------------------------------------------------------

"""
    with_click_feedback!(button) -> button

Flash `button` (and play the click sound) after every click.
"""
function with_click_feedback!(b)
    # The button's own colours, restored after each flash.
    normal = (b.buttoncolor[], b.buttoncolor_hover[], b.labelcolor[], b.labelcolor_hover[])
    stop = Ref{Any}(() -> nothing)
    Makie.on(b.clicks) do _
        click_sound()
        CLICK_FLASH[] || return
        stop[]()                         # a new click restarts the flash
        b.buttoncolor[] = b.buttoncolor_hover[] = FLASH_COLOR
        b.labelcolor[] = b.labelcolor_hover[] = RGBf(1, 1, 1)
        stop[] = after_frames(b.blockscene; seconds=FLASH_TIME) do
            b.buttoncolor[], b.buttoncolor_hover[], b.labelcolor[], b.labelcolor_hover[] = normal
        end
    end
    return b
end

# --- ripple in the 3D view -------------------------------------------------------

"""Overlay for click ripples, above the 3D view (created before any popup menu, so menus stay on top)."""
function install_ripple!(v::Viewer)
    scene = Makie.Scene(v.figure.scene; viewport=v.figure.scene.viewport, clear=false)
    Makie.campixel!(scene)
    points = Makie.Observable([Makie.Point2f(NaN, NaN)])
    color = Makie.Observable(Makie.RGBAf(0, 0, 0, 0))
    Makie.lines!(scene, points; color, linewidth=3.5, fxaa=false, overdraw=true)
    v.dialogs[:ripple] = (; scene, points, color, stop=Ref{Any}(() -> nothing))
    return v
end

"""
    click_ripple!(viewer, xy, kind=:empty)

An expanding, fading ring at window pixel `xy` (`kind` :atom, :empty or
:menu picks the colour).
"""
function click_ripple!(v::Viewer, xy, kind::Symbol=:empty)
    CLICK_FLASH[] && haskey(v.dialogs, :ripple) || return nothing
    r = v.dialogs[:ripple]
    r.stop[]()
    c = RIPPLE_COLORS[kind]
    # From just outside a selected atom's halo outwards, fading.
    circle(ρ) = [Makie.Point2f(xy[1] + ρ * cos(t), xy[2] + ρ * sin(t)) for t in range(0, 2π; length=49)]
    frame(f) = (r.points[] = circle(18 + 40f); r.color[] = Makie.RGBAf(c.r, c.g, c.b, 0.95 * (1 - f)))
    frame(0.0)
    n = Ref(0)
    r.stop[] = on_frames(r.scene) do t
        f = t / RIPPLE_TIME
        (n[] += 1) >= 3 && f >= 1 || (frame(min(f, 1.0)); return false)
        r.points[] = [Makie.Point2f(NaN, NaN)]
        return true
    end
    return nothing
end
