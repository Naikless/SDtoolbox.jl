"""
Shock and Detonation Toolbox
"ZND" module

Calculates ZND explosions.

This module defines the following functions:

    znd!
    zndsolve
    getThermicity
    soundspeed_fr


################################################################################
Theory, numerical methods and applications are described in the following report:

    SDToolbox Numerical Tools for Shock and Detonation Wave Modeling,
    Explosion Dynamics Laboratory, Contributors:
    S. Browne, J. Ziegler, N. Bitter, B. Schmidt, J. Lawson and J. E. Shepherd,
    GALCIT Technical Report FM2018.001 Revised January 2021.
    California Institute of Technology, Pasadena, CA USA

Please cite this report and the website if you use these routines.

Please refer to LICENCE.txt or the above report for copyright and disclaimers.

http://shepherd.caltech.edu/EDL/PublicResources/sdt/

"""
module ZND

export zndsolve

using OrdinaryDiffEq
using PyCall


const ct = PyNULL()

function __init__()
    copy!(ct, pyimport_conda("cantera","cantera","cantera"))
    global R̄ = ct.gas_constant
end


U1 = nothing
gas1 = nothing
ρ₁ = nothing
Mᵢ = nothing

"""
    soundspeed_fr(gas::PyObject)

Computes the frozen sound speed by using a forward finite
difference approximation and evaluating frozen composition states on the
isentrope passing through the reference (S, V) state supplied by the gas
object passed to the function.

FUNCTION SYNTAX:
    a_frz =  soundspeed_fr(gas)

INPUT:
    gas = working gas object (restored to original state at end of function)

OUTPUT:
    a_frz = frozen sound speed = sqrt({d P/d rho)_{s,x0})

"""
function soundspeed_fr(gas::PyObject)

    ρ₀ = gas.density
    p₀ = gas.P
    s₀ = gas.entropy_mass
    ρ₁ = 1.001 * ρ₀
    X₀ = gas.X
    # workaround to avoid unsized object error when only one species in a .cti file
    # (flagged to be fixed in future Cantera version)
    if gas.n_species > 1
        gas.SVX =  s₀, 1 ./ ρ₁, X₀
    else
        gas.SV = s₀, 1 ./ ρ₁
    end
    p₁ = gas.P
    dpdρₛ = (p₁ - p₀)/(ρ₁ - ρ₀)
    a_frz = √dpdρₛ

    # Restore gas object to original state
    if gas.n_species > 1
        gas.SVX =  s₀, 1 ./ ρ₀, X₀
    else
        gas.SV = s₀, 1 ./ ρ₀
    end

    return a_frz
end

"""
    znd!(dy::Vector{Float64},y::Vector{Float64},params,t::Real)::Vector{Float64}

Set of ODEs to solve ZND Detonation Problem.

INPUT:
    t = time
    y = solution array [pressure, density, position, species mass 1, 2, ..]
    gas = working gas object
    U₁ = shock velocity (m/s)
    ρ₁ = initial density (kg/m^3)

OUTPUT:
    An array containing time derivatives of:
        pressure, density, distance and species mass fractions,
    formatted in a way that the integrator in zndsolve can recognize.

"""
function znd!(dy::Vector{Float64},y::Vector{Float64},params,t::Real)::Vector{Float64}

    # current gas state and constant parameters
    gas::PyObject,U₁::Real,ρ₁::Real,Mᵢ::Vector{Float64} = params[1:4]

    # make sure physical bounds are respected to avoid cantera errors
    if !any(y[1:2].<=0)
        ρ, p, Y = y[2], y[1], y[4:end]
        gas.DPY = ρ, p, Y
    else # if unphysical, just continue with latest working values and hope step size is adjusted
        ρ, p, Y = gas.DPY::Tuple{Float64, Float64, Vector{Float64}}
        println("out of bounds for ρ or p")
    end

    T = gas.T::Float64
    M̄ = gas.mean_molecular_weight::Float64

    # c = soundspeed_fr(gas)
    # calculating frozen sound speed directly using finite differences saves some PyCalls
    s = gas.entropy_mass::Float64
    ρ′ = 1.001 * ρ
    gas.SVY =  s, 1/ρ′, Y
    p′ = gas.P::Float64
    dpdρₛ = (p′ - p)/(ρ′ - ρ)
    # Restore gas object to original state
    gas.SVY =  s, 1/ρ, Y
    c = √dpdρₛ

    # # alternative frozen sound speed using ideal gas law
    # c = √(gas.cp_mass/gas.cv_mass * T * R̄ / M̄)

    U = U₁*ρ₁/ρ
    Ma = U/c
    η = 1 - Ma^2
    ω̇  = gas.net_production_rates::Vector{Float64}

    Ẏ = ω̇ .* Mᵢ ./ ρ

    # σ̇  = getThermicity(gas)
    # calculating thermicity directly saves some PyCalls
    hₛ_RT = gas.standard_enthalpies_RT
    hₛ = hₛ_RT .* R̄ .* T ./ Mᵢ
    cₚ = gas.cp_mass::Float64
    σ̇  = sum((M̄./Mᵢ .- hₛ./(cₚ * T)) .* Ẏ)

    # Pressure and density derivatives for ODE system
    ṗ  = -ρ * U^2 * σ̇ /η
    ρ̇  = -ρ * σ̇ /η

    # Save current Ma number to be used as a termination criterion.
    # This avoids running into the inevitable singularity at Ma = 1.
    params[end] = Ma

    dy[1] = ṗ
    dy[2] = ρ̇
    dy[3] = U
    dy[4:end] = Ẏ

    return dy
end

"""
    getThermicity(gas::PyObject)

Returns the thermicity = sum ( (w/wi-hsi/(cp*T))*dyidt ).

FUNCTION SYNTAX:
    thermicity = getThermicity(gas)

INPUT:
    gas = Cantera gas object (not modified by this function)

OUTPUT:
    thermicity (1/s)
"""
function getThermicity(gas::PyObject)

    # Mᵢ = gas.molecular_weights
    M̄ = gas.mean_molecular_weight
    T, ρ = gas.TD
    hₛ_RT = gas.standard_enthalpies_RT
    hₛ = hₛ_RT * R̄ * T ./ Mᵢ
    ω̇  = gas.net_production_rates
    cₚ = gas.cp_mass
    dydt = ω̇  .* Mᵢ / ρ
    σ̇ = sum((M̄ ./ Mᵢ - hₛ/(cₚ*T)).*dydt)
end

"""
zndsolve(gas::PyObject,gas₁::PyObject,U₁::Real;
             t_end::Real=1e-3,max_step::Real=1e-4,t_eval=nothing,
             relTol::Real=1e-5,absTol::Real=1e-8,
             advanced_output::Bool=false,solver_algorithm=Rosenbrock23,)

ZND Model Detonation Computation
Solves the set of ODEs defined in znd!().

FUNCTION SYNTAX:
output = zndsolve(gas,gas1,U1,**kwargs)

INPUT
    gas = Cantera gas object - postshock state
    gas₁ = Cantera gas object - initial state
    U₁ = shock velocity (m/s)

OPTIONAL INPUT:
    t_end = end time for integration, in sec
    max_step = maximum time step for integration, in sec
    t_eval = array of time values to evaluate the solution at.
                If left as "None", solver will select values.
                Sometimes these may be too sparse for good-looking plots.
    relTol = relative tolerance
    absTol = absolute tolerance
    advanced_output = calculates optional extra parameters such as induction lengths


OUTPUT:
    output = a dictionary containing the following results:
        time = time array
        distance = distance array

        T = temperature array
        P = pressure array
        rho = density array
        U = velocity array
        thermicity = thermicity array
        species = species mass fraction array

        M = Mach number array
        af = frozen sound speed array
        g = gamma (cp/cv) array
        wt = mean molecular weight array
        sonic = sonic parameter (c^2-U^2) array

        tfinal = final target integration time
        xfinal = final distance reached

        gas₁ = a copy of the input initial state
        U₁ = shock velocity

        and, if advanced_output=True:
        ind_time_ZND = time to maximum thermicity gradient
        ind_len_ZND = distance to maximum thermicity gradient
        exo_time_ZND = pulse width (in secs) of thermicity  (using 1/2 max)
        ind_time_ZND = pulse width (in meters) of thermicity (using 1/2 max)
        max_thermicity_width_ZND = according to Ng et al definition
"""
function zndsolve(gas::PyObject,gas₁::PyObject,U₁::Real;
             t_end::Real=1e-3,max_step::Real=1e-4,t_eval=nothing,
             relTol::Real=1e-5,absTol::Real=1e-8,
             advanced_output::Bool=false,solver_algorithm=Rosenbrock23,)

    ###########################################################
    # Define initial information
    ###########################################################
    println("starting znd solver")
    global gas1 = gas₁
    global U1 = U₁
    global ρ₁ = gas₁.density::Float64

    x_start = 0.
    y₀ = vcat([gas.P,gas.density,x_start],gas.Y)::Vector{Float64}

    global Mᵢ = gas.molecular_weights::Vector{Float64}

    params = [gas,U₁,ρ₁,Mᵢ,0.]

    tel = [0.,t_end] # Timespan

    function approaches_singularity(y::Vector{Float64},t::Float64,integrator)
        """ Returns true if Ma >= 0.99 """
        if integrator.p[end] >= 0.99
            println("Ma = $(integrator.p[end]): Terminating to avoid singularity!")
            return true
        else
            return false
        end
    end

    # Discrete Callback to terminate ODE solver before reaching singularity at Ma = 1
    cb = DiscreteCallback(approaches_singularity,terminate!)
    prob = ODEProblem(znd!,y₀,tel,params)

    @time begin
        # Benchmarks needed for: abstol, reltol, Algos: Rosenbrock23, RadauIIA5, Rodas4
        out = solve(prob,solver_algorithm(autodiff=false),progress=true,callback=cb,abtol=absTol,reltol=relTol)
    end

    create_output_dict(out,gas,advanced_output)
end

function create_output_dict(ode_output::ODESolution,gas::PyObject,advanced_output::Bool=false)

    output = Dict() # could probably be optimized, but for now for convinience

    output["time"] = ode_output.t
    output["P"] = [u[1] for u in ode_output.u]
    output["rho"] = [u[2] for u in ode_output.u]
    output["distance"] = [u[3] for u in ode_output.u]
    output["species"] = [u[4:end] for u in ode_output.u]

    output["tfinal"] = ode_output.t[end]
    output["xfinal"] = output["distance"][end]

    # Initialize additional output matrices where needed
    b = length(output["time"])
    output["T"] = zeros(0)
    output["U"] = zeros(0)
    output["thermicity"] = zeros(0)
    output["af"] = zeros(0)
    output["g"] = zeros(0)
    output["wt"] = zeros(0)
    output["dTdt"] = zeros(0)
    if advanced_output
        output["ind_len_ZND"] = 0
        output["ind_time_ZND"] = 0
        output["exo_len_ZND"] = 0
        output["exo_time_ZND"] = 0
    end

    #############################################################################
    # Extract TEMPERATURE, WEIGHT, GAMMA, SOUND SPEED, VELOCITY, MACH NUMBER,
    # c^2-U^2, THERMICITY, and TEMPERATURE GRADIENT
    #############################################################################

    # Have to loop for operations involving the working gas object
    for (ρ,p,Y) in zip(output["rho"],output["P"],output["species"])
        gas.DPY = ρ,p,Y
        af = soundspeed_fr(gas)
        U = U1*ρ₁/ρ

        push!(output["T"], gas.T)
        push!(output["U"], U)
        push!(output["thermicity"], getThermicity(gas))
        push!(output["af"], af)
        push!(output["g"], gas.cp/gas.cv)
        push!(output["wt"], gas.mean_molecular_weight)
        # push!(output["dTdt"], getTempDeriv(gas,r1,U1))
    end


    #%%

    # Vectorize operations where possible
    output["M"] = output["U"]./output["af"]
    eta = 1 .- output["M"].^2
    output["sonic"] = eta.*output["af"].^2

    if advanced_output
        ################################################################################################
        # Find INDUCTION TIME and LENGTH based on MAXIMUM THERMICITY
        ################################################################################################
        n = argmax(output["thermicity"])

        output["ind_time_ZND"] = output["time"][n]
        output["ind_len_ZND"] = output["distance"][n]
        output["max_thermicity_ZND"] = maximum(output["thermicity"]) # required for Ng et al Chi parameter

        #######################################################
        # Check for eigenvalue detonation
        #######################################################

        if n == b
            println("Error: Maximum thermicity occurs at the end of the reaction zone")
            println("       You may have an eigenvalue detonation, your final integration length may be too short,")
            println("       your mixture may be too rich/lean, or something else may be wrong")
            println(" ")
            println("Mach Number (end of reaction): $(output["M"][b]) - if close to 1, check for eigenvalue detonation")
            output["ind_time_ZND"] = output["time"][b]
            output["ind_len_ZND"] = output["distance"][b]
            output["exo_time_ZND"] = 0
            output["exo_len_ZND"] = 0
            println("Induction Time: $(output["ind_time_ZND"])")
            println("Exothermic Pulse Time: $(output["exo_time_ZND"])")
            return output

        elseif n == 0
            println("Error: Maximum thermicity occurs at the beginning of the reaction zone")
            println("       You may have an eigenvalue detonation, your final integration length may be too short,")
            println("       your mixture may be too rich/lean, or something else may be wrong")
            println(" ")
            println("Mach Number (end of reaction): $(output["M"][b]) - if close to 1, check for eigenvalue detonation")
            output["ind_time_ZND"] = output["time"][0]
            output["ind_len_ZND"] = output["distance"][0]
            output["exo_time_ZND"] = 0
            output["exo_len_ZND"] = 0
            println("Induction Time: $(output["ind_time_ZND"])")
            println("Exothermic Pulse Time: $(output["exo_time_ZND"])")
            return output

        else
            max_sigmadot = maximum(output["thermicity"])
            half_sigmadot_flag1 = 0
            half_sigmadot_flag2 = 0
            # Go into a loop to find two times when sigma_dot is half its maximum
            tstep1 = 1
            tstep2 = 1 # JML temporary
            for (j,thermicity) in enumerate(output["thermicity"])
                # global half_sigmadot_flag1, half_sigmadot_flag2, tstep2, tstep1
                if half_sigmadot_flag1 == 0
                    if thermicity > 0.5*max_sigmadot
                        half_sigmadot_flag1 = 1
                        tstep1 = j
                    end

                elseif half_sigmadot_flag2 == 0
                    if thermicity < 0.5*max_sigmadot
                        half_sigmadot_flag2 = 1
                        tstep2 = j
                    else
                        tstep2 = 1
                    end
                end
            end
        end


        if tstep2 == 1
            print("Error: No pulse in the thermicity")
            print("       You may have an eigenvalue detonation, your final integration length may be too short,")
            print("       your mixture may be too rich/lean, or something else may be wrong")
            output["exo_time_ZND"] = 0
            output["exo_len_ZND"] = 0
        else
            output["exo_time_ZND"] = output["time"][tstep2] - output["time"][tstep1]
            output["exo_len_ZND"] = output["distance"][tstep2] - output["distance"][tstep1]
        end

    end # if advanced_output


    #################################################################
    # Append extra data used to make output file (via znd_fileout)
    output["gas1"] = gas1
    output["U1"] = U1

    return output
end

end #module
