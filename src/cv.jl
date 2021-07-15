module CV
"""
Shock and Detonation Toolbox
"cv" module

Calculates constant-volume explosions.

This module defines the following functions:

    cvsolve

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


################################################################################
Transfered to Julia language in June 2021 by Niclas Garan, TU Berlin
Tested with:
    Julia 1.6.1 and Cantera 2.5.1
Under these operating systems:
    Windows 10, Linux (CentOS)
"""

export cvsolve

using OrdinaryDiffEq
using PyCall


const ct = PyNULL()

function __init__()
    copy!(ct, pyimport_conda("cantera","cantera","cantera"))
    global R̄ = ct.gas_constant
end


function cv!(dy::Vector{Float64},y::Vector{Float64},params,t::Real)::Vector{Float64}
    """
    Evaluates the system of ordinary differential equations for an adiabatic,
    constant-volume, zero-dimensional reactor.
    It assumes that the "gas" object represents a reacting ideal gas mixture.

    INPUT:
        t = time
        y = solution array [temperature, species mass 1, 2, ...]
        gas = working gas object

    OUTPUT:
        An array containing time derivatives of:
            temperature and species mass fractions,
        formatted in a way that the integrator in cvsolve can recognize.

    """
    # current gas state and constant parameters
    gas::PyObject,ρ₁::Real,Mᵢ::Vector{Float64} = params

    # make sure physical bounds are respected to avoid cantera errors
    if y[1] >= 0
        T, Y = y[1], y[2:end]
        gas.TDY = T, ρ₁, Y
    else # if unphysical, just continue with latest working values and hope step size is adjusted
        T, _, Y = gas.TDY::Tuple{Float64, Float64, Vector{Float64}}
        println("out of bounds for T")
    end

    ω̇  = gas.net_production_rates::Vector{Float64}
    hₛ_RT = gas.standard_enthalpies_RT
    cᵥ = gas.cv_mass::Float64
    Ṫ = -T * R̄ /(ρ₁*cᵥ) * sum((hₛ_RT .- 1) .* ω̇ )
    Ẏ = ω̇  .* Mᵢ ./ ρ₁

    dy[1] = Ṫ
    dy[2:end] = Ẏ

    return dy
end


function cvsolve(gas::PyObject;t_end::Real=1e-6,max_step::Real=1e-5,
                t_eval=nothing,relTol::Real=1e-5,absTol::Real=1e-8)
    """
    Solves the ODE system defined in cv!, taking the gas object input as the
    initial state.


    FUNCTION SYNTAX:
        output = cvsolve(gas,**kwargs)

    INPUT:
        gas = working gas object

    OPTIONAL INPUT:
        t_end = end time for integration, in sec
        max_step = maximum time step for integration, in sec
        t_eval = array of time values to evaluate the solution at.
                    If left as "None", solver will select values.
                    Sometimes these may be too sparse for good-looking plots.
        relTol = relative tolerance
        absTol = absolute tolerances

    OUTPUT:
        output = a dictionary containing the following results:
            time = time array
            T = temperature profile array
            P = pressure profile array
            speciesY = species mass fraction array
            speciesX = species mole fraction array

            gas = working gas object

            exo_time = pulse width (in secs) of temperature gradient (using 1/2 max)
            ind_time = time to maximum temperature gradient
            ind_time_10 = time to 10% of maximum temperature gradient
            ind_time_90 = time to 90% of maximum temperature gradient

    """
    global ρ₁ = gas.density::Float64
    y₀ = vcat(gas.T::Float64,gas.Y::Vector{Float64})
    global Mᵢ = gas.molecular_weights::Vector{Float64}

    tel = [0.,t_end] # Timespan

    # determine adiabatic flame temperature
    old_state = gas.state
    gas.equilibrate("UV")
    T_ad = gas.T::Float64
    gas.state = old_state

    function approaches_equilibrium(y::Vector{Float64},t::Real,integrator)
        """ Returns true if T >= 0.999 * T_ad """
        if y[1] >= 0.999*T_ad && get_du(integrator)[1] <= 1e-6
            println("Reached thermal equilibrium. Terminating!")
            return true
        else
            return false
        end
    end

    params = [gas,ρ₁,Mᵢ]

    # Discrete Callback to terminate ODE solver when approaching equilibrium
    cb = DiscreteCallback(approaches_equilibrium,terminate!)
    prob = ODEProblem(cv!,y₀,tel,params)

    @time begin
        # Benchmarks needed for: abstol, reltol, Algos: Rosenbrock23, RadauIIA5, Rodas4
        out = solve(prob,Rosenbrock23(autodiff=false),progress=true,callback=cb,abtol=absTol,reltol=relTol)
    end

    output = create_output_dict(out,gas)

end


function create_output_dict(ode_output::ODESolution,gas::PyObject)

    output = Dict() # could probably be optimized, but for now for convinience

    output["time"] = ode_output.t
    output["T"] = [u[1] for u in ode_output.u]
    output["speciesY"] = [u[2:end] for u in ode_output.u]

    # Initialize additional output matrices where needed
    b = length(output["time"])
    output["P"] = zeros(0)
    output["dTdt"] = zeros(0)
    output["speciesX"] = Vector{Float64}[]
    output["ind_time"] = 0
    output["ind_time_90"] = 0
    output["ind_time_10"] = 0
    output["exo_time"] = 0
    temp_grad = zeros(0)

    #############################################################################
    # Extract PRESSSURE and TEMPERATURE GRADIENT
    #############################################################################

    # Have to loop for operations involving the working gas object

    for (T,Y) in zip(output["T"],output["speciesY"])
        gas.TDY = T,ρ₁,Y
        M̄ = gas.mean_molecular_weight
        dTdt = cv!(zeros(length(Y)+1),vcat(T,Y),[gas,ρ₁,Mᵢ],0.)[1]::Float64

        push!(output["P"],gas.P)
        push!(output["speciesX"],gas.X)
        push!(output["dTdt"],dTdt)
    end

    n = argmax(output["dTdt"])

    if n == b
        println("Error: Maximum temperature gradient occurs at the end of the reaction zone")
        println("       Your final integration length may be too short,")
        println("       your mixture may be too rich/lean, or something else may be wrong")
        println(" ")
        output["ind_time"] = output["time"][b]
        output["ind_time_10"] = output["time"][b]
        output["ind_time_90"] = output["time"][b]
        output["exo_time"] = 0
        println("Induction Time: $(output["ind_time"])")
        println("Exothermic Pulse Time: $(output["exo_time"])")
        return output
    elseif n == 0
        println("Error: Maximum temperature gradient occurs at the beginning of the reaction zone")
        println("       Your final integration length may be too short,")
        println("       your mixture may be too rich/lean, or something else may be wrong")
        println(" ")
        output["ind_time"] = output["time"][1]
        output["ind_time_10"] = output["time"][1]
        output["ind_time_90"] = output["time"][1]
        output["exo_time"] = 0
        println("Induction Time: $(output["ind_time"])")
        println("Exothermic Pulse Time: $(output["exo_time"])")
        return output
    else
        output["ind_time"] = output["time"][n]

        k = 1
        MAX10 = 0.1*maximum(output["dTdt"])
        d = output["dTdt"][1]
        while d < MAX10 && k < n
            k = k + 1
            d = output["dTdt"][k]
        end
        output["ind_time_10"] = output["time"][k]

        k = 1
        MAX90 = 0.9*maximum(output["dTdt"])
        d = output["dTdt"][1]
        while d < MAX90 && k < n
            k = k + 1
            d = output["dTdt"][k]
        end
        output["ind_time_90"] = output["time"][k]

        # find exothermic time
        tstep1 = 1
        tstep2 = 1
        half_T_flag1 = 0
        half_T_flag2 = 0
        # Go into a loop to find two times when temperature is half its maximum
        for (j,tgrad) in enumerate(output["dTdt"])
            if half_T_flag1 == 0
                if tgrad > 0.5*maximum(output["dTdt"])
                    half_T_flag1 = 1
                    tstep1 = j
                end

            elseif half_T_flag2 == 0
                if tgrad < 0.5*maximum(output["dTdt"])
                    half_T_flag2 = 1
                    tstep2 = j
                else
                    tstep2 = 1
                end
            end
        end
    end


    # Exothermic time for CV explosion
    if tstep2 == 1
        println("Error: No pulse in the temperature gradient")
        println("       Your final integration length may be too short,")
        println("       your mixture may be too rich/lean, or something else may be wrong")
        output["exo_time"] = 0
    else
        output["exo_time"] = output["time"][tstep2] - output["time"][tstep1]
    end

    output["gas"] = gas
    return output
end

end #module
