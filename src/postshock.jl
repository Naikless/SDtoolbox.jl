module Postshock
"""
Shock and Detonation Toolbox
"Postshock" module

Calculates CJ detonation speed and post-shock states for frozen and equilibrium cases.

This module defines the following functions:

    LSQ_CJspeed
    hug_fr
    hug_eq
    FHFP
    CJ_calc
    CJspeed
    PostShock_fr
    PostShock_eq
    shock_calc
    shk_eq_calc

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

export CJspeed, PostShock_eq, PostShock_fr

using PyCall
const ct = PyNULL()

function __init__()
    copy!(ct, pyimport_conda("cantera","cantera","cantera"))
    global R̄ = ct.gas_constant
end

# constants from sdtoolbox/config
const ERRFT = 1e-4
const ERRFV = 1e-4
const volumeBoundRatio = 5


function LSQ_CJspeed(x,y)
    """
    Determines least squares fit of parabola to input data

    FUNCTION SYNTAX:
    [a,b,c,R2,SSE,SST] = LSQ_CJspeed(x,y)

    INPUT:
        x = independent data points
        y = dependent data points

    OUTPUT:
        a,b,c = coefficients of quadratic function (ax^2 + bx + c = 0)
        R2 = R-squared value
        SSE = sum of squares due to error
        SST = total sum of squares

    """
    # Calculate Sums
    X = 0.0; X2 = 0.0; X3 = 0.0; X4 = 0.0;
    Y = 0.0; Y1 = 0.0; Y2 = 0.0;
    a = 0.0; b = 0.0; c = 0.0; R2 = 0.0
    n = length(x)

    for k = 1:n
        X += x[k]
        X2 += x[k]^2
        X3 += x[k]^3
        X4 += x[k]^4
        Y += y[k]
        Y1 += y[k]*x[k]
        Y2 += y[k]*x[k]^2
    end
    m = Y/n

    den = X3*n - X2*X
    temp = den*(X*X2-X3*n) + X2*X2*(X*X-n*X2) - X4*n*(X*X-X2*n)
    temp2 = den*(Y*X2-Y2*n) + (Y1*n-Y*X)*(X4*n-X2*X2)

    b = temp2/temp
    a = 1.0/den*(n*Y1 - Y*X - b*(X2*n-X*X))
    c = 1/n*(Y - a*X2 - b*X)

    SSE = 0.0; SST = 0.0;

    for k = 1:n
        f = a*x[k]^2 + b*x[k] + c
        SSE += (y[k] - f)^2
        SST += (y[k] - m)^2
    end
    R2 = 1 - SSE/SST

    return a,b,c,R2,SSE,SST
end


function hug_fr(x,vb,h1,P1,v1,gas)
    """
    Computes difference in enthalpy computed from assumed (T, V)
    state and fixed composition (frozen) hugoniot evaluation.  Used
    with root solver such as 'fsolve' to compute frozen hugoniot as a function
    of volume.  Input gas object is modified to correspond to input (T, V).

    FUNCTION SYNTAX:
        diff = hug_fr(x,vb,h1,P1,v1,gas)

    USAGE:
        fval = fsolve(hug_fr,Ta,args=(vb,h1,P1,v1,gas))
            = frozen Hugoniot temperature (K) corresponding to vb

    INPUT:
        Ta = initial guess for frozen Hugoniot temperature (K)
        vb = desired frozen Hugoniot specific volume (m^3/kg)
        h1 = enthalpy at state 1 (J/kg)
        P1 = pressure at state 1 (Pa)
        v1 = specific volume at state 1 (m^3/kg)
        gas = working gas object

    OUTPUT:
        enthalpy difference

    """
    gas.TD = x, 1.0/vb
    hb1 = gas.enthalpy_mass
    Pb = R̄ * x/(gas.mean_molecular_weight*vb)

    hb2 = h1 + 0.5*(Pb-P1)*(vb+v1)
    return hb2-hb1
end


function hug_eq(x,vb,h1,P1,v1,gas)
    """
    Computes difference in enthalpy computed from assumed (T, V)
    state and equilibrium composition state hugoniot evaluation.  Used
    with root solver such as 'fsolve' to compute equilibrium hugoniot as a function
    of volume.  Input gas object is modified to correspond to input (T, V) and
    an equilibrium composition.

    FUNCTION SYNTAX:
        diff = hug_eq(x,vb,h1,P1,v1,gas)

    USAGE:
        fval = fsolve(hug_eq,Ta,args=(vb,h1,P1,v1,gas))
            = equilibrium Hugoniot temperature (K) corresponding to vb

    INPUT:
        Ta = initial guess for equilibrium Hugoniot temperature (K)
        vb = desired equilibrium Hugoniot specific volume (m^3/kg)
        h1 = enthalpy at state 1 (J/kg)
        P1 = pressure at state 1 (Pa)
        v1 = specific volume at state 1 (m^3/kg)
        gas = working gas object

    OUTPUT:
        enthalpy difference

    """
    gas.TD = x, 1.0/vb
    gas.equilibrate("TV")
    hb1 = gas.enthalpy_mass
    Pb = R̄ * x / (gas.mean_molecular_weight*vb)
    hb2 = h1 + 0.5*(Pb-P1)*(vb+v1)
    return hb2-hb1
end


function FHFP(w1,gas2,gas1)
    """
    Uses the momentum and energy conservation equations to calculate
    error in pressure and enthalpy given shock speed, upstream (gas1)
    and downstream states (gas2).  States are not modified by these routines.

    FUNCTION SYNTAX:
        [FH,FP] = FHFP(w1,gas2,gas1)

    INPUT:
        w1 = shock speed (m/s)
        gas2 = gas object at working/downstream state
        gas1 = gas object at initial/upstream state

    OUTPUT:
        FH,FP = error in enthalpy and pressure

    """
    P1 = gas1.P
    H1 = gas1.enthalpy_mass
    r1 = gas1.density
    P2 = gas2.P
    H2 = gas2.enthalpy_mass
    r2 = gas2.density
    w1s = w1^2
    w2s = w1s*(r1/r2)^2
    FH = H2 + 0.5*w2s - (H1 + 0.5*w1s)
    FP = P2 + r2*w2s - (P1 + r1*w1s)
    return FH, FP
end


function CJ_calc(gas, gas1, ERRFT, ERRFV, x)
    """
    Calculates the Chapman-Jouguet wave speed using Reynolds' iterative method.

    FUNCTION SYNTAX:
        [gas,w1] = CJ_calc(gas,gas1,ERRFT,ERRFV,x)

    INPUT:
        gas = working gas object
        gas1 = gas object at initial state
        ERRFT,ERRFV = error tolerances for iteration
        x = density ratio

    OUTPUT:
        gas = gas object at equilibrium state
        w1 = initial velocity to yield prescribed density ratio

    """
    T = 2000; r1 = gas1.density
    V1 = 1/r1
    i = 0; DT = 1000; DW = 1000;
    #PRELIMINARY GUESS
    V = V1/x; r = 1/V; w1 = 2000
    gas.TD = T, r
    gas.equilibrate("TV")
    H, P = gas.HP
    # START LOOP
    while (abs(DT) > ERRFT*T) || (abs(DW) > ERRFV*w1)
        i += 1
        if i == 500
            println("CJ_calc did not converge")
            return
        end
        # CALCULATE FH & FP FOR GUESS 1
        FH,FP = FHFP(w1,gas,gas1)
        # TEMPERATURE PERTURBATION
        DT = T*0.02; Tper = T + DT
        Vper = V; Rper = 1/Vper
        Wper = w1
        gas.TD = Tper, Rper
        gas.equilibrate("TV")
        # CALCULATE FHX & FPX FOR "IO" STATE
        FHX,FPX = FHFP(Wper,gas,gas1)
        #ELEMENTS OF JACOBIAN
        DFHDT = (FHX-FH)/DT; DFPDT = (FPX-FP)/DT
        # VELOCITY PERTURBATION
        DW = 0.02*w1; Wper = w1 + DW
        Tper = T; Rper = 1/V
        gas.TD = Tper, Rper
        gas.equilibrate("TV")
        # CALCULATE FHX & FPX FOR "IO" STATE
        FHX,FPX = FHFP(Wper,gas,gas1)
        #ELEMENTS OF JACOBIAN
        DFHDW = (FHX-FH)/DW; DFPDW = (FPX-FP)/DW
        # INVERT MATRIX
        J = DFHDT*DFPDW - DFPDT*DFHDW
        b = DFPDW, -DFHDW, -DFPDT, DFHDT
        a = -FH, -FP
        DT = (b[1]*a[1]+b[2]*a[2])/J; DW = (b[3]*a[1]+b[4]*a[2])/J
        # CHECK & LIMIT CHANGE VALUES
        # VOLUME
        DTM = 0.2*T
        if abs(DT) > DTM
            DT = DTM*DT/abs(DT)
        end
        # MAKE THE CHANGES
        T += DT; w1 += DW
        gas.TD = T, r
        gas.equilibrate("TV")
    end
    return w1
end


function CJspeed(P1, T1, q, mech; method = "CEA", kwargs...)
    """ defaults to CEA Algorithm """
    if method == "CEA"
        _CJspeed_CEA(P1, T1, q, mech)
    elseif method == "umin"
        _CJspeed_umin(P1, T1, q, mech; kwargs...)
    elseif method == "aeq"
        _CJspeed_aeq(P1, T1, q, mech)
    end
end


function _CJspeed_umin(P1, T1, q, mech; fullOutput=false)
    """
    Calculates CJ detonation velocity for a given pressure, temperature, and
    composition.

    FUNCTION SYNTAX:
        If only CJ speed required:
        cj_speed = CJspeed(P1,T1,q,mech)
        If full output required:
        [cj_speed,R2,plot_data] = CJspeed(P1,T1,q,mech,fullOutput=true)

    INPUT:
        P1 = initial pressure (Pa)
        T1 = initial temperature (K)
        q = reactant species mole fractions in one of Cantera's recognized formats
        mech = cti file containing mechanism data (e.g. 'gri30.cti')

    OPTIONAL INPUT:
        fullOutput = set true for R-squared value and pre-formatted plot data
                    (the latter for use with sdtoolbox.utilities.CJspeed_plot)

    OUTPUT
        cj_speed = CJ detonation speed (m/s)
        R2 = R-squared value of LSQ curve fit (optional)
        plot_data = tuple (rr,w1,dnew,a,b,c)
                    rr = density ratio
                    w1 = speed
                    dnew = minimum density
                    a,b,c = quadratic fit coefficients

    """
    #DECLARATIONS
    numsteps = 20; maxv = 2.0; minv = 1.5
    w1 = zeros(Float64, numsteps+1)
    rr = zeros(Float64,numsteps+1)
    gas1 = ct.Solution(mech)
    gas  = ct.Solution(mech)
    #INTIAL CONDITIONS
    gas.TPX  = T1, P1, q
    gas1.TPX = T1, P1, q
    #INITIALIZE ERROR VALUES & CHANGE VALUES
    ERRFT = 1e-4;  ERRFV = 1e-4

    T1 = gas1.T; P1 = gas1.P
    counter = 1; R2 = 0.0; cj_speed = 0.0
    a = 0.0; b = 0.0; c = 0.0; dnew = 0.0
    while (counter <= 4) || (R2 < 0.99999)
        step = (maxv-minv)/numsteps; i = 1; x = minv
        while x <= maxv
            gas.TPX = T1, P1, q
            w1[i] = CJ_calc(gas, gas1, ERRFT, ERRFV, x)
            rr[i] = gas.density/gas1.density
            i += 1; x += step
        end
        a,b,c,R2,SSE,SST = LSQ_CJspeed(rr,w1)
        dnew = -b/(2.0*a)
        minv = dnew * (1 - 0.001)
        maxv = dnew * (1 + 0.001)
        counter += 1
        cj_speed = a*dnew^2 + b*dnew + c
    end

    if fullOutput
        # Optional output data for plotting (with sdtoolbox.utilities.CJspeed_plot)
        plot_data = rr,w1,dnew,a,b,c
        return cj_speed,R2,plot_data
    else
        return cj_speed
    end
end


function PostShock_fr(U1, P1, T1, q, mech)
    """
    Calculates frozen post-shock state for a specified shock velocity and pre-shock state.

    FUNCTION SYNTAX:
        gas = PostShock_fr(U1,P1,T1,q,mech)

    INPUT:
        U1 = shock speed (m/s)
        P1 = initial pressure (Pa)
        T1 = initial temperature (K)
        q = reactant species mole fractions in one of Cantera's recognized formats
        mech = cti file containing mechanism data (e.g. 'gri30.cti')

    OUTPUT:
        gas = gas object at frozen post-shock state

    """

    gas1 = ct.Solution(mech)
    gas  = ct.Solution(mech)
    # INTIAL CONDITIONS
    gas.TPX = T1, P1, q
    gas1.TPX = T1, P1, q
    # CALCULATES POST-SHOCK STATE
    gas = shk_calc(U1, gas, gas1, ERRFT, ERRFV)
end


function PostShock_eq(U1, P1, T1, q, mech)
    """
    Calculates equilibrium post-shock state for a specified shock velocity and pre-shock state.

    FUNCTION SYNTAX:
        gas = PostShock_eq(U1,P1,T1,q,mech)

    INPUT:
        U1 = shock speed (m/s)
        P1 = initial pressure (Pa)
        T1 = initial temperature (K)
        q = reactant species mole fractions in one of Cantera's recognized formats
        mech = cti file containing mechanism data (e.g. 'gri30.cti')

    OUTPUT:
        gas = gas object at equilibrium post-shock state

    """

    gas1 = ct.Solution(mech);
    gas  = ct.Solution(mech);
    # INTIAL CONDITIONS
    # workaround to avoid unsized object error when only one species in a .cti file
    # (flagged to be fixed in future Cantera version)
    if length(q) > 1
        gas.TPX = T1, P1, q
        gas1.TPX= T1, P1, q
    else
        gas.TP = T1, P1
        gas1.TP= T1, P1
    end
    # CALCULATES POST-SHOCK STATE
    gas = shk_eq_calc(U1, gas, gas1, ERRFT, ERRFV)
end


function shk_calc(U1, gas, gas1, ERRFT, ERRFV)
    """
    Calculates frozen post-shock state using Reynolds' iterative method.

    FUNCTION SYNTAX:
        gas = shk_calc(U1,gas,gas1,ERRFT,ERRFV)

    INPUT:
        U1 = shock speed (m/s)
        gas = working gas object
        gas1 = gas object at initial state
        ERRFT,ERRFV = error tolerances for iteration

    OUTPUT:
        gas = gas object at frozen post-shock state

    """

    r1 = gas1.density; V1 = 1/r1
    P1 = gas1.P; T1 = gas1.T
    i = 0
    deltaT = 1000; deltaV = 1000
    # PRELIMINARY GUESS
    Vg = V1/volumeBoundRatio; rg = 1/Vg
    Pg = P1 + r1*(U1^2)*(1-Vg/V1); Tg = T1*Pg*Vg/(P1*V1)
    gas.TD = Tg, rg
    Hg, Pg = gas.HP
    # SAVE STATE
    V = Vg; r = rg; P = Pg
    T = Tg; H = Hg
    # START LOOP
    while (abs(deltaT) > ERRFT*T) || (abs(deltaV) > ERRFV*V)
        i += 1
        if i == 500
            println("shk_calc did not converge for U = $U1")
            return gas
        end
        # CALCULATE FH & FP FOR GUESS 1
        FH,FP = FHFP(U1,gas,gas1)

        # TEMPERATURE PERTURBATION
        DT = T*0.02; Tper = T + DT;
        Vper = V; Rper = 1/Vper;
        gas.TD = Tper, Rper
        # CALCULATE FHX & FPX FOR "IO" STATE
        FHX,FPX = FHFP(U1,gas,gas1)
        # ELEMENTS OF JACOBIAN
        DFHDT = (FHX-FH)/DT; DFPDT = (FPX-FP)/DT

        # VOLUME PERTURBATION
        DV = 0.02*V; Vper = V + DV
        Tper = T; Rper = 1/Vper
        gas.TD = Tper, Rper
        # CALCULATE FHX & FPX FOR "IO" STATE
        FHX,FPX = FHFP(U1,gas,gas1)
        # ELEMENTS OF JACOBIAN
        DFHDV = (FHX-FH)/DV; DFPDV = (FPX-FP)/DV

        # INVERT MATRIX
        J = DFHDT*DFPDV - DFPDT*DFHDV
        b = DFPDV, -DFHDV, -DFPDT, DFHDT
        a = -FH, -FP
        deltaT = (b[1]*a[1]+b[2]*a[2])/J; deltaV = (b[3]*a[1]+b[4]*a[2])/J

        # CHECK & LIMIT CHANGE VALUES
        # TEMPERATURE
        DTM = 0.2*T
        if abs(deltaT) > DTM
            deltaT = DTM*deltaT/abs(deltaT)
        end
        # VOLUME
        V2X = V + deltaV
        if V2X > V1
            DVM = 0.5*(V1 - V)
        else
            DVM = 0.2*V
        end
        if abs(deltaV) > DVM
            deltaV = DVM*deltaV/abs(deltaV)
        end
        # MAKE THE CHANGES
        T += deltaT; V += deltaV; r = 1/V
        gas.TD = T, r
    end
    return gas
end


function shk_eq_calc(U1, gas, gas1, ERRFT, ERRFV)
    """
    Calculates equilibrium post-shock state using Reynolds' iterative method.

    FUNCTION SYNTAX:
        gas = shk_calc(U1,gas,gas1,ERRFT,ERRFV)

    INPUT:
        U1 = shock speed (m/s)
        gas = working gas object
        gas1 = gas object at initial state
        ERRFT,ERRFV = error tolerances for iteration

    OUTPUT:
        gas = gas object at equilibrium post-shock state

    """

    r1 = gas1.density; V1 = 1/r1
    P1 = gas1.P; T1 = gas1.T
    i = 0
    deltaT = 1000; deltaV = 1000
    # PRELIMINARY GUESS
    V = V1/volumeBoundRatio; r = 1/V
    P = P1 + r1*(U1^2)*(1-V/V1); T = T1*P*V/(P1*V1)
    gas.TD = T,r
    H,P = gas.HP
    # START LOOP
    while (abs(deltaT) > ERRFT*T) || (abs(deltaV) > ERRFV*V)
        i += 1
        if i == 500
            println("shk_calc did not converge for U = $U1")
            return gas
        end
        # CALCULATE FH & FP FOR GUESS 1
        FH,FP = FHFP(U1,gas,gas1)

        # TEMPERATURE PERTURBATION
        DT = T*0.02; Tper = T + DT
        Vper = V; Rper = 1/Vper
        gas.TD = Tper, Rper
        gas.equilibrate("TV")
        # CALCULATE FHX & FPX FOR "IO" STATE
        FHX,FPX = FHFP(U1,gas,gas1)
        # ELEMENTS OF JACOBIAN
        DFHDT = (FHX-FH)/DT; DFPDT = (FPX-FP)/DT

        # VOLUME PERTURBATION
        DV = 0.02*V; Vper = V + DV
        Tper = T; Rper = 1/Vper
        gas.TD = Tper, Rper
        gas.equilibrate("TV")
        # CALCULATE FHX & FPX FOR "IO" STATE
        FHX,FPX = FHFP(U1,gas,gas1)
        # ELEMENTS OF JACOBIAN
        DFHDV = (FHX-FH)/DV; DFPDV = (FPX-FP)/DV

        # INVERT MATRIX
        J = DFHDT*DFPDV - DFPDT*DFHDV
        b = DFPDV, -DFHDV, -DFPDT, DFHDT
        a = -FH, -FP
        deltaT = (b[1]*a[1]+b[2]*a[2])/J; deltaV = (b[3]*a[1]+b[4]*a[2])/J

        # CHECK & LIMIT CHANGE VALUES
        # TEMPERATURE
        DTM = 0.2*T
        if abs(deltaT) > DTM
            deltaT = DTM*deltaT/abs(deltaT)
        end
        # VOLUME
        V2X = V + deltaV
        if V2X > V1
            DVM = 0.5*(V1 - V)
        else
            DVM = 0.2*V
        end
        if abs(deltaV) > DVM
            deltaV = DVM*deltaV/abs(deltaV)
        end
        # MAKE THE CHANGES
        T += deltaT; V += deltaV; r = 1/V
        gas.TD = T,r
        gas.equilibrate("TV")
    end

    return gas
end


function _CJspeed_aeq(P1, T1, q, mech)

    """

    CJspeed_aeq
    Calculates CJ detonation velocity and CJ state based on equilibrium sound speed

    FUNCTION
    SYNTAX
    cj_speed,gas = CJspeed(P1,T1,q,mech)

    INPUT
    P1 = initial pressure (Pa)
    T1 = initial temperature (K)
    q = string of reactant species mole fractions
    mech = cti file containing mechanism data (i.e. 'gri30.cti')

    OUTPUT
    cj_speed = CJ detonation speed (m/s)

    """

    gas2 = ct.Solution(mech)
    gas1 = ct.Solution(mech)
    gas  = ct.Solution(mech)

    #INTIAL CONDITIONS
    gas.TPX = T1, P1, q
    gas1.TPX = T1, P1, q

    #INITIALIZE ERROR VALUES & CHANGE VALUES
    ERRFT = 1e-4;  ERRFV = 1e-4

    r1 = gas1.density; V1 = 1/r1
    P1 = gas1.P; T1 = gas1.T
    i = 0
    #PRELIMINARY GUESS
    Vg = V1/5; rg = 1/Vg

    gas.TD = T1,rg
    gas.equilibrate("UV")
    Tg = gas.T
    gas2.TDX = Tg, rg, gas.X

    #SAVE STATE
    V = Vg; r = rg
    T = Tg
    deltaT = 1000; deltaV = 1000; cj_speed = 0
    #START LOOP
    while (abs(deltaT) > ERRFT*T) || (abs(deltaV) > ERRFV*V)
        i += 1
        if i == 500
            println("CJspeed_sound did not converge")
            return
        end

        #CALCULATE FH & FP FOR GUESS 1
        FH,FP,cj_speed = FHFP_CJ2(gas,gas1,gas2)

        #TEMPERATURE PERTURBATION
        DT = T*0.01; Tper = T + DT
        Vper = V; Rper = 1/Vper

        gas.TD = Tper, Rper
        gas.equilibrate("TV")
        gas2.TDX = Tper, Rper, gas.X

        #CALCULATE FHX & FPX FOR "IO" STATE
        FHX,FPX,cj_speed = FHFP_CJ2(gas,gas1,gas2)
        #ELEMENTS OF JACOBIAN
        DFHDT = (FHX-FH)/DT; DFPDT = (FPX-FP)/DT

        #VOLUME PERTURBATION
        DV = 0.01*V; Vper = V + DV
        Tper = T; Rper = 1/Vper

        gas.TD = Tper, Rper
        gas.equilibrate("TV")
        gas2.TDX = Tper, Rper, gas.X

        #CALCULATE FHX & FPX FOR "IO" STATE
        FHX,FPX,cj_speed = FHFP_CJ2(gas,gas1,gas2)
        #ELEMENTS OF JACOBIAN
        DFHDV = (FHX-FH)/DV; DFPDV = (FPX-FP)/DV

        #INVERT MATRIX
        J = DFHDT*DFPDV - DFPDT*DFHDV
        b = DFPDV, -DFHDV, -DFPDT, DFHDT
        a = -FH, -FP
        deltaT = (b[1]*a[1]+b[2]*a[2])/J; deltaV = (b[3]*a[1]+b[4]*a[2])/J

        #CHECK & LIMIT CHANGE VALUES
        #TEMPERATURE
        DTM = 0.2*T
        if abs(deltaT) > DTM
            deltaT = DTM*deltaT/abs(deltaT)
        end
        #VOLUME
        V2X = V + deltaV
        if V2X > V1
            DVM = 0.5*(V1 - V)
        else
            DVM = 0.2*V
        end
        if abs(deltaV) > DVM
            deltaV = DVM*deltaV/abs(deltaV)
        end
        #MAKE THE CHANGES
        T += deltaT; V += deltaV; r = 1/V
        gas.TD = T, r
        gas.equilibrate("TV")
        gas2.TDX = T, r, gas.X
    end

    FH,FP,cj_speed = FHFP_CJ2(gas,gas1,gas2)

    return cj_speed
end


function FHFP_CJ2(gas,gas1,gas2)

    """

    FHFP_CJ2
    Uses the momentum and energy conservation equations and the equilibrium sound speed to calculate error in current pressure and enthalpy guesses.  In this case, state 2 is in equilibrium.

    FUNCTION
    SYNTAX
    [FH,FP,cj_speed] = FHFP_CJ2(gas,gas1,gas2)

    INPUT
    gas = working gas object
    gas1 = gas object at initial state
    gas2 = dummy gas object (for calculating numerical derivatives)

    OUTPUT
    FH,FP = error in enthalpy and pressure
    cj_speed = CJ detonation speed (m/s)

    """

    P1 = gas1.P
    H1 = gas1.enthalpy_mass
    r1 = gas1.density
    P2 = gas.P
    H2 = gas.enthalpy_mass
    r2 = gas.density

    speeds = equilSoundSpeeds(gas2)
    w2s = speeds[1]^2
    w1s = w2s*(r2/r1)^2
    FH = H2 + 0.5*w2s - (H1 + 0.5*w1s)
    FP = P2 + r2*w2s - (P1 + r1*w1s)
    return FH, FP, sqrt(w1s)
end

function equilSoundSpeeds(gas)

    """

    equilSoundSpeeds
    Calculates equilibrium and frozen sound speeds. For the equilibrium sound speed, the gas is equilibrated holding entropy and specific volume constant.

    FUNCTION
    SYNTAX
    [aequil,afrozen] = equilSoundSpeeds(gas)

    INPUT
    gas = working gas object (modified inside function)

    OUTPUT
    aequil = equilibrium sound speed (m/s)
    afrozen = frozen sound speed (m/s)

    """

    # set the gas to equilibrium at its current T and P
    gas.equilibrate("TP")

    # save properties
    s0 = gas.entropy_mass
    p0 = gas.P
    r0 = gas.density

    # perturb the density
    r1 = r0*1.0001

    # set the gas to a state with the same entropy and composition but
    # the perturbed density
    gas.SV = s0, 1.0/r1

    # save the pressure for this case for the frozen sound speed
    pfrozen = gas.P

    # now equilibrate the gas holding S and V constant
    gas.equilibrate("SV")

    p1 = gas.P

    # equilibrium sound speed
    aequil = sqrt((p1 - p0)/(r1 - r0))

    # frozen sound speed
    afrozen = sqrt((pfrozen - p0)/(r1 - r0));
    return aequil, afrozen
end

function _CJspeed_CEA(P1,T1,X1,mech)
    """
    This function calculates the CJ-detonation velocity for given initial
    conditions.

    It is based on (i.e. mostly copied from) the CJ-detonation part of the
    NASA CEA Fortran95 code, available at
    http://www.grc.nasa.gov/WWW/CEAWeb/ceaHome.htm,
    which itself is based on "Calculation of Detonation Properties and Effect
    of Independent Parameters on Gaseous Detonations" by Frank J. Zeleznik
    and Sanford Gordon, ARS Journal, April 1962.

    INPUT:
    P1 = pressure in (Pa)
    T1 = temperature in (K)
    X1 = mole fraction column vector
    mech = file name of kinetic mechanism, e.g. 'gri30.xml'

    OUTPUT:
    cj_speed = CJ detonation velocity in (m/s)
    """


    gas = ct.Solution(mech)
    gas.TPX = T1,P1,X1


    M1 = gas.mean_molecular_weight
    h1 = gas.enthalpy_mass

    P = 15 * P1

    # Acquire initial guess for T2 by finding flame temperature according to
    # chapter 8.3 of CEA documentation
    h_init = h1 + 0.75 * R̄ * T1/M1 * P/P1
    gas.HP = h_init, P
    gas.equilibrate("HP")
    T_init = gas.T
    T = T_init

    # Further improving initial guesses for T and P by recursion formulas
    for i = 1:3
        gas.TPX = T,P,X1
        gas.equilibrate("TP")
        M2 = gas.mean_molecular_weight
        α = T1/T * M2/M1
        γₛ, cₚ_eq, Dlvtp, Dlvpt = partDeriv(T,P,X1,gas)
        P = P1 * (1+γₛ) / (2*γₛ*α) * (1 + (1- (4*γₛ*α)/(1+γₛ)^2)^0.5)
        r = α * P/P1
        T = T1 * T_init/T1 - 0.75 * R̄/(M1*cₚ_eq) * 15 + (R̄*γₛ)/(2*M1*cₚ_eq) * (r^2-1)/r * P/P1
    end

    # Main iteration loop
    cj_speed = 0
    I = 0
    while true
        I += 1
        gas.TPX = T,P,X1
        gas.equilibrate("TP")
        M2 = gas.mean_molecular_weight
        γₛ, cₚ_eq, Dlvtp, Dlvpt = partDeriv(T,P,X1,gas)
        α = T1/T * M2/M1
        r = α * P/P1
        a₁₁ = P1/P + γₛ*r*Dlvpt
        a₁₂ = γₛ*r*Dlvtp
        a₂₁ = 0.5*γₛ*(r^2 - 1 - Dlvpt*(1+r^2)) + Dlvtp - 1
        a₂₂ = -0.5*γₛ*Dlvtp*(r^2+1) - M2*cₚ_eq/R̄
        b₁ = P1/P - 1 + γₛ*(r-1)
        b₂ = M2*(gas.enthalpy_mass-h1)/(R̄*T) - 0.5*γₛ*(r^2-1)
        d = a₁₁*a₂₂ - a₁₂*a₂₁
        x₁ = (a₂₂*b₁-a₁₂*b₂)/d
        x₂ = (a₁₁*b₂-a₂₁*b₁)/d

        δ = 1
        temp = maximum((abs(x₁), abs(x₂)))
        if temp > 0.4054652
            δ = 0.4054652/temp
        end
        P = P*exp(x₁*δ)
        T = T*exp(x₂*δ)
        soundspeed_equil = sqrt(R̄*γₛ*T/M2)
        cj_speed = r*soundspeed_equil

        # Convergence test
        if I > 8
            println("Could not converge in 8 iterations")
            break
        elseif I <= 8 && temp < 0.5e-4
            break
        end
    end
    return cj_speed
end

function partDeriv(T,P,X,gas)
    """
    Calculates all needed variables related to partial derivatives, i.e
    γₛ (isentropic exponent), cₚ_eq (equilibrium specific heat at
    constant pressure), Dlvtp (d ln(v) / d ln(T) at constant pressure),
    Dlvpt (d ln(v) / d ln(P) at constant temperature). Derivatives are
    approximated via central differences.
    """

    v = 1/gas.density
    ΔT = T*1e-4
    ΔP = P*1e-4

    v_P = zeros(2)
    v_T = zeros(2)
    h = zeros(2)

    sign = -1, 1

    for i = 1:2

        gas.TPX = T, P + sign[i]*ΔP, X
        gas.equilibrate("TP")
        v_P[i] = 1/gas.density

        gas.TPX = T + sign[i]*ΔT, P, X
        gas.equilibrate("TP")
        v_T[i] = 1/gas.density
        h[i] = gas.enthalpy_mass
    end

    Dlvtp = T/v * (v_T[2] - v_T[1]) / (2*ΔT)
    Dlvpt = P/v * (v_P[2] - v_P[1]) / (2*ΔP)
    cₚ_eq = (h[2]-h[1])/(2*ΔT)

    gas.TPX = T,P,X
    gas.equilibrate("TP")

    γₛ = -1/(Dlvpt + R̄ * Dlvtp^2 / (cₚ_eq * gas.mean_molecular_weight))

    return γₛ, cₚ_eq, Dlvtp, Dlvpt
end

end #module
