using PyCall
ct = pyimport("cantera")

P1 = 1e5
T1 = 300
X1 = "H2:42, O2:21,N2:79"
mech = "gri30.xml"

gas = ct.Solution(mech)

T = T1
print("new\n")

function set_T(gas)
    global T *= 1.01
    gas.TPX = T,P1,X1
    print("$T\n")
    print("$(gas.T)\n")
end

for i = 1:100

    # global T *= 1.01
    set_T(gas)

end
