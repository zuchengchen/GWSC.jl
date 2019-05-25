
abstract type Detector end

# spedd of light in vaccum (m/s)
const c0 = ustrip(float(c_0))

# Julian year in seconds
const YEAR = ustrip(uconvert(Unitful.s, 1*u"yr"))

# Day in seconds
const DAY = YEAR/365.25

# Hubble constant (s^-1)
const H0 = ustrip(uconvert(u"s"^-1, 67.66 * u"km" * u"Mpc"^-1 * u"s"^-1))

# Astronomical unit (meters)
# const AU = ustrip(uconvert(Unitful.m, 1*u"AU"))

const Float = Float64

"""
    backup(array1, array2, fileName::String)

Write the contents of two arrays to a file.

`array1`: the first column of the file;

`array2`: the second column of the file;

`fileName`: name of the file to be written to.
"""
function backup(array1, array2, fileName::String)
    
    if length(array1) != length(array2)
        print("The two arrays do not have the same number of elements.
            I will quit now!")
        return
    end
    
    len = length(array1)
    open(fileName, "w") do io
        for i in 1:len
            write(io, string(array1[i])*"  "*string(array2[i])*"\n")
        end
    end
end