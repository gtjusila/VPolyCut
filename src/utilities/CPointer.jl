# A type which emulates C style pointers
# The idea is in C such calls are possible:
# SCIP* x;
# SCIPcreate(&x);
#
# We emulate this by
# x = CPtr(SCIP); 
# SCIPcreate(address_of(x))

struct CPtr{T}
    pointer::Ref{Ptr{T}}
end

function CPtr(::Type{T}) where {T}
    CPtr(Ref{Ptr{T}}(C_NULL))
end

# Automatic conversion when needed
Base.convert(::Type{Ptr{T}}, cptr::CPtr{T}) where {T} = cptr.pointer[]

# Needed for passing `CPtr` to C functions expecting `Ptr{T}`
Base.unsafe_convert(::Type{Ptr{T}}, cptr::CPtr{T}) where {T} = cptr.pointer[]

function address_of(pt::CPtr)
    return pt.pointer
end