
# Deprecate RodriguesVec => RotationVec
struct RodriguesVec{T} end
@deprecate RodriguesVec(x, y, z) RotationVec(x, y, z)
@deprecate RodriguesVec(t::NTuple{9}) RotationVec(t::NTuple{9})
@deprecate RodriguesVec(r::Rotation) RotationVec(r::Rotation)
function RodriguesVec{T}(x, y, z) where T
    Base.depwarn("`RodriguesVec` is deprecated, use `RotationVec` instead", :Type)
    RotationVec{T}(x,y,z)
end
function Base.one(::Type{RodriguesVec})
    Base.depwarn("`one(::Type{RodriguesVec})` is deprecated, use `one(::Type{RotationVec})` instead.", :one)
    one(RotationVec)
end
function Base.one(::Type{RodriguesVec{T}}) where T
    Base.depwarn("`one(::Type{RodriguesVec{T}})` is deprecated, use `one(::Type{RotationVec{T}})` instead.", :one)
    one(RotationVec{T})
end
