
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
function Base.rand(::Type{RodriguesVec})
    Base.depwarn("`RodriguesVec` is deprecated, use `RotationVec` instead", :Type)
    rand(RotationVec)
end
function Base.rand(::Type{RodriguesVec{T}}) where T
    Base.depwarn("`RodriguesVec` is deprecated, use `RotationVec` instead", :Type)
    rand(RotationVec{T})
end


# Deprecate Quat => UnitQuaternion
struct Quat{T} end
@deprecate Quat(w,x,y,z, normalize::Bool=true) UnitQuaternion(w,x,y,z, normalize)
@deprecate Quat(t::NTuple{9}) UnitQuaternion(t::NTuple{9})
@deprecate Quat(r::Rotation) UnitQuaternion(r::Rotation)
function Quat{T}(r::Rotation) where T
    Base.depwarn("`Quat` is deprecated, use `UnitQuaternion` instead.", :Type)
    UnitQuaternion{T}(r)
end
function Quat{T}(w,x,y,z, normalize::Bool=true) where T
    Base.depwarn("`Quat` is deprecated, use `UnitQuaternion` instead.", :Type)
    UnitQuaternion{T}(w,x,y,z, normalize)
end
function Base.one(::Type{Quat})
    Base.depwarn("`one(::Type{Quat})` is deprecated, use `one(::Type{UnitQuternion})` instead.", :one)
    one(UnitQuaternion)
end
function Base.one(::Type{Quat{T}}) where T
    Base.depwarn("`one(::Type{Quat{T}})` is deprecated, use `one(::Type{UnitQuternion{T}})` instead.", :one)
    one(UnitQuaternion{T})
end
function Base.rand(::Type{Quat})
    Base.depwarn("`Quat` is deprecated, use `UnitQuaternion` instead.", :Type)
    rand(UnitQuaternion)
end
function Base.rand(::Type{Quat{T}}) where T
    Base.depwarn("`Quat` is deprecated, use `UnitQuaternion` instead.", :Type)
    rand(UnitQuaternion{T})
end
Base.@deprecate_binding RodriguesVec RotationVec true
