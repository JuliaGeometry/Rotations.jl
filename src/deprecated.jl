
# Deprecate RodriguesVec => RotationVec
abstract type RodriguesVec end
@deprecate RodriguesVec(x, y, z) RotationVec(x, y, z)
@deprecate RodriguesVec(t::NTuple{9}) RotationVec(t::NTuple{9})
@deprecate RodriguesVec(r::Rotation) RotationVec(r::Rotation)
# QUESTION: not sure how to deprecate Base.one(::Type{RodriguesVec})
