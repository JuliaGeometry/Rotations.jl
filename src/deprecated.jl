
# Deprecate RodriguesVec => RotationVec
Base.@deprecate_binding RodriguesVec RotationVec true

# Deprecate Quat => UnitQuaternion
Base.@deprecate_binding Quat UnitQuaternion true

# Deprecate SPQuat => MRP
Base.@deprecate_binding SPQuat MRP true
