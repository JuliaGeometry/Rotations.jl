# 3d
function Base.log(R::RotationVec)
    x, y, z = params(R)
    return @SMatrix [0 -z  y
                     z  0 -x
                    -y  x  0]
end

function Base.log(R::Rotation{3})
    log(RotationVec(R))
end


# 2d
function Base.log(R::Angle2d)
    θ, = params(R)
    return @SMatrix [0 -θ
                     θ  0]
end

function Base.log(R::Rotation{2})
    log(Angle2d(R))
end
