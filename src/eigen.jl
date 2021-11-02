# 3D
function LinearAlgebra.eigvals(R::Rotation{3})
    θ = rotation_angle(R)
    return SVector(exp(-θ*im), exp(θ*im), 1)
end

function LinearAlgebra.eigvecs(R::Rotation{3})
    n3 = normalize(rotation_axis(R))
    n1 = normalize(perpendicular_vector(n3))
    n2 = normalize(n3×n1)
    v1 = normalize(n1 + n2*im)
    v2 = normalize(n1*im + n2)
    v3 = n3
    return hcat(v1,v2,v3)
end

function LinearAlgebra.eigen(R::Rotation{3})
    return eigen(AngleAxis(R))
end

function LinearAlgebra.eigen(R::AngleAxis)
    λs = eigvals(R)
    vs = eigvecs(R)
    return Eigen(λs, vs)
end


# 2D
function LinearAlgebra.eigvals(R::Rotation{2})
    θ = rotation_angle(R)
    return SVector(exp(-θ*im), exp(θ*im))
end

function LinearAlgebra.eigvecs(R::Rotation{2})
    return @SMatrix [1/√2 1/√2;im/√2 -im/√2]
end

function LinearAlgebra.eigen(R::Rotation{2})
    λs = eigvals(R)
    vs = eigvecs(R)
    return Eigen(λs, vs)
end
