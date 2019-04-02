Option{T} = Union{Missing, T}
Vec3 = SVector{3}

unit_vector(v::Vec3) = v ./ norm(v)

function sphere_rand()
    p = Vec3(1.0, 1.0, 1.0)
    while norm(p) >= 1.0
        p = 2.0 .* Vec3(rand(3)...) - Vec3(1.0, 1.0, 1.0)
    end
    p
end

function reflect(v::Vec3, n::Vec3)
    v .- 2 .* (v ⋅ n) .* n
end

function refract(v::Vec3, n::Vec3, nint::Float64)::Option{Vec3}
    uv = unit_vector(v)
    dt = uv ⋅ n
    discriminant = 1.0 - (1-dt^2)*nint^2
    if discriminant > 0
        nint .* (uv .- n .* dt) .- n .* sqrt(discriminant)
    else
        missing
    end
end

function schlick(cosine::Float64, refidx::Float64)
    r0 = ((1-refidx) / (1+refidx)) ^ 2
    r0 + (1-r0)*(1-cosine)^5
end
