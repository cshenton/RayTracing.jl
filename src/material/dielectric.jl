struct Dielectric <: Material
    refidx::Float64
end

function scatter(d::Dielectric, rayin::Ray, rec::HitRecord)::Scatter
    reflected = reflect(rayin.direction, rec.normal)
    if rayin.direction ⋅ rec.normal > 0.0
        normal = -rec.normal
        nint = d.refidx
        cosine = d.refidx * (rayin.direction ⋅ rec.normal) / norm(rayin.direction)
    else
        normal = rec.normal
        nint = 1.0 / d.refidx
        cosine = -(rayin.direction ⋅ rec.normal) / norm(rayin.direction)
    end
    direction = refract(rayin.direction, normal, nint)
    if !ismissing(direction)
        reflect_prob = schlick(cosine, d.refidx)
    else
        reflect_prob = 1.0
    end
    if rand() < reflect_prob
        Scatter(Ray(rec.p, reflected), Vec3(1.0, 1.0, 1.0), false)
    else
        Scatter(Ray(rec.p, direction), Vec3(1.0, 1.0, 1.0), true)
    end
end
