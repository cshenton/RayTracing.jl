struct Metal <: Material
    albedo::Vec3
    fuzz::Float64
    Metal(albedo, fuzz) = new(albedo, fuzz < 1 ? fuzz : 1.0)
end

function scatter(m::Metal, rayin::Ray, rec::HitRecord)::Scatter
    reflected = reflect(unit_vector(rayin.direction), rec.normal)
    Scatter(Ray(rec.p, reflected .+ m.fuzz .* sphere_rand()), m.albedo, (reflected ⋅ rec.normal) > 0.0)
end
