struct Lambert <: Material
    albedo::Vec3
end

function scatter(l::Lambert, rayin::Ray, rec::HitRecord)::Scatter
    target = rec.p .+ rec.normal .+ sphere_rand()
    Scatter(Ray(rec.p, target .- rec.p), l.albedo, true)
end
