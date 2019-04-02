struct Ray
    origin::Vec3
    direction::Vec3
end

function point(r::Ray, t::Float64)::Vec3
    r.origin .+ t .* r.direction
end
