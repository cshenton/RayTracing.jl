module RayTracing

using Distributed
using LinearAlgebra
using StaticArrays


Option{T} = Union{Missing, T}
Vec3 = SVector{3}

abstract type Hitable end
abstract type Material end

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


struct Ray
    origin::Vec3
    direction::Vec3
end

function point(r::Ray, t::Float64)::Vec3
    r.origin .+ t .* r.direction
end

struct HitRecord
    t::Float64
    p::Vec3
    normal::Vec3
    mat::Material
end


struct Sphere <: Hitable
    center::Vec3
    radius::Float64
    mat::Material
end

function hit(s::Sphere, r::Ray, tmin::Float64, tmax::Float64)::Option{HitRecord}
    oc = r.origin - s.center
    a = (r.direction ⋅ r.direction)
    b = (oc ⋅ r.direction)
    c = (oc ⋅ oc) - s.radius^2
    discriminant = b^2 - a*c
    if discriminant > 0
        temp = -(b + sqrt(discriminant)) / a
        if tmin < temp < tmax
            t = temp
            p = point(r, t)
            n = (p - s.center) / s.radius
            return HitRecord(t, p, n, s.mat)
        end
        temp = (-b + sqrt(discriminant)) / a
        if tmin < temp < tmax
            t = temp
            p = point(r, t)
            n = (p - s.center) / s.radius
            return HitRecord(t, p, n, s.mat)
        end
    end
    return missing
end


struct HitableList <: Hitable
    list::Vector{Hitable}
end

function hit(h::HitableList, r::Ray, tmin::Float64, tmax::Float64)::Option{HitRecord}
    closest = tmax
    rec = missing
    for el in h.list
        temprec = hit(el, r, tmin, closest)
        if !ismissing(temprec)
            rec = temprec
            closest = rec.t
        end
    end
    rec
end


struct Scatter
    ray::Ray
    attenuation::Vec3
    reflect::Bool
end


struct Lambert <: Material
    albedo::Vec3
end

function scatter(l::Lambert, rayin::Ray, rec::HitRecord)::Scatter
    target = rec.p .+ rec.normal .+ sphere_rand()
    Scatter(Ray(rec.p, target .- rec.p), l.albedo, true)
end


struct Metal <: Material
    albedo::Vec3
    fuzz::Float64
    Metal(albedo, fuzz) = new(albedo, fuzz < 1 ? fuzz : 1.0)
end

function scatter(m::Metal, rayin::Ray, rec::HitRecord)::Scatter
    reflected = reflect(unit_vector(rayin.direction), rec.normal)
    Scatter(Ray(rec.p, reflected .+ m.fuzz .* sphere_rand()), m.albedo, (reflected ⋅ rec.normal) > 0.0)
end


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


struct Camera
    lower_left_corner::Vec3
    horizontal::Vec3
    vertical::Vec3
    origin::Vec3
    u::Vec3
    v::Vec3
    w::Vec3
    radius::Float64
end

function Camera(from::Vec3, to::Vec3, vup::Vec3, vfov::Float64, aspect::Float64, aperture::Float64, focus::Float64)
    theta = vfov * π / 180
    half_height = tan(theta/2)
    half_width = aspect * half_height
    w = unit_vector(from-to)
    u = unit_vector(cross(vup, w))
    v = cross(w, u)
    Camera(
        from .- half_width.*u .- half_height.*v .- w,
        2 .* half_width .* u,
        2 .* half_height .* v,
        from, u, v, w,
        aperture / 2,
    )
end

function Ray(c::Camera, u::Float64, v::Float64)
    rd = c.radius .* sphere_rand()
    offset = u .* rd[1] .+ v .* rd[2]
    Ray(c.origin .+ offset, c.lower_left_corner .+ u .* c.horizontal .+ v .* c.vertical .- c.origin)
end



function color(r::Ray, world::Hitable, depth::Int)::Vec3
    rec = hit(world, r, 0.0, typemax(Float64))
    if !ismissing(rec)
        s = scatter(rec.mat, r, rec)
        if s.reflect && depth < 20
            return s.attenuation .* color(s.ray, world, depth+1)
        else
            return Vec3(0.0, 0.0, 0.0)
        end
    else
        unit_direction = Vec3(r.direction)
        t = 0.5 * (unit_direction[2] + 1.0)
        (1.0 - t) .* Vec3(1.0, 1.0, 1.0) .+ t.*Vec3(0.5, 0.7, 1.0)
    end
end

function scene()::Hitable
    n = 500
    spheres = Sphere[]
    push!(spheres, Sphere(Vec3(0.0, -1000.0, 0.0), 1000, Lambert(Vec3(0.5, 0.5, 0.5))))
    for a=-11:11
        for b=-11:11
            choose_mat = rand()
            center = Vec3(a+0.9*rand(), 0.2, b+0.9*rand())
            if norm(center .- Vec3(4.0, 0.2, 0.0)) > 0.9
                if choose_mat < 0.8
                    push!(spheres, Sphere(center, 0.2, Lambert(Vec3(rand()*rand(), rand()*rand(), rand()*rand()))))
                elseif choose_mat < 0.95
                    push!(spheres, Sphere(center, 0.2, Metal(Vec3(0.5*(1+rand()), 0.5*(1+rand()), 0.5*(1+rand())), 0.5*rand())))
                else
                    push!(spheres, Sphere(center, 0.2, Dielectric(1.5)))
                end
            end
        end
    end

    push!(spheres, Sphere(Vec3(0.0, 1.0, 0.0), 1.0, Metal(Vec3(0.7, 0.6, 0.5), 0.0)))
    push!(spheres, Sphere(Vec3(-4.0, 1.0, 0.0), 1.0, Lambert(Vec3(0.4, 0.2, 0.1))))
    push!(spheres, Sphere(Vec3(4.0, 1.0, 0.0), 1.0, Dielectric(1.5)))
    HitableList(spheres)
end


function main()
    nx = 400
    ny = 200
    ns = 200

    println("P3")
    println("$nx $ny")
    println("255")

    world = scene()

    from = Vec3(13.0, 2.0, 3.0)
    to = Vec3(0.0, 0.0, 0.0)
    focus = 10.0
    aperture = 0.1
    cam = Camera(from, to, Vec3(0.0, 1.0, 0.0), 20.0, nx/ny, aperture, focus)

    for j in reverse(1:ny)
        for i in 1:nx
            cols = Vector{Vec3}(undef, ns)
            for k in 1:ns
                r = Ray(cam, (i+rand())/nx, (j+rand())/ny)
                cols[k] = color(r, world, 0)
            end
            col = sum(cols) ./ ns
            col = sqrt.(col)
            ir = trunc(Int, 255.99*col[1])
            ig = trunc(Int, 255.99*col[2])
            ib = trunc(Int, 255.99*col[3])
            println("$ir $ig $ib")
        end
    end
end

main()

end # module
