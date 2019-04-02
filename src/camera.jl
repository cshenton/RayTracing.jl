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

function Camera()
    from = Vec3(13.0, 2.0, 3.0)
    to = Vec3(0.0, 0.0, 0.0)
    focus = 10.0
    aperture = 0.1
    Camera(from, to, Vec3(0.0, 1.0, 0.0), 20.0, 2.0, aperture, focus)
end

function Ray(c::Camera, u::Float64, v::Float64)
    rd = c.radius .* sphere_rand()
    offset = u .* rd[1] .+ v .* rd[2]
    Ray(c.origin .+ offset, c.lower_left_corner .+ u .* c.horizontal .+ v .* c.vertical .- c.origin)
end
