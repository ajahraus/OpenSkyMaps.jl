module OpenSkyMaps
using Meshes
using LinearAlgebra: dot

export OpenSkyMap, polarFromCartesian, project_to_polar


struct OpenSkyMap
    data::BitMatrix
    OpenSkyMap(data::BitMatrix) = new(data)
    OpenSkyMap(size::Integer) = OpenSkyMap(BitMatrix(undef, (size, size)))
    OpenSkyMap(rows::Integer, cols::Integer) = OpenSkyMap(BitMatrix(undef, (rows, cols)))
    OpenSkyMap() = new(BitMatrix(hypot(indx[1], indx[2]) < 512 for indx in CartesianIndices((-512:511, -512:511))))
end

function polarFromCartesian(point)
    all(iszero, point.coords[1:2]) && return Point(0, 0)
    ρ = hypot(point.coords...)
    φ = acos(point.coords[3] / ρ)
    return Point(point.coords[1:2] .* (2φ / (π * ρ))...)
end

function filter_normals(geometry)
    normals = geometry |>
              elements .|>
              simplexify .|>
              first .|>
              normal

    sample_coordinates = geometry .|>
                         vertices .|>
                         first .|>
                         coordinates
    facing_bit_vector = dot.(normals, sample_coordinates) .< 0.0
    return facing_bit_vector
end

function filter_height(geometry)
    below_horizon = geometry .|>
                    vertices .|>
                    (v -> coordinates.(v)) .|>
                    (v -> getindex.(v, 3)) .|>
                    splat(max) .> 0.0
    return below_horizon
end

function collect_valid_faces(geometry, observer)
    geometry_t = geometry |>
                 Translate(-observer.coords...)
    filter_bv = filter_height(geometry_t) .&& filter_normals(geometry_t)
    [f for (f, b) in zip(elements(geometry_t), filter_bv) if b]
end

function clamp_z(point)
    c = coordinates(point)
    Point(c[1], c[2], max(c[3], 0.0))
end

function clamp_vertices_to_horizon(remaining_faces)
    faces_needing_clamping = remaining_faces .|>
                             vertices .|>
                             (v -> coordinates.(v)) .|>
                             (v -> getindex.(v, 3)) .|>
                             splat(min) .< 0.0
    clamped_faces = [b ? Ngon(clamp_z.(vertices(f))...) : f for (f, b) in zip(remaining_faces, faces_needing_clamping)]
    return clamped_faces
end

function satellite_unit_vector(φ, θ)
    z, r = sincos(φ)
    y, x = r .* sincos(θ)
    Point3f(x, y, z)
end

function project_to_polar(observer, geometry)
    remaining_faces = collect_valid_faces(geometry, observer)
    clamped_faces = clamp_vertices_to_horizon(remaining_faces)

    sphere_faces = [Ngon(polarFromCartesian.(vertices(f))...) for f in clamped_faces]
    return sphere_faces
end

end # Module

