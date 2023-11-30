module OpenSkyMaps
using Meshes
using LinearAlgebra: dot
using Base.Iterators: take
using Base.Iterators: map as imap

export OpenSkyMap, polarFromCartesian, project_to_polar, filter_height, filter_normal, sample_normal, sample_coordinate, clamp_z


struct OpenSkyMap
    data::BitMatrix
    OpenSkyMap(data::BitMatrix) = new(data)
    OpenSkyMap(size::Integer) = OpenSkyMap(BitMatrix(undef, (size, size)))
    OpenSkyMap(rows::Integer, cols::Integer) = OpenSkyMap(BitMatrix(undef, (rows, cols)))
    OpenSkyMap() = new(BitMatrix(hypot(indx[1], indx[2]) < 512 for indx in CartesianIndices((-512:511, -512:511))))
end

function polarFromCartesian(point::Point3f)::Point2f
    all(iszero, point.coords[1:2]) && return Point(0, 0)
    ρ = hypot(point.coords...)
    φ = acos(point.coords[3] / ρ)
    return Point2f(point.coords[1:2] .* (2φ / (π * ρ))...)
end

sample_normal = normal ∘ first ∘ simplexify
sample_coordinate = coordinates ∘ first

filter_normal(el, vert) = sample_normal(el) ⋅ sample_coordinate(vert) < 0.0

filter_height(v::NTuple{3,Point3f}) = any(c[3] > 0.0f0 for c in coordinates.(v))

function collect_valid_faces(geometry::Mesh)
    els = elements(geometry)
    verts = vertices.(geometry)
    filter_bv = @. filter_height(verts) && filter_normal(els, verts)
    [f for (f, b) in zip(els, filter_bv) if b]
end

function clamp_z(point::Point3f)::Point3f
    c = coordinates(point)
    c[3] ≥ 0.0 ? point : Point3f(c[1], c[2], 0.0)
end

face_needs_clamp(face) = any(coordinates(v)[3] < 0.0 for v in vertices.(face))

clamp_vertices_to_horizon(f) = face_needs_clamp(f) ? Ngon(clamp_z.(vertices(f))...) : f

function project_to_polar(observer::Point3f, geometry::Mesh)
    geometry = geometry |> Translate(-observer.coords...)
    remaining_faces = collect_valid_faces(geometry)
    clamped_faces = clamp_vertices_to_horizon.(remaining_faces)
    sphere_faces = [Ngon(polarFromCartesian.(vertices(f))...) for f in clamped_faces]
    return sphere_faces
end

end # Module

