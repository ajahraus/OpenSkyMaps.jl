includet("OpenSkyMaps.jl")
using .OpenSkyMaps
using Meshes
using BenchmarkTools


observer() = Point3f(rand(Float32) + 15.0f0, rand(Float32) + -5.0f0, rand(Float32) + 0.0f0)

geometry() = Box((-0.5f0, 0.0f0, 0.0f0), (0.5f0, 1.0f0, 10.0f0)) |>
             boundary |>
             (g -> refine(g, TriRefinement())) |>
             (g -> refine(g, TriRefinement()))

geometries() = [geometry() |> Translate((rand(Float32) + 3.0f0 * x, rand(Float32) + -10.0f0 * y, 0.0f0)) for x in 0:10 for y in 0:1]

polar_faces(observer, geometries) = [project_to_polar(observer, geometry) for geometry in geometries]

function main()
    o = observer()
    g = geometries()

    @benchmark polar_faces($o, $g)
end
