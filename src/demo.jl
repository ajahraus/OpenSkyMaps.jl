includet("OpenSkyMaps.jl")
using .OpenSkyMaps
using Meshes
import GLMakie

function main()
    observer = Point3f(15.0f0, -5.0f0, 0.0f0)

    geometry1 = Box((-0.5f0, 0.0f0, 0.0f0), (0.5f0, 1.0f0, 10.0f0)) |>
                boundary |>
                (g -> refine(g, TriRefinement())) |>
                (g -> refine(g, TriRefinement()))

    geometries = [geometry1 |> Translate((rand(Float32) + 3.0f0 * x, rand(Float32) + -10.0f0 * y, 0.0f0)) for x in 0:10 for y in 0:1]

    polar_faces = [convexhull(project_to_polar(observer, geometry)) for geometry in geometries]
    geometry = reduce(merge, geometries)

    if false
        fig = GLMakie.Figure(resolution=(1800, 900))
        viz(fig[1, 1], Sphere(observer, 0.1f0), color=:red, aspect=1)
        viz!(fig[1, 1], geometry, showfacets=true)


        viz(fig[1, 2], Sphere((0.0f0, 0.0f0), 1.0f0), color=:red, aspect=1)
        viz!(fig[1, 2], Sphere((0.0f0, 0.0f0), 0.01f0))
        viz!(fig[1, 2], reduce(merge, polar_faces), showfacets=true)
        display(fig)
    end
    return polar_faces
end
main()