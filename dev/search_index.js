var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ClimFlowsPlots","category":"page"},{"location":"#ClimFlowsPlots","page":"Home","title":"ClimFlowsPlots","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ClimFlowsPlots.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ClimFlowsPlots, ClimFlowsPlots.VoronoiSphere, ClimFlowsPlots.SphericalInterpolations]","category":"page"},{"location":"#ClimFlowsPlots.VoronoiSphere.plot_2D-Tuple{Any, Observables.Observable}","page":"Home","title":"ClimFlowsPlots.VoronoiSphere.plot_2D","text":"data = Makie.Observable(data) # must be on primal mesh\nfig = plot_2D(sphere, data ; resolution=1.0)\ndisplay(fig)\ndata[] = new_data # updates the plot\ndisplay(fig)\n\nPlots scalar data as a function of longitude and latitude. Data must be a scalar defined at primal cells of a spherical Voronoi mesh and passed as a Makie.Observable. Data is interpolated linearly to a regular lon-lat grid of resolution resolution before plotting.\n\n\n\n\n\n","category":"method"},{"location":"#ClimFlowsPlots.VoronoiSphere.plot_native_3D-Tuple{Any, Observables.Observable}","page":"Home","title":"ClimFlowsPlots.VoronoiSphere.plot_native_3D","text":"data = Makie.Observable(data) # must be on primal mesh\nfig = plot_native_3D(sphere, data ; zoom=1.6)\ndisplay(fig)\ndata[] = new_data # updates the plot\ndisplay(fig)\n\nPlots scalar data in a stereographic view. Data must be a scalar defined at primal cells of a spherical Voronoi mesh and passed as a Makie.Observable. Data on\nnative mesh is used without interpolation, but plotting interpolates linearly between\nprimal cell centers (which coincides with triangle vertices).\n\n\n\n\n\n","category":"method"},{"location":"#ClimFlowsPlots.VoronoiSphere.plot_orthographic-Tuple{Any, Observables.Observable}","page":"Home","title":"ClimFlowsPlots.VoronoiSphere.plot_orthographic","text":"data = Makie.Observable(data) # must be on primal mesh\nfig = plot_voronoi_orthographic(sphere, data ; zoom=1.6)\ndisplay(fig)\ndata[] = new_data # updates the plot\ndisplay(fig)\n\nPlots scalar data in a stereographic view. Data must be a scalar defined at primal cells of a spherical Voronoi mesh and passed as a Makie.Observable. Data is interpolated linearly to a regular lon-lat grid of resolution resolution before plotting.\n\n\n\n\n\n","category":"method"},{"location":"#ClimFlowsPlots.SphericalInterpolations","page":"Home","title":"ClimFlowsPlots.SphericalInterpolations","text":"Acceptably fast computation of interpolation weights from an unstructured spherical triangulation to arbitrary points on the unit sphere. The mesh is first recursively decomposed into submeshes using Metis, resulting in a tree. Bounding spheres are constructed starting from the leaves. The tree structure is used to search for the triangle containing a given point in logarithmic time. Interpolation itsef is linear with positive weights.\n\n\n\n\n\n","category":"module"},{"location":"#ClimFlowsPlots.SphericalInterpolations.SubGraph","page":"Home","title":"ClimFlowsPlots.SphericalInterpolations.SubGraph","text":"Subset of a parent graph. Keeps the mapping from indices in the subgraph to indices in the parent graph.\n\n\n\n\n\n","category":"type"},{"location":"#ClimFlowsPlots.SphericalInterpolations.Tree","page":"Home","title":"ClimFlowsPlots.SphericalInterpolations.Tree","text":"Tree with data of type Node\n\n\n\n\n\n","category":"type"},{"location":"#ClimFlowsPlots.SphericalInterpolations.lonlat_interp-Tuple{Any, Any, Any}","page":"Home","title":"ClimFlowsPlots.SphericalInterpolations.lonlat_interp","text":"Returns a function that interpolates data given at mesh nodes onto a regular lon-lat mesh. Lons and Lats are ine degrees.\n\n\n\n\n\n","category":"method"},{"location":"#ClimFlowsPlots.SphericalInterpolations.meshtree-Tuple{Metis.Graph, Any}","page":"Home","title":"ClimFlowsPlots.SphericalInterpolations.meshtree","text":"Turns a graph into a tree by recursive partitioning into parts subgraphs.\n\n\n\n\n\n","category":"method"},{"location":"#ClimFlowsPlots.SphericalInterpolations.spheretree","page":"Home","title":"ClimFlowsPlots.SphericalInterpolations.spheretree","text":"Given bounding spheres for leaves, constructs bounding spheres at each level of the tree.\n\n\n\n\n\n","category":"function"},{"location":"#ClimFlowsPlots.SphericalInterpolations.subgraphs-Tuple{ClimFlowsPlots.SphericalInterpolations.SubGraph, Any}","page":"Home","title":"ClimFlowsPlots.SphericalInterpolations.subgraphs","text":"Partitions a (sub)graph into parts subgraphs using Metis.\n\n\n\n\n\n","category":"method"},{"location":"#ClimFlowsPlots.SphericalInterpolations.traverse-Union{Tuple{Fun}, Tuple{Fun, Any, Any}} where Fun","page":"Home","title":"ClimFlowsPlots.SphericalInterpolations.traverse","text":"Applies function fun to x, v for all triangles v of tree tree whose bounding sphere contains x     until fun returns anything different from nothing.\n\n\n\n\n\n","category":"method"}]
}
