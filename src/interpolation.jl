"""
Acceptably fast computation of interpolation weights from an unstructured spherical triangulation to arbitrary points on the unit sphere.
The mesh is first recursively decomposed into submeshes using Metis, resulting in a tree. Bounding spheres are constructed starting from the leaves.
The tree structure is used to search for the triangle containing a given point in logarithmic time.
Interpolation itsef is linear with positive weights.
"""
module SphericalInterpolations

macro fast(code)
    debug = haskey(ENV, "GF_DEBUG") && (ENV["GF_DEBUG"]!="")
    return debug ? esc(code) : esc(quote @inbounds $code end)
end

using Metis, LightGraphs, BoundingSphere, StaticArrays

point(lon,lat) = SVector( map(Float64, (cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat))) )

"""
    interp = lonlat_interp(mesh, lons, lats)
    data_lonlat = interp(data_mesh)
Return `interp`, which can be called to interpolate data
given at mesh nodes onto a regular lon-lat mesh.
*Lons and Lats are in degrees*.

"""
function lonlat_interp(mesh, lons, lats)
    ijrange = eachindex(mesh.Ai)
    x = [ point((pi/180)*lon,(pi/180)*lat) for lat in lats, lon in lons]
    pts = [ point(mesh.lon_i[ij], mesh.lat_i[ij]) for ij in ijrange]
    g = SimpleGraph([ Edge(mesh.edge_down_up[1,i], mesh.edge_down_up[2,i]) for i in eachindex(mesh.le)])
    sphtree = spheretree(meshtree(Metis.graph(g),4), dual_spheres(mesh))
    ww = [ traverse(xx, sphtree) do x, v
        return weights(x, v, pts, mesh.dual_vertex)
    end for xx in x ]
    return Interpolator(ijrange, ww)
#    return data -> interpolate(ww, data)
end

struct Interpolator{IJ, W}
    ijrange::IJ
    ww::W
end
function (interp::Interpolator)(data::AbstractMatrix)
    if axes(data,1) == interp.ijrange
        return interpolate_HV(interp.ww, data)
    else
        @assert axes(data,2) == interp.ijrange
        return interpolate_VH(interp.ww, data)
    end
end
function (interp::Interpolator)(data::AbstractVector)
    @assert axes(data, 1) == interp.ijrange
    return interpolate(interp.ww, data)
end

# Linearly interpolated value from indices (i,j,k) with weights (a,b,c)
@inline interpolate( ((i,j,k), (a,b,c)), data::AbstractVector) = @fast a*data[i]+b*data[j]+c*data[k]
# Linearly interpolated values from a matrix of tuples of indices (i,j,k) and weights (a,b,c)
interpolate(w::Matrix, data::AbstractVector) = [interpolate(ww, data) for ww in w]
interpolate(w::Matrix, datas::Tuple) = Tuple(interpolate(w,data) for data in datas)

function interpolate_VH(w::Matrix, data::AbstractMatrix)
    result = similar(data, size(data,1), size(w, 1), size(w,2))
    @fast for j in axes(w,1), i in axes(w, 2)
        (ii,jj,kk), (a,b,c) = w[j,i]
        @simd for k in axes(data,1)
            result[k,j,i] = muladd(a , data[k,ii], muladd(b, data[k,jj], c*data[k,kk]))
        end
    end
    return result
end

# Returns interpolation weights if point `x` inside triangle `( pts[vtx[i,v]] for i=1:3 )`, else nothing
function weights(x, v, pts, vtx)
    i, j, k = ( vtx[i,v] for i in 1:3 )
    ijk = i, j, k # tuple
    a, b, c = ( pts[i]   for i in ijk )
    xab = triprod(x,a,b)
    xbc = triprod(x,b,c)
    xca = triprod(x,c,a)
    if xab*xbc>0 && xab*xca>0 # all triple products have same sign => point is inside triangle
        w = inv(xab+xbc+xca)
        return ijk, (w*xbc, w*xca, w*xab)
    else
        return nothing
    end
end

@inline triprod(p,q,r) = dprod(p, vprod(q,r))
@inline dprod((x,y,z),(a,b,c)) = @fastmath a*x+b*y+c*z
@inline vprod((x,y,z),(a,b,c)) = @fastmath (y*c-z*b, z*a-x*c, x*b-y*a)

"""
Tree with data of type Node
"""
struct Tree{Node}
    node :: Node
    child :: Vector{Tree{Node}}
end

"""
Subset of a parent graph. Keeps the mapping from indices in the subgraph
to indices in the parent graph."""
struct SubGraph
    glob  :: Vector{Int32} # index of each cell in parent graph
    graph :: Metis.Graph
end

"""
Applies function `fun` to `x, v` for all triangles `v` of tree `tree` whose bounding sphere contains `x`
    until `fun` returns anything different from nothing.
"""
function traverse(fun::Fun, x, tree) where Fun
    node = tree.node
    center, radius = node.center, node.radius
    result = nothing
    if isinsphere(x, center, radius)
        child=tree.child
        if length(child)>1
            for c in child
                result = traverse(fun, x, c)
                result === nothing || break
            end
        else
            for v in node.glob
                result = fun(x, v)
                result === nothing || break
            end
        end
    end
    return result
end
@inline isinsphere((x,y,z),(a,b,c), radius) = @fastmath (x-a)^2+(y-b)^2+(z-c)^2<=radius*radius*1.000001

"""
Given bounding spheres for leaves, constructs bounding spheres at each level of the tree.
"""
function spheretree(sg::Tree, spheres, glob=eachindex(spheres))
    glob = glob[sg.node.glob]
    if length(sg.child)>1
        # recursively compute bounding spheres of children
        child = [ spheretree(c, spheres, glob) for c in sg.child ]
        center, radius = bounding_sphere(c.node for c in child)
        return Tree( (glob=Int[], center=center, radius=radius), child)
    else
        # This is a leaf, use user-provided bounding spheres
        center, radius = bounding_sphere(spheres[v] for v in glob)
        return Node( (glob=glob, center=center, radius=radius) )
    end
end
Node(x::T) where T = Tree{T}(x,Tree{T}[])

# Find sphere containg a set of spheres. The algorithm is not sharp.
# We find a sphere containing all sphere centers, then add the radius of the largest sphere.
# TODO : reduce bounding sphere, to the sphere containing its intersection with the unit sphere
# see https://mathworld.wolfram.com/Sphere-SphereIntersection.html
@inline function bounding_sphere(spheres)
    center, radius = BoundingSphere.boundingsphere([ sphere.center for sphere in spheres ])
    return center, radius + maximum(sphere.radius for sphere in spheres)
end

function dual_sphere(pts, v)
    pts = [pts[vv] for vv in v]
    center, radius = BoundingSphere.boundingsphere(pts)
    return (center=center, radius=radius)
end

function dual_spheres(mesh)
    pts = [ point(mesh.lon_i[ij], mesh.lat_i[ij]) for ij in eachindex(mesh.Ai) ]
    return [dual_sphere(pts, @view mesh.dual_vertex[:,v]) for v in eachindex(mesh.Av)]
end

"""
Turns a graph into a tree by recursive partitioning into `parts` subgraphs.
"""
meshtree(graph::Metis.Graph, parts) = meshtree(SubGraph(1:graph.nvtxs, graph), parts)

function meshtree(sg::SubGraph, parts)
    graph=sg.graph
    if graph.nvtxs>parts*parts
        return Tree(sg, [meshtree(sub, parts) for sub in subgraphs(graph, parts)])
    else
        return Tree(sg, Tree{SubGraph}[])
    end
end

"""
Partitions a (sub)graph into `parts` subgraphs using Metis.
"""
subgraphs(sg::SubGraph, npart) = subgraphs(sg.graph, npart)

function subgraphs(graph::Metis.Graph, npart)
    owner = Metis.partition(graph, npart)
    return [subgraph(p, owner, graph) for p=1:npart]
end

function subgraph(part, owner, graph)
    xadj=similar(graph.xadj)
    adjncy=similar(graph.adjncy)
    glob=similar(graph.xadj)

    ncell=0
    nedge=0
    for i in eachindex(owner)
        if owner[i]==part
            ncell = ncell+1
            glob[ncell] = i
            xadj[ncell] = nedge+1
            # now filter those neighbors that belong to the same partition
            for edge = graph.xadj[i]:graph.xadj[i+1]-1
                j = graph.adjncy[edge]
                if owner[j]==part
                    nedge = nedge+1
                    adjncy[nedge] = j
                end
            end
        end
    end
    xadj[ncell+1] = nedge+1

    # translate global index to local index
    loc = Dict( (glob[i], Int32(i)) for i=1:ncell )
    adjncy = [ loc[i] for i in adjncy[1:nedge] ]

    return SubGraph(glob[1:ncell], Metis.Graph(ncell, xadj[1:ncell+1], adjncy))
end

end # module
