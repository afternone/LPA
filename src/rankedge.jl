function similarity{V}(g::AbstractGraph{V})
    similarities = zeros(num_edges(g))
    for e in edges(g)
        e_idx = edge_index(e, g)
        u = source(e, g)
        v = target(e, g)
        u_nei = out_neighbors(u, g)
        v_nei = out_neighbors(v, g)
        uv_common_nei = intersect(u_nei, v_nei)
        uv_all_nei = union(u_nei, v_nei)
        similarities[e_idx] = (length(uv_common_nei) + 2) / length(uv_all_nei)
    end
    similarities
end

function triangle{V}(g::AbstractGraph{V})
    △ = zeros(Int, num_edges(g))
    for e in edges(g)
        e_idx = edge_index(e, g)
        u = source(e, g)
        v = target(e, g)
        u_nei = out_neighbors(u, g)
        v_nei = out_neighbors(v, g)
        uv_nei = intersect(u_nei, v_nei)
        △[e_idx] = length(uv_nei)
    end
    △
end

function quadrangle{V}(g::AbstractGraph{V})
    □ = zeros(Int, num_edges(g))
    for e in edges(g)
        e_idx = edge_index(e, g)
        u = source(e, g)
        v = target(e, g)
        u_nei = out_neighbors(u, g)
        v_nei = out_neighbors(v, g)
        u_unique_nei = setdiff(u_nei, v_nei)
        u_unique_nei = setdiff(u_unique_nei, v)
        v_unique_nei = setdiff(v_nei, u_nei)
        v_unique_nei = setdiff(v_unique_nei, u)
        for w in u_unique_nei
            w_nei = out_neighbors(w, g)
            vw_nei = intersect(v_unique_nei, w_nei)
            □[e_idx] += length(vw_nei)
        end
    end
    □
end

function triquadrsim{V}(g::AbstractGraph{V})
    △ = zeros(Int, num_edges(g))
    □ = zeros(Int, num_edges(g))
    similarities = zeros(num_edges(g))
    for e in edges(g)
        e_idx = edge_index(e, g)
        u = source(e, g)
        v = target(e, g)
        u_nei = out_neighbors(u, g)
        v_nei = out_neighbors(v, g)
        uv_nei = intersect(u_nei, v_nei)
        uv_all_nei = union(u_nei, v_nei)
        similarities[e_idx] = (length(uv_nei) + 2) / length(uv_all_nei)
        △[e_idx] = length(uv_nei)
        u_unique_nei = setdiff(u_nei, uv_nei)
        u_unique_nei = setdiff(u_unique_nei, v)
        v_unique_nei = setdiff(v_nei, uv_nei)
        v_unique_nei = setdiff(v_unique_nei, u)
        for w in u_unique_nei
            w_nei = out_neighbors(w, g)
            vw_nei = intersect(v_unique_nei, w_nei)
            □[e_idx] += length(vw_nei)
        end
    end
    △, □, similarities
end
