# Neighborhood Strength Driven Label Propagation Algorithm
function nsdlpa{V,T<:Real}(g::AbstractGraph{V}, λ::T=1.0, CN::Vector{Int}=triangle(g))
    @graph_requires g edge_map incidence_list vertex_map

    # calculate number of common neighbors
    N = num_vertices(g)
    membership = [1:N]
    label_counters = zeros(T, N)
    dominant_labels = Int[]
    nonzero_labels = Int[]

    running = true
    STEP = 0
    while running
        STEP += 1
        running = false

        # Shuffle the node ordering vector
        X = shuffle(collect(vertices(g)))

        # In the prescribed order, loop over the vertices and reassign labels
        for v in X
            v_idx = vertex_index(v, g)

            # Clear dominant_labels and nonzero_labels
            dominant_labels = Int[]
            nonzero_labels = Int[]

            # recount
            max_count = zero(T)

            for w in out_edges(v, g)
                k = membership[vertex_index(target(w, g), g)]

                # skip if it has no label yet
                k != 0 || continue

                was_zero = label_counters[k] == zero(T)
                label_counters[k] += 1 + λ*CN[edge_index(w, g)]
                if was_zero && label_counters[k] != zero(T)
                    # counter just became nonzero
                    push!(nonzero_labels, k)
                end
                if max_count < label_counters[k]
                    max_count = label_counters[k]
                    resize!(dominant_labels, 1)
                    dominant_labels[1] = k
                elseif max_count == label_counters[k]
                    push!(dominant_labels, k)
                end
            end

            if !isempty(dominant_labels)

                # Select randomly from the dominant labels
                k = dominant_labels[rand(1:length(dominant_labels))]

                # Check if the current label of the node is also dominant
                if label_counters[membership[v_idx]] != max_count
                    # Nope, we need at least one more iteration
                    running = true
                end

                # Update label of the current node
                membership[v_idx] = k
            end
            # Clear the nonzero elements in label_counters
            for i in nonzero_labels
                label_counters[i] = zero(T)
            end
        end
    end
    permute_labels!(membership), STEP
end

function nsdlpa1{V,T<:Real}(g::AbstractGraph{V}, λ::T=1.0, CN::Vector{Int}=triangle(g), DS::Vector{Int}=quadrangle(g))
    @graph_requires g edge_map incidence_list vertex_map

    # calculate number of common neighbors
    N = num_vertices(g)
    membership = [1:N]
    label_counters = zeros(T, N)
    proximity_sum = zeros(T, N)

    running = true
    STEP = 0
    while running
        STEP += 1
        running = false

        # Shuffle the node ordering vector
        X = shuffle(collect(vertices(g)))

        # In the prescribed order, loop over the vertices and reassign labels
        for v in X
            v_idx = vertex_index(v, g)

            # Clear dominant_labels and nonzero_labels
            dominant_labels = Int[]
            nonzero_labels = Int[]
            similarity = Float64[]

            # recount
            max_count = zero(T)

            for w in out_edges(v, g)
                k = membership[vertex_index(target(w, g), g)]

                # skip if it has no label yet
                k != 0 || continue

                was_zero = label_counters[k] == zero(T)
                label_counters[k] += 1 + λ*CN[edge_index(w, g)]
                proximity_sum[k] += DS[edge_index(w, g)]
                if was_zero && label_counters[k] != zero(T)
                    # counter just became nonzero
                    push!(nonzero_labels, k)
                end
                if max_count < label_counters[k]
                    max_count = label_counters[k]
                    resize!(dominant_labels, 1)
                    resize!(similarity, 1)
                    dominant_labels[1] = k
                    similarity[1] = proximity_sum[k]
                elseif max_count == label_counters[k]
                    push!(dominant_labels, k)
                    push!(similarity, proximity_sum[k])
                end
            end

            if !isempty(dominant_labels)

                max_similarity = maximum(similarity)
                dominant_labels = dominant_labels[similarity.<=max_similarity]

                # Select randomly from the dominant labels
                k = dominant_labels[rand(1:length(dominant_labels))]

                # Check if the current label of the node is also dominant
                if label_counters[membership[v_idx]] != max_count
                    # Nope, we need at least one more iteration
                    running = true
                end

                # Update label of the current node
                membership[v_idx] = k
            end
            # Clear the nonzero elements in label_counters
            for i in nonzero_labels
                label_counters[i] = zero(T)
                proximity_sum[i] = zero(T)
            end
        end
    end
    permute_labels!(membership), STEP
end
