module TrophicTuning
    """
    TROPHIC TUNING
    
    Generate directed networks with a specified target trophic incoherence,
    using Preferential Preying model [4] and Generalised Preferential Preying model [5].

    If making use of this code please cite [1].

    REFS:
     [1] Kissack, P., Drysdale, C. and Johnson, S. (202x).
     [4] Johnson, S., Domínguez-García, V., Donetti, L. & Muñoz, M. A. (2014). Trophic coherence determines food-web stability. Proc. Natl Acad. Sci. USA 111, 17923–17928. https://doi.org/10.1073/pnas.1409077111
     [5] Klaise, J., Johnson, S. (2016). From neurons to epidemics: How trophic coherence affects spreading processes. Chaos 26 (6): 065310. https://doi.org/10.1063/1.4953160    
    """
    using Graphs
    using StatsBase

    include("./TrophicAnalysis.jl")
    using .TrophicAnalysis

    function trophic_tuned_network(B::Int64, 
                                    N::Int64,
                                    L::Int64, 
                                    T::Float64; 
                                    algorithm::String="gen_pref_preying",
                                    weakly_connected::Bool=true,
                                    rseed::Union{Int64, Nothing} = nothing)

        if !isnothing(rseed)
            Random.seed!(rseed)
        end

        if algorithm == "pref_preying"
            while true
                G = pref_preying(B, N, L, T)
                if !weakly_connected || is_weakly_connected(G)
                    return G
                end
            end
        elseif algorithm == "gen_pref_preying"
            while true
                G = gen_pref_preying(B, N, L, T)
                if !weakly_connected || is_weakly_connected(G)
                    return G
                end
            end
        else
            error("Unknown algorithm: $algorithm")
        end
    end

    function pref_preying(B::Int64, 
                            N::Int64, 
                            L::Int64, 
                            T::Float64)

                        
        # Input validation
        if B >= N
            error("Number of basal nodes B must be less than total nodes N.")
        end
        if L >= N*(N-1)
            error("Number of links L must be less than total possible directed links.")
        end
        if T <= 0
            error("Trophic specialisation parameter T must be positive.")
        end
        if B <= 0 || N <= 0 || L <= 0
            error("Number of basal nodes, total nodes, and links must be positive integers.")
        end
        if L < (N - B)
            error("Number of links L must be at least equal to the number of non-basal nodes.")
        end

        # initialize directed graph
        G = DiGraph(N)

        # beta distribution parameters for number of prey links
        α = 1
        β = (N^2 - B^2)/(2*L) - 1

        temp_trophic_levels = zeros(Float64,N)

        # add trophic levels of basal nodes
        for i in 1:B
            temp_trophic_levels[i] = 1.0
        end

        # add initial links and temp trophic levels for non-basal nodes
        for i in B+1:N
            candidate_prey = collect(1:i-1)
            first_prey = rand(candidate_prey)
            add_edge!(G, first_prey, i)

            # update temp trophic level
            temp_trophic_levels[i] = 1.0 + temp_trophic_levels[first_prey]

            # get probabilities for additional prey
            P = exp.(-abs.(temp_trophic_levels[i] .- temp_trophic_levels[candidate_prey])./T)

            κ = round(Int,rand(Beta(α,β))*length(candidate_prey))

            additional_prey = sample(candidate_prey, Weights(P), κ; replace=false)

            for prey in additional_prey
                add_edge!(G, prey, i)
            end
        end
        
        return G
    end

    function gen_pref_preying(B::Int64, 
                                N::Int64, 
                                L::Int64, 
                                T::Float64)
        
        # Input validation
        if B >= N
            error("Number of basal nodes B must be less than total nodes N.")
        end
        if L >= N*(N-1)
            error("Number of links L must be less than total possible directed links.")
        end
        if T <= 0
            error("Trophic specialisation parameter T must be positive.")
        end
        if B <= 0 || N <= 0 || L <= 0
            error("Number of basal nodes, total nodes, and links must be positive integers.")
        end
        if L < (N - B)
            error("Number of links L must be at least equal to the number of non-basal nodes.")
        end

        # initialize directed graph
        G = DiGraph(N)

        # add initial links for non-basal nodes
        for i in B+1:N
            candidate_prey = collect(1:i-1)
            first_prey = rand(candidate_prey)
            add_edge!(G, first_prey, i)
        end

        # get preliminary trophic levels
        trophic_prelim = TrophicAnalysis.trophic_levels_levine(G)
        trophic_array = zeros(Float64, 1,N)
        for i in 1:N
            trophic_array[1,i] = trophic_prelim[i]
        end

        # compute tentative trophic distances
        xij = trophic_array .- trophic_array'

        # compute connection probabilities
        P = exp.(- (xij .- 1.0).^2 ./ (2.0 * T))

        # get candidate edges 
        candidate_edges = [(prey,predator) for prey in B+1:N
                                for predator in 1:N
                                if prey != predator && !has_edge(G,prey,predator)]

        # select additional edges
        n_additional_edges = L - (N - B)
        edge_probs = [P[prey,predator] for (prey,predator) in candidate_edges]
        selected_indices = sample(1:length(candidate_edges), 
                                    Weights(edge_probs), 
                                    n_additional_edges; 
                                    replace=false)

        for idx in selected_indices
            prey, predator = candidate_edges[idx]
            add_edge!(G, prey, predator)
        end

        return G
    end

end