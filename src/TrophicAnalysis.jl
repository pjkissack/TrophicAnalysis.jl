module TrophicAnalysis
    """
    TROPHIC ANALYSIS

    Adapted from Bazil Sansom's Trophic Analysis Toolbox python code 

    This Julia module [1] is for conducting trophic analysis as introduced by MacKay, Johnson and Sansom in [2], and also provides tools developed for network visualisation using these methods.
        Equation numbers refered to in annotation are from [2].
        
        If you make use of code provided here please cite [1,2].
        

    Functions included are:
    - adapted from trophic_tools.py:
        - trophic_levels : returns trophic levels as per [2]
        - trophic_incoherence : returns trophic incoherence as per [2] (where trophic coherence is 1-incoherence)
        - HHD: Calculate the network Helmholtz-Hodge decomposition
        - trophic_layout : returns a layout where y-possitions are given by trophic levels, and x-possitions based on a modified force-directed graph drawing algorithm to spread nodes out on the x-axis. For reproducibility, user can save and specify seed.
        - trophic_plot : plots network according to trophic_layout
    - additions:
        - trophic_levels_levine : returns trophic levels as per [3]; using either iterative or linear algebra method.
        - calculate_trophic_levels_iter : calculates trophic levels as per [3]; using iterative method.
        - calculate_trophic_levels_linear : calculates trophic levels as per [3]; using linear algebra method.

    REFS
     [1] Kissack, P., Drysdale, C. and Johnson, S. (202x).
     [2] MacKay, R. S., Johnson, S., Sansom, B. (2020). How directed is a directed network?. R Soc Open Sci. 7 (9): 201138. https://doi.org/10.1098/rsos.201138
     [3] Levine, S. (1980). Several measures of trophic structure applicable to complex food webs. J Theor Biol. 83 (2): 195-207. https://doi.org/10.1016/0022-5193(80)90288-X.
    """

    using Graphs
    using SimpleWeightedGraphs
    using LinearAlgebra
    using SparseArrays
    using LinearSolve
    using Random
    using GraphPlot
    using Plots
    using StatsBase

    function trophic_levels(G)
        """
        This function takes networkx graph object as input, and calculates trophic levels
            INPUTS:
                G networkx graph object
            OUTPUTS:
                h array of trophic levels
        """

        # Checks:
        if !isa(G,SimpleDiGraph) && !isa(G,SimpleWeightedDiGraph)
            throw(ArgumentError("This function takes graphs.jl SimpleDiGraph object or SimpleWeightedGraphs.jl SimpleWeightedDiGraph object."))
        end

        if !is_weakly_connected(G)
            throw(ArgumentError("Network must be weakly connected."))
        end

        # Check network weakly connected
        W = adjacency_matrix(G)

        in_weight = vec(sum(W,dims=1))
        out_weight = vec(sum(W,dims=2))
        u = float(in_weight + out_weight)
        v = float(in_weight - out_weight)

        L_init = Symmetric(diagm(u) - (W + W'))
        L = Matrix(L_init)
        L[1,1] = 0                              # what is this doing here?????

        h = solve(LinearProblem(L,v))
        return round.(h .- minimum(h),digits=10)

    end

    function trophic_incoherence(G)
        # This function takes networkx graph object as input, and returns 
        # trophic levels and coherence
        #     INPUTS
        #       G networkx graph object
        #     OUTPUTS:
        #       F_0 trophic incoherence
        #       h   array of trophic levels

        h = trophic_levels(G)
        W = adjacency_matrix(G)
        
        h_diff = h' .- h
        H = (h_diff .- 1).^2
        F_0 = sum(W.*H)/sum(W)

        return F_0, h
    end

    function HHD(G)
        """
        This function takes networkx graph object as input, and returns 
        network HHD: the HHD enables us to break down the flow on a directed network into 
        two flow components: potential flow (also directonal or gradient flow) and 
        circular flow. The potential flow between a pair of nodes is given by the difference 
        of their potentials obtained by the Helmholtz-Hodge decomposition [ref] or 
        equivalently trophic-levels (shown in [1]) and runs from a node with higher potential 
        to a node with lower potential. On the other hand, the circular flow component is 
        balanced so contributes nothing to node imbalance vector w_in-w_out (Eq.2.3 (the 
        difference between the flow into and out of each node)).
        
        INPUTS
          G networkx graph object
        OUTPUTS:
          F_c circular flow network
          F_p potential flow network
        """

        # Get trophic levels
        h = trophic_levels(G)

        # Obtain flow decomposition
        W = adjacency_matrix(G)
        h_diff = h' .- h

        F_p = W.*h_diff
        F_c = W - F_p

        # Directed network with positive sign ordered pairs
        for i in 1:size(W,1)
            for j in 1:size(W,2)
                if F_p[i,j] < 0
                    F_p[j,i] += abs(F_p[i,j])
                    F_p[i,j] = 0.0
                end
                if F_c[i,j] < 0
                    F_c[j,i] += abs(F_c[i,j])
                    F_c[i,j] = 0.0
                end
            end
        end

        return F_c, F_p
    end

    function trophic_layout(G;k=nothing,
                            ypos = nothing,
                            iterations = 50,
                            seed = nothing,
                            threshold = 1e-4,
                            SCC = false)
        """
        This function position nodes in network G using modified Fruchterman-Reingold layout. 
            The layout spreads the nodes on the x-axis taking y-possitions as given. 
            By default the function uses for y-possitions trophic levels as defined in [1], 
            but may also be passed any y-possitions user chooses to define.
            
            REQUIRED INPUTS:
            G   : networkx graph object.
                 Positions will be assigned to every node in G.
            
            OPTIONAL INPUTS
            k   : integer or None (default=None). If None the distance is set to 1/sqrt(nnodes) 
                where nnodes is the number of nodes.  Increase this value to spread nodes farther apart on x-axis .
            ypos: array or None (default=None). Initial y-positions for nodes. If None, then use
                trophic levels as defined by [1]. Alternatively user can pass any desired possitions as ypos. 
            iterations : integer (default=50)
                Maximum number of iterations taken.
            seed : integer, tuple, or None (default=None). For reproducible layouts, specify seed for random number generator.
                   This can be an integer; or the output seedState obtained from a previous run of trophic_layout or trophic_plot in order to reporduce
                   the same layout result.
                   If you run:
                       pos1, seedState = trophic_layout(G)
                       pos2, _ = trophic_layout(G,seed=seedState)
                    then pos1==pos2
            threshold: float (default = 1e-4)
                Threshold for relative error in node position changes.
                The iteration stops if the error is below this threshold.
            
            OUTPUTS
                pos a dictionary of possitions for each node in G
                seedState : (tuple) the seedState needed to reproduce layout obtained (e.g. if you run:
                                pos1, seedState = trophic_layout(G)
                                pos2, _ = trophic_layout(G,seed=seedState)
                          then pos2==pos1
        """

        # RANDOM STATE - unfinished
        if isnothing(seed) || isa(seed,Int) || isa(seed, Xoshiro)
            Random.seed!(seed)
            seedState = Random.GLOBAL_SEED
        else
            throw("Seed should be None (default); integer or Xoshiro (use existing RNG which is output of trophic_layout).")
        end
        
        if !isa(G,SimpleDiGraph) && !isa(G,SimpleWeightedDiGraph)
            throw(ArgumentError("This function takes graphs.jl SimpleDiGraph object or SimpleWeightedGraphs.jl SimpleWeightedDiGraph object."))
        end

        if !is_weakly_connected(G)
            throw(ArgumentError("Network must be weakly connected."))
        end

        A = adjacency_matrix(G)
        nnodes = nv(G)

        A_sym = A + A'

        if isnothing(ypos)
            h = trophic_levels(G)
            pos = hcat(rand(nnodes), h)
        else
            pos = hcat(rand(nnodes), ypos)
        end

        if isnothing(k)
            k = sqrt(1.0/nnodes)
        end

        t = maximum([(maximum(pos[:,1]) - minimum(pos[:,1])), (maximum(pos[:,2]) - minimum(pos[:,2]))]) * 0.1

        dt = t/(iterations + 1.0)
        
        delta_pos = zeros(nnodes, 2)
        displacement = zeros(nnodes, 2)

        for iteration in 1:iterations

            @inbounds for i in 1:nnodes
                for j in 1:nnodes
                    if i != j
                        dx = pos[i,1] - pos[j,1]
                        dy = pos[i,2] - pos[j,2]

                        distance = sqrt(dx^2 + dy^2)
                        distance = max(distance,0.01)

                        force = (k * k / distance^2 - A_sym[i,j] * distance / k)

                        displacement[i,1] += dx*force
                        displacement[i,2] += dy*force
                    end
                end
            end


            @inbounds for i in 1:nnodes
                length = sqrt(displacement[i,1]^2 + displacement[i,2]^2)
                length = max(length,0.01)

                delta_pos[i,1] = displacement[i,1] * t/length
                delta_pos[i,2] = 0.0
            end


            pos .+= delta_pos

            t -= dt

            err = norm(delta_pos) / nnodes

            if err < threshold
                break
            end
        end

        # pos_end = Dict(vertices(G) .=> (pos[:,1],pos[:,2]))
        if SCC
            SCCs = strongly_connected_components(G)
            n_SCCs = length(SCCs)

            for (i,SCC) in enumerate(SCCs)
                pos_SCC = copy(pos[:,1][SCC])

                pos_SCC .-= minimum(pos_SCC)
                if length(SCC) > 1
                    pos_SCC ./= (maximum(pos_SCC))
                end

                pos_SCC .* 100

                for (j,v) in enumerate(SCC)
                    pos[v,1] = (pos_SCC[j]/(n_SCCs*1.1)  + ((i-1)/n_SCCs))
                end
            end
        end

        return vertices(G), pos[:,1], .-pos[:,2], seedState
    end
        

    function trophic_plot(G;
                            k=nothing,
                            ypos=nothing,
                            iterations=50,
                            seed=nothing,
                            threshold=1e-4,
                            title="",
                            node_color=[0 0.451 0.7412])
        """
        This is just a wrapper for trophic_layout that automates some plotting decissions.
            Has same options and defaults as trophic_layout.
            Will add some plotting options.
            
            REQUIRED INPUTS:
            
            G   : networkx graph object.
                 Positions will be assigned to every node in G.
            
            OPTIONAL INPUTS
            
            Layout options:
            k   : integer or None (default=None). If None the distance is set to 1/sqrt(nnodes) 
                where nnodes is the number of nodes.  Increase this value to spread nodes farther apart on x-axis .
            ypos: array or None (default=None). Initial y-positions for nodes. If None, then use
                trophic levels as defined by [1]. Alternatively user can pass any desired possitions as ypos. 
            iterations : integer (default=50)
                Maximum number of iterations taken.
            seed : integer, tuple, or None (default=None). For reproducible layouts, specify seed for random number generator.
                   This can be an integer; or the output seedState obtained from a previous run of trophic_layout or trophic_plot in order to reporduce
                   the same layout result.
                   If you run:
                       pos1, seedState = trophic_layout(G)
                       pos2, _ = trophic_layout(G,seed=seedState)
                    then pos1==pos2
            threshold: float (default = 1e-4)
                Threshold for relative error in node position changes.
                The iteration stops if the error is below this threshold.
            
            Plotting options (a selection of draw_networkx options):
                title : (str) otionally provide a title for the chart
                node_color : color string, or array of floats, (default RGB triplet [0, 0.451, 0.7412])
                             Node color. Can be a single color format string, or a  sequence of colors with
                             the same length as nodelist. If numeric values are specified they will be mapped 
                             to colors using the cmap and vmin,vmax parameters.  See matplotlib.scatter for 
                             more details.
            
            OUTPUTS
                plots input network G
                seedState : (tuple) the seedState needed to reproduce layout obtained (e.g. if you run:
                                seedState = trophic_plot(G)
                                then this exact plot can be replicated by running
                                _ = trophic_plot(G,seed=seedState)
        """

        F_0, _ = trophic_incoherence(G)
        verts, xpos, ypos = trophic_layout(G,k=k,ypos=ypos,iterations=iterations,threshold=threshold,seed=seed)

        nnodes = nv(G)

        scale = 1/nnodes

        #pyplot()
        
        plot(xaxis=nothing)
        display(gplot(G,xpos,ypos))
        



    end

    
    ## END OF ORIGINAL TROPHIC TOOLBOX.


    ## Additional functionality:
    # Implementation of trophic levels definition from Levine, 1980;
    # i.e. as in python's networkx.trophic_levels.

    # Refs:
    #  [1] Levine, S. (1980). "Several measures of trophic structure applicable to complex food webs". Journal of Theoretical Biology. 83 (2): 195–207.
    

    function trophic_levels_levine(G::DiGraph; method::String="iterative")
        """
        Wrapper function to calculate trophic levels using specified method.
            INPUTS:
                G       : graphs.jl DiGraph object
                method  : String specifying method to use; options are:
                            "iterative" : iterative method
                            "linear"    : linear algebra method
            OUTPUTS:
                trophic_levels : array of trophic levels
        """
        if method == "iterative"
            return calculate_trophic_levels_iter(G)
        elseif method == "linear"
            return calculate_trophic_levels_linear(G)
        else
            error("Unknown method: $method")
        end
    end

    function calculate_trophic_levels_iter(G::DiGraph)
        """
        This function calculates trophic levels as per Levine, 1980;
            using iterative method.
            INPUTS:
                G : graphs.jl DiGraph object
            OUTPUTS:
                trophic_levels : array of trophic levels
        """
        N = nv(G)
        trophic_levels = ones(Float64, N)
        
        # Iteratively calculate trophic levels for non-basal nodes
        changed = true
        while changed
            changed = false
            for v in vertices(G)
                indeg = indegree(G, v)
                if indeg > 0
                    prey_levels = [trophic_levels[u] for u in inneighbors(G, v) if trophic_levels[u] > 0]
                    if length(prey_levels) == indeg
                        new_level = 1.0 + mean(prey_levels)
                        if trophic_levels[v] != new_level
                            trophic_levels[v] = new_level
                            changed = true
                        end
                    end
                end
            end
        end
        
        return trophic_levels
    end

    function calculate_trophic_levels_linear(G::DiGraph)
        """
        This function calculates trophic levels as per Levine, 1980;
            using linear algebra method.
            INPUTS:
                G : graphs.jl DiGraph object
            OUTPUTS:
                trophic_levels : array of trophic levels
        """
        N = nv(G)
        A = adjacency_matrix(G)
                
        # Create laplacian-like matrix 
        Λ = zeros(Float64, N, N)
        v = ones(Float64, N)
        for i in 1:N
            if indegree(G, i) > 0
                for j in 1:N
                    if A[j, i] == 1
                        Λ[i, j] = - 1.0
                    end
                end
                Λ[i, i] = indegree(G, i)
                v[i] = indegree(G, i)
            else
                Λ[i, i] = 1.0
            end
        end

        prob = LinearProblem(Λ, v)
        linsolve = init(prob)
        trophic_levels = solve!(linsolve).u
        
        return trophic_levels
    end

end
