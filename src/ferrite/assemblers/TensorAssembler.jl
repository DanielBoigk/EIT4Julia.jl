function calculate_bilinear_map(a::AbstractVector, b::AbstractVector, fe::FerriteFESpace, M)
    cellvalues = fe.cellvalues
    dh = fe.dh
    n = fe.n
    rhs = zeros(n)
    n_basefuncs = getnbasefunctions(cellvalues)
    qpoints = getnquadpoints(cellvalues)
    re = zeros(n_basefuncs)
    
    for cell in CellIterator(dh)
        dofs = celldofs(cell)
        reinit!(cellvalues, cell)
        fill!(re, 0.0)
        
        ae = a[dofs]
        be = b[dofs]
        
        for q in 1:qpoints
            dΩ = getdetJdV(cellvalues, q)
            
            ∇a_q = zero(Vec{2,Float64})
            ∇b_q = zero(Vec{2,Float64})
            
            for j in 1:n_basefuncs
                ∇ϕⱼ = shape_gradient(cellvalues, q, j)
                ∇a_q += ae[j] * ∇ϕⱼ
                ∇b_q += be[j] * ∇ϕⱼ
            end
            
            grad_dot_product = ∇a_q ⋅ ∇b_q
            
            for i in 1:n_basefuncs
                ϕᵢ = shape_value(cellvalues, q, i)
                re[i] += grad_dot_product * ϕᵢ * dΩ
            end
        end
        assemble!(rhs, dofs, re)
    end
    
    return M \ rhs
end