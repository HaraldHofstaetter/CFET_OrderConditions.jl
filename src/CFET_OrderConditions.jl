module CFET_OrderConditions

using Giac
using Combinatorics

export commutator, LCCC
export gen_exp_tBt_derivatives_coeffs
export gen_C_derivatives_at_0
export gen_C_derivatives_at_0_in_terms_of_A
export gen_exp_tBt_derivatives_at_0_in_terms_of_A
export expand_commutators
export gen_CFET_order_conditions

export number_commutators_of_grade_equal_to
export number_commutators_of_grade_equal_or_less_than

export integer_vector_weighted, lyndon_words, gen_hall_basis_elements_of_grade
export gen_CFET_order_conditions_without_redundancies

typealias Commutator Array{Int64,1}
typealias LCCC{T} Dict{Array{Array{Int64,1},1},T} # Linear Combination of Commutator Chains



function gen_exp_tBt_derivatives_coeffs(N::Integer)
    A = [Dict{Array{Int64,1},Int64}([]=>1),    #Id
         Dict{Array{Int64,1},Int64}([0]=>1)]
    for n = 2:N
        a0 = A[n]
        # a1 = a0*C:
        a1 = [vcat(key, 0)=>val for (key,val) in a0]
        # a1 = a1 + (d/dt)a0:
        for (p,c) in a0
            for j = 1:length(p)
                q = copy(p)
                q[j] += 1
                d = get(a1, q, 0)
                a1[q] = d+c
            end            
        end                
        push!(A, a1)
    end
    A
end


function gen_C_derivatives_at_0(N::Integer)
    A = Array{Dict{Array{Int64,1},Int64},1}[[Dict(vcat(zeros(Int64,k),1)=>1)] for k=1:N-1]
    for k = 1:N-1        
        for l = 1:N-k-1
            a0 = A[k][l]
            a1 = Dict{Array{Int64,1},Int64}()
            for (p,c) in a0
                for j=1:length(p)
                    q = copy(p)
                    q[j] += 1
                    if q[end] != q[end-1]
                        d = get(a1, q, 0)
                        a1[q] = d+c
                    end    
                end            
            end                           
            push!(A[k], a1)
        end
    end
    C = [Dict{Array{Int64,1},Int64}([n]=>(n+1)) for n=0:N]
    for n = 2:N
        a1 = C[n+1]
        for k=1:n-1
            f = binomial(n,k+1)
            a0 = A[k][n-k-1+1]
            for (p,c) in a0
                d = get(a1, p, 0//1)
                a1[p] = d+f*c
            end                                   
        end        
        C[n+1] = a1
    end    
    C
end


function gen_C_derivatives_at_0_in_terms_of_A{T}(N::Integer, b::Array{T,1})
    C_derivatives_at_0=gen_C_derivatives_at_0(N)
    C1 = [Dict{Array{Int64,1},T}() for n=0:N]
    for n = 0:N
        for (p,coeff) in C_derivatives_at_0[n+1]
            C1[n+1][p] = coeff * prod([b[p[j]+1]  for j =1:length(p)])
        end
    end
    C1
end

Id(T::Type; f=one(T)) = LCCC{T}(Array[Int64[]]=>f)
is_c_number(x) = length(x)==1 && length(first(keys(x))[1])==0

function mult{T}(x::LCCC{T}, y::Dict{Array{Int64,1},T})
    r = LCCC{T}()
    for (p,c) in x
        for (q, d) in y
            pq = copy(p)
            push!(pq, q)
            r[pq] = c*d
        end
    end     
    r
end


function mult{T}(x::LCCC{T}, y::LCCC{T})
    if is_c_number(x)
        c = first(values(x))
        return LCCC{T}([p=>c*d for (p,d) in y])
    elseif is_c_number(y)
        c = first(values(y))
        return LCCC{T}([p=>c*d for (p,d) in x])
    else
        r = LCCC{T}()
        for (p,c) in x
            for (q, d) in y
                pq = vcat(p,q)
                r[pq] = c*d
            end
        end     
        return r
    end
end


function add{T}(x::LCCC{T}, y::LCCC{T}; 
                fx=one(T), fy=one(T))
    r = copy(x)
    for (p,c) in y
        d = get(r, p, 0)
        r[p] = fy*c + fx*d
    end
    r
end

function gen_exp_tBt_derivatives_at_0_in_terms_of_A{T}(N::Integer, b::Array{T,1}) 
    C_derivatives_at_0_in_terms_of_A=gen_C_derivatives_at_0_in_terms_of_A(N-1, b)
    exp_tBt_derivatives_coeffs=gen_exp_tBt_derivatives_coeffs(N)
    r = [LCCC{T}(Array[Int64[]]=>one(T))] # Id    
    for n = 1:N
        y = LCCC{T}()
        for (p,d) in exp_tBt_derivatives_coeffs[n+1]
            x = LCCC{T}(
            [Array[p1]=>d*d1 for (p1,d1) in C_derivatives_at_0_in_terms_of_A[p[1]+1]])
            for i = 2:length(p)
                x = mult(x, C_derivatives_at_0_in_terms_of_A[p[i]+1])
            end  
            y = add(y, x)
        end
        push!(r,y)
    end
    r
end


function expand_commutators{T}(a::Array{Int64,1}, f::T; coeff_type=T)
    if length(a)>=3
        b = LCCC{coeff_type}(Array[[a[1]]]=>one(T))
        c = expand_commutators(a[2:end], f, coeff_type=coeff_type)
        return add(mult(b,c), mult(c,b), fy=-one(T))
    elseif length(a)==2
        return LCCC{coeff_type}(Array[[a[1]], [a[2]]]=>f, Array[[a[2]], [a[1]]]=>-f)
    else
        return LCCC{coeff_type}(Array[a]=>f)
    end
end


function expand_commutators{T}(a::LCCC{T})    
    b = LCCC{T}()
    for (p,c) in a        
        q = expand_commutators(p[1], one(T), coeff_type=T)
        for i = 2:length(p)
            q = mult(q, expand_commutators(p[i], one(T), coeff_type=T))
        end
        b = add(b, q, fy=c)
    end
    b
end


multinomial_coeff(q::Int, k::Array{Int,1}) = div(factorial(q), prod([factorial(i) for i in k]))


function gen_CFET_order_conditions(N::Integer,J::Integer)
    b = [giac[giac(string("b",j,"_",n)) for n=0:N] for j=1:J]
    C = gen_exp_tBt_derivatives_at_0_in_terms_of_A(N, b[1])
    C = [expand_commutators(c) for c in C]
    C = [LCCC{giac}([key=>factor(val) for (key,val) in c]) for c in C]
    CC = Array[C]
    for j=2:J
        push!(CC, [LCCC{giac}([key=>subst(val,b[1],b[j]) for (key,val) in c]) for c in C])
    end    
    Y = Dict{Array{Int64,1},giac}[]
    for q=0:N-1
        y = LCCC{giac}()
        for k0 in Combinatorics.WithReplacementCombinations(1:J,q+1)
            k = zeros(Int,J)
            for j in k0
                k[j]+=1
            end
            coeff = multinomial_coeff(q+1, k) 
            x = Id(giac, f=giac(coeff))
            for j = J:-1:1
                x = mult(x, CC[j][k[j]+1])
            end
            y = add(y, x)
        end   
        for k0 in Combinatorics.WithReplacementCombinations(1:J+1,q)
            k = zeros(Int,J+1)
            for j in k0
                k[j]+=1
            end
            coeff = multinomial_coeff(q, k) 
            x = LCCC{giac}(Array[[k[J+1]]]=>coeff)
            for j = J:-1:1
                x = mult(x, CC[j][k[j]+1])
            end
            y = add(y, x, fy=-one(giac))
            
        end   
        push!(Y, Dict{Array{Int64,1},giac}([[k[1] for k in key]=>factor(val) for (key,val) in y]))
    end
    Y
end

function number_commutators_of_grade_equal_to(n::Integer)
    if n==1
        return 1
    end   
    f = factor(n)
    p = collect(keys(f))
    k = collect(values(f))
    nu = 0
    for l = 0:(2^length(p)-1)
        lb = bin(l,length(p))
        mu = (-1)^count_ones(l)
        pot = prod([p[j]^(k[j]-(lb[j]=='0'?0:1)) for j=1:length(p)])
        nu += mu*2^pot    
    end
    div(nu,n)
end

number_commutators_of_grade_equal_or_less_than(n::Integer) =
    sum([num_of_eqs(k) for k=1:n])



# The following functions are stolen from MuPAD-Combinat,
# see http://mupad-combinat.sourceforge.net

# mupad/COMBINAT/words.mu/evaluationTable:
function evaluation_table(word::Vector)
    t = Dict{Any, Int}()
    for l in word
        if haskey(t, l)
            t[l] += 1
        else
            t[l] = 1
        end
    end
    t
end 

# mupad/COMBINAT/words.mu/standard
function standard(word::Vector)
    t = evaluation_table(word)
    kk = sort(collect(keys(t)))
    offset = 0
    tmp = 0
    for k in kk
        tmp = t[k]
        t[k] = offset
        offset += tmp
    end
    r = zeros(Int, length(word))
    i = 1
    for l in word
        t[l] += 1
        r[i] = t[l]
        i += 1
    end
    r
end

# mupad/COMBINAT/integerVectorsWeighted.mu/
# (There it is a local function inside function list
function recfun(n::Int, l1::Vector{Int})
    result = []
    l = copy(l1)
    w = pop!(l)
    if length(l)==0
        d = div(n, w)
        if mod(n, w)==0
            return Vector{Int}[[d]]
        else
            return Vector{Int}[]
        end
    end
    for d = div(n, w):-1:0
        result = vcat(result, [vcat(x,d) for x in recfun(n-d*w, l)])
    end
    result
end

# mupad/COMBINAT/permutations.mu/mult2
function mult2(x::Vector{Int}, y::Vector{Int})
    nx = length(x)
    ny = length(y)
    # complete with fixed points if the permutations are not
    # of the same length
    x = vcat(x, nx+1:ny)
    y = vcat(y, ny+1:nx)
    result = x
    for i=1:length(x)
        result[y[i]] = x[i]
    end
    result
end


# mupad/COMBINAT/integerVectorsWeighted.mu/list
function integer_vector_weighted(n::Int, l::Vector{Int})
    if length(l)==0
        if n==0
            return Vector{Int}[[]]
        else
            return Vector{Int}[]
        end
    end

    # Heuristic: we sort the weights increasingly to limit
    # the backtracking, and reorder back the result afterward
    perm = standard(l)
    sort!(l)

    [mult2(x, perm) for x in recfun(n,l)]
end



# For the following 2 functions, see: 
#   Sawada, Joe.  "A fast algorithm to generate necklaces with fixed content", 
#   Source Theoretical Computer Science archive Volume 301 , Issue 1-3 (May 2003)
# and sage/src/sage/combinat/necklace.py


function simple_fixed_content1(res::Vector{Vector{Int}}, a::Vector{Int}, 
            content::Vector{Int}, t::Int, p::Int, k::Int)
    n = length(a)
    if t>n
        if n==p 
            push!(res, a)
        end
    else
        for j=a[t-p]:k-1
            if content[j+1]!=0
                a[t] = j
                content[j+1] -= 1
                if j==a[t-p]
                    simple_fixed_content1(res, copy(a), content, t+1, p, k) 
                else
                    simple_fixed_content1(res, copy(a), content, t+1, t, k) 
                end
                content[j+1] += 1
            end
        end
    end
end

function simple_fixed_content(content::Vector{Int})
    a = zeros(Int, sum(content))
    content[1] -= 1
    k = length(content)
    res = Vector{Int}[]
    simple_fixed_content1(res, a, content, 2, 1, k)
    res
end

#The following function is stolen from sage/combinat/necklace.py

function lyndon_words(e::Vector{Int})
    k=1
    while e[k]==0
        k += 1
    end
    Vector{Int}[[i+k for i in z] for z in simple_fixed_content(e[k:end])]
end

#The following is stolen from http://mathoverflow.net/questions/97703/list-of-hall-basis

function gen_hall_basis_elements_of_grade(g::Int)
    L = integer_vector_weighted(g, collect(1:g))
    res = Vector{Int}[]
    for l in L
        for w in lyndon_words(l)
            push!(res, w-1)
        end
    end
    res
end

function gen_CFET_order_conditions_without_redundancies(N::Integer,J::Integer)
    C0 = gen_CFET_order_conditions(N, J)
    Dict{Array{Int, 1}, giac}[Dict{Array{Int, 1}, giac}(
        [b=>C0[k][b] for b in gen_hall_basis_elements_of_grade(k)]) for k=1:length(C0)]
end


end # module
