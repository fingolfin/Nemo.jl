################################################################################
#
#  fmpz_mod_mat.jl: flint fmpz_mod_mat (matrices over Z/nZ, large n)
#
################################################################################

################################################################################
#
#  Data type and parent object methods
#
################################################################################

dense_matrix_type(::Type{ZZModRingElem}) = ZZModMatrix

###############################################################################
#
#   Similar
#
###############################################################################

function similar(::MatElem, R::ZZModRing, r::Int, c::Int)
  z = ZZModMatrix(R, r, c)
  return z
end

################################################################################
#
#  Manipulation
#
################################################################################

@inline function getindex(a::T, i::Int, j::Int) where T <: Zmod_fmpz_mat
  @boundscheck _checkbounds(a, i, j)
  u = ZZRingElem()
  ccall((:fmpz_mod_mat_get_entry, libflint), Nothing,
        (Ref{ZZRingElem}, Ref{T}, Int, Int, Ref{fmpz_mod_ctx_struct}),
        u, a, i - 1 , j - 1, base_ring(a).ninv)
  return ZZModRingElem(u, base_ring(a)) # no reduction needed
end

# as above, but as a plain ZZRingElem, no bounds checking
function getindex_raw(a::T, i::Int, j::Int) where T <: Zmod_fmpz_mat
  u = ZZRingElem()
  ccall((:fmpz_mod_mat_get_entry, libflint), Nothing,
        (Ref{ZZRingElem}, Ref{T}, Int, Int, Ref{fmpz_mod_ctx_struct}),
        u, a, i - 1, j - 1, base_ring(a).ninv)
  return u
end

@inline function setindex!(a::T, u::ZZRingElem, i::Int, j::Int) where T <: Zmod_fmpz_mat
  @boundscheck _checkbounds(a, i, j)
  R = base_ring(a)
  setindex_raw!(a, _reduce(u, R.ninv), i, j)
end

@inline function setindex!(a::T, u::ZZModRingElem, i::Int, j::Int) where T <: Zmod_fmpz_mat
  @boundscheck _checkbounds(a, i, j)
  (base_ring(a) != parent(u)) && error("Parent objects must coincide")
  setindex_raw!(a, u.data, i, j) # no reduction needed
end

function setindex!(a::T, u::Integer, i::Int, j::Int) where T <: Zmod_fmpz_mat
  setindex!(a, ZZRingElem(u), i, j)
end

# as per setindex! but no reduction mod n and no bounds checking
@inline function setindex_raw!(a::T, u::ZZRingElem, i::Int, j::Int) where T <: Zmod_fmpz_mat
  ccall((:fmpz_mod_mat_set_entry, libflint), Nothing,
        (Ref{T}, Int, Int, Ref{ZZRingElem}, Ref{Nothing}),
        a, i - 1, j - 1, u, C_NULL) # ctx is not needed here
end

function setindex!(a::ZZModMatrix, b::ZZModMatrix, r::UnitRange{Int64}, c::UnitRange{Int64})
  _checkbounds(a, r, c)
  size(b) == (length(r), length(c)) || throw(DimensionMismatch("tried to assign a $(size(b, 1))x$(size(b, 2)) matrix to a $(length(r))x$(length(c)) destination"))
  A = view(a, r, c)
  ccall((:fmpz_mod_mat_set, libflint), Nothing,
        (Ref{ZZModMatrix}, Ref{ZZModMatrix}, Ref{Nothing}), A, b, C_NULL) # ctx not used
end

function deepcopy_internal(a::ZZModMatrix, dict::IdDict)
  z = ZZModMatrix(nrows(a), ncols(a), base_ring(a).ninv)
  if isdefined(a, :base_ring)
    z.base_ring = a.base_ring
  end
  ccall((:fmpz_mod_mat_set, libflint), Nothing,
        (Ref{ZZModMatrix}, Ref{ZZModMatrix}, Ref{fmpz_mod_ctx_struct}),
        z, a, base_ring(a).ninv)
  return z
end

number_of_rows(a::T) where T <: Zmod_fmpz_mat = a.r

number_of_columns(a::T) where T <: Zmod_fmpz_mat = a.c

base_ring(a::T) where T <: Zmod_fmpz_mat = a.base_ring

function one(a::ZZModMatrixSpace)
  (nrows(a) != ncols(a)) && error("Matrices must be square")
  z = a()
  ccall((:fmpz_mod_mat_one, libflint), Nothing,
        (Ref{ZZModMatrix}, Ref{fmpz_mod_ctx_struct}), z, base_ring(a).ninv)
  return z
end

function iszero(a::T) where T <: Zmod_fmpz_mat
  r = ccall((:fmpz_mod_mat_is_zero, libflint), Cint,
            (Ref{T}, Ref{fmpz_mod_ctx_struct}), a, base_ring(a).ninv)
  return Bool(r)
end

################################################################################
#
#  Comparison
#
################################################################################

==(a::T, b::T) where T <: Zmod_fmpz_mat = (a.base_ring == b.base_ring) &&
Bool(ccall((:fmpz_mod_mat_equal, libflint), Cint,
           (Ref{T}, Ref{T}, Ref{fmpz_mod_ctx_struct}),
           a, b, base_ring(a).ninv))

isequal(a::T, b::T) where T <: Zmod_fmpz_mat = ==(a, b)

################################################################################
#
#  Transpose
#
################################################################################

function transpose(a::T) where T <: Zmod_fmpz_mat
  z = similar(a, ncols(a), nrows(a))
  ccall((:fmpz_mod_mat_transpose, libflint), Nothing,
        (Ref{T}, Ref{T}, Ref{fmpz_mod_ctx_struct}), z, a, base_ring(a).ninv)
  return z
end

function transpose!(a::T) where T <: Zmod_fmpz_mat
  !is_square(a) && error("Matrix must be a square matrix")
  ccall((:fmpz_mod_mat_transpose, libflint), Nothing,
        (Ref{T}, Ref{T}, Ref{fmpz_mod_ctx_struct}), a, a, base_ring(a).ninv)
end

###############################################################################
#
#   Row and column swapping
#
###############################################################################

#= Not currently implemented in Flint

function swap_rows!(x::T, i::Int, j::Int) where T <: Zmod_fmpz_mat
ccall((:fmpz_mod_mat_swap_rows, libflint), Nothing,
(Ref{T}, Ptr{Nothing}, Int, Int), x, C_NULL, i - 1, j - 1)
return x
end

function swap_rows(x::T, i::Int, j::Int) where T <: Zmod_fmpz_mat
(1 <= i <= nrows(x) && 1 <= j <= nrows(x)) || throw(BoundsError())
y = deepcopy(x)
return swap_rows!(y, i, j)
end

function swap_cols!(x::T, i::Int, j::Int) where T <: Zmod_fmpz_mat
ccall((:fmpz_mod_mat_swap_cols, libflint), Nothing,
(Ref{T}, Ptr{Nothing}, Int, Int), x, C_NULL, i - 1, j - 1)
return x
end

function swap_cols(x::T, i::Int, j::Int) where T <: Zmod_fmpz_mat
(1 <= i <= ncols(x) && 1 <= j <= ncols(x)) || throw(BoundsError())
y = deepcopy(x)
return swap_cols!(y, i, j)
end

function reverse_rows!(x::T) where T <: Zmod_fmpz_mat
ccall((:fmpz_mod_mat_invert_rows, libflint), Nothing,
(Ref{T}, Ptr{Nothing}), x, C_NULL)
return x
end

reverse_rows(x::T) where T <: Zmod_fmpz_mat = reverse_rows!(deepcopy(x))

function reverse_cols!(x::T) where T <: Zmod_fmpz_mat
ccall((:fmpz_mod_mat_invert_cols, libflint), Nothing,
(Ref{T}, Ptr{Nothing}), x, C_NULL)
return x
end

reverse_cols(x::T) where T <: Zmod_fmpz_mat = reverse_cols!(deepcopy(x))

=#

################################################################################
#
#  Unary operators
#
################################################################################

function -(x::T) where T <: Zmod_fmpz_mat
  z = similar(x)
  ccall((:fmpz_mod_mat_neg, libflint), Nothing,
        (Ref{T}, Ref{T}, Ref{fmpz_mod_ctx_struct}), z, x, base_ring(x).ninv)
  return z
end

################################################################################
#
#  Binary operators
#
################################################################################

function +(x::T, y::T) where T <: Zmod_fmpz_mat
  check_parent(x,y)
  z = similar(x)
  ccall((:fmpz_mod_mat_add, libflint), Nothing,
        (Ref{T}, Ref{T}, Ref{T}, Ref{fmpz_mod_ctx_struct}),
        z, x, y, base_ring(x).ninv)
  return z
end

function -(x::T, y::T) where T <: Zmod_fmpz_mat
  check_parent(x,y)
  z = similar(x)
  ccall((:fmpz_mod_mat_sub, libflint), Nothing,
        (Ref{T}, Ref{T}, Ref{T}, Ref{fmpz_mod_ctx_struct}),
        z, x, y, base_ring(x).ninv)
  return z
end

function *(x::T, y::T) where T <: Zmod_fmpz_mat
  (base_ring(x) != base_ring(y)) && error("Base ring must be equal")
  (ncols(x) != nrows(y)) && error("Dimensions are wrong")
  z = similar(x, nrows(x), ncols(y))
  ccall((:fmpz_mod_mat_mul, libflint), Nothing,
        (Ref{T}, Ref{T}, Ref{T}, Ref{fmpz_mod_ctx_struct}),
        z, x, y, base_ring(x).ninv)
  return z
end


################################################################################
#
#  Unsafe operations
#
################################################################################

function mul!(a::T, b::T, c::T) where T <: Zmod_fmpz_mat
  ccall((:fmpz_mod_mat_mul, libflint), Nothing,
        (Ref{T}, Ref{T}, Ref{T}, Ref{fmpz_mod_ctx_struct}),
        a, b, c, base_ring(b).ninv)
  return a
end

function add!(a::T, b::T, c::T) where T <: Zmod_fmpz_mat
  ccall((:fmpz_mod_mat_add, libflint), Nothing,
        (Ref{T}, Ref{T}, Ref{T}, Ref{fmpz_mod_ctx_struct}),
        a, b, c, base_ring(b).ninv)
  return a
end

function zero!(a::T) where T <: Zmod_fmpz_mat
  ccall((:fmpz_mod_mat_zero, libflint), Nothing,
        (Ref{T}, Ref{fmpz_mod_ctx_struct}), a, base_ring(a).ninv)
  return a
end

function mul!(z::Vector{ZZRingElem}, a::T, b::Vector{ZZRingElem}) where T <: Zmod_fmpz_mat
  ccall((:fmpz_mod_mat_mul_fmpz_vec_ptr, libflint), Nothing,
        (Ptr{Ref{ZZRingElem}}, Ref{T}, Ptr{Ref{ZZRingElem}}, Int, Ref{fmpz_mod_ctx_struct}),
        z, a, b, length(b), base_ring(a).ninv)
  return z
end

function mul!(z::Vector{ZZRingElem}, a::Vector{ZZRingElem}, b::T) where T <: Zmod_fmpz_mat
  ccall((:fmpz_mod_mat_fmpz_vec_mul_ptr, libflint), Nothing,
        (Ptr{Ref{ZZRingElem}}, Ptr{Ref{ZZRingElem}}, Int, Ref{T}, Ref{fmpz_mod_ctx_struct}),
        z, a, length(a), b, base_ring(b).ninv)
  return z
end

function Generic.add_one!(a::ZZModMatrix, i::Int, j::Int)
  @boundscheck _checkbounds(a, i, j)
  GC.@preserve a begin
    x = mat_entry_ptr(a, i, j)
    ccall((:fmpz_mod_add_si, libflint), Nothing,
          (Ptr{ZZRingElem}, Ptr{ZZRingElem}, Int, Ref{fmpz_mod_ctx_struct}),
          x, x, 1, base_ring(a).ninv)
  end
  return a
end

################################################################################
#
#  Ad hoc binary operators
#
################################################################################

function *(x::T, y::Int) where T <: Zmod_fmpz_mat
  z = similar(x)
  ccall((:fmpz_mod_mat_scalar_mul_si, libflint), Nothing,
        (Ref{T}, Ref{T}, Int, Ref{fmpz_mod_ctx_struct}),
        z, x, y, base_ring(x).ninv)
  return z
end

*(x::Int, y::T) where T <: Zmod_fmpz_mat = y*x

function *(x::T, y::ZZRingElem) where T <: Zmod_fmpz_mat
  z = similar(x)
  ccall((:fmpz_mod_mat_scalar_mul_fmpz, libflint), Nothing,
        (Ref{T}, Ref{T}, Ref{ZZRingElem}, Ref{fmpz_mod_ctx_struct}),
        z, x, y, base_ring(x).ninv)
  return z
end

*(x::ZZRingElem, y::T) where T <: Zmod_fmpz_mat = y*x

*(x::T, y::Integer) where T <: Zmod_fmpz_mat = x*flintify(y)

*(x::Integer, y::T) where T <: Zmod_fmpz_mat = y*x

function *(x::ZZModMatrix, y::ZZModRingElem)
  (base_ring(x) != parent(y)) && error("Parent objects must coincide")
  return x*y.data
end

*(x::ZZModRingElem, y::ZZModMatrix) = y*x

################################################################################
#
#  Powering
#
################################################################################

#= Not implemented in Flint yet

function ^(x::T, y::Int) where T <: Zmod_fmpz_mat
if y < 0
x = inv(x)
y = -y
end
z = similar(x)
ccall((:fmpz_mod_mat_pow, libflint), Nothing,
(Ref{T}, Ref{T}, Int), z, x, y)
return z
end
=#

function ^(x::T, y::ZZRingElem) where T <: Zmod_fmpz_mat
  (y > ZZRingElem(typemax(Int))) &&
  error("Exponent must be smaller than ", ZZRingElem(typemax(Int)))
  (y < ZZRingElem(typemin(Int))) &&
  error("Exponent must be bigger than ", ZZRingElem(typemin(Int)))
  return x^(Int(y))
end

################################################################################
#
#  Strong echelon form and Howell form
#
################################################################################

function strong_echelon_form!(a::T) where T <: Zmod_fmpz_mat
  ccall((:fmpz_mod_mat_strong_echelon_form, libflint), Nothing,
        (Ref{T}, Ref{fmpz_mod_ctx_struct}), a, base_ring(a).ninv)
end

@doc raw"""
    strong_echelon_form(a::ZZModMatrix)

Return the strong echeleon form of $a$. The matrix $a$ must have at least as
many rows as columns.
"""
function strong_echelon_form(a::ZZModMatrix)
  (nrows(a) < ncols(a)) &&
  error("Matrix must have at least as many rows as columns")
  z = deepcopy(a)
  strong_echelon_form!(z)
  return z
end

function howell_form!(a::T) where T <: Zmod_fmpz_mat
  ccall((:fmpz_mod_mat_howell_form, libflint), Nothing,
        (Ref{T}, Ref{fmpz_mod_ctx_struct}), a, base_ring(a).ninv)
end

@doc raw"""
    howell_form(a::ZZModMatrix)

Return the Howell normal form of $a$. The matrix $a$ must have at least as
many rows as columns.
"""
function howell_form(a::ZZModMatrix)
  (nrows(a) < ncols(a)) &&
  error("Matrix must have at least as many rows as columns")

  z = deepcopy(a)
  howell_form!(z)
  return z
end

################################################################################
#
#  Trace
#
################################################################################

function tr(a::T) where T <: Zmod_fmpz_mat
  !is_square(a) && error("Matrix must be a square matrix")
  R = base_ring(a)
  r = ZZRingElem()
  ccall((:fmpz_mod_mat_trace, libflint), Nothing,
        (Ref{ZZRingElem}, Ref{T}, Ref{fmpz_mod_ctx_struct}),
        r, a, R.ninv)
  return elem_type(R)(r, R)
end

################################################################################
#
#  Determinant
#
################################################################################

#= Not implemented in Flint yet

function det(a::ZZModMatrix)
!is_square(a) && error("Matrix must be a square matrix")
if is_prime(a.n)
r = ccall((:fmpz_mod_mat_det, libflint), UInt, (Ref{ZZModMatrix}, ), a)
return base_ring(a)(r)
else
try
return AbstractAlgebra.det_fflu(a)
catch
return AbstractAlgebra.det_df(a)
end
end
end

=#

################################################################################
#
#  Rank
#
################################################################################

#= Not implemented in Flint yet

function rank(a::T) where T <: Zmod_fmpz_mat
r = ccall((:fmpz_mod_mat_rank, libflint), Int, (Ref{T}, ), a)
return r
end

=#

################################################################################
#
#  Inverse
#
################################################################################

function inv(a::ZZModMatrix)
  !is_square(a) && error("Matrix must be a square matrix")
  if is_probable_prime(modulus(base_ring(a)))
    X, d = pseudo_inv(a)
    if !is_unit(d)
      error("Matrix is not invertible")
    end
    return divexact(X, d)
  else
    b = map_entries(x -> x.data, a)
    c, d = pseudo_inv(b)
    R = base_ring(a)
    if !isone(gcd(d, modulus(R)))
      error("Matrix not invertible")
    end
    return change_base_ring(R, c) * inv(R(d))
  end
end

#= Not implemented in Flint yet

function inv(a::T) where T <: Zmod_fmpz_mat
!is_square(a) && error("Matrix must be a square matrix")
z = similar(a)
r = ccall((:fmpz_mod_mat_inv, libflint), Int,
(Ref{T}, Ref{T}), z, a)
!Bool(r) && error("Matrix not invertible")
return z
end

=#

################################################################################
#
#  Linear solving
#
################################################################################

function Solve._can_solve_internal_no_check(::Solve.LUTrait, A::ZZModMatrix, b::ZZModMatrix, task::Symbol; side::Symbol = :left)
  if side === :left
    fl, sol, K = Solve._can_solve_internal_no_check(Solve.LUTrait(), transpose(A), transpose(b), task, side = :right)
    return fl, transpose(sol), transpose(K)
  end

  x = similar(A, ncols(A), ncols(b))
  # This is probably only correct if the characteristic is prime
  fl = ccall((:fmpz_mod_mat_can_solve, libflint), Cint,
             (Ref{ZZModMatrix}, Ref{ZZModMatrix}, Ref{ZZModMatrix}, Ref{fmpz_mod_ctx_struct}), x, A, b, base_ring(x).ninv)
  if task === :only_check || task === :with_solution
    return Bool(fl), x, zero(A, 0, 0)
  end
  return Bool(fl), x, kernel(A, side = :right)
end

# For _can_solve_internal_no_check(::HowellFormTrait, ...) we use generic
# AbstractAlgebra code

################################################################################
#
#  LU decomposition
#
################################################################################

#= Not implemented in Flint yet

function lu!(P::Perm, x::T) where T <: Zmod_fmpz_mat
P.d .-= 1

rank = Int(ccall((:fmpz_mod_mat_lu, libflint), Cint, (Ptr{Int}, Ref{T}, Cint),
P.d, x, 0))

P.d .+= 1

# flint does x == PLU instead of Px == LU (docs are wrong)
inv!(P)

return rank
end

function lu(x::T, P = SymmetricGroup(nrows(x))) where T <: Zmod_fmpz_mat
m = nrows(x)
n = ncols(x)
P.n != m && error("Permutation does not match matrix")
p = one(P)
R = base_ring(x)
U = deepcopy(x)

L = similar(x, m, m)

rank = lu!(p, U)

for i = 1:m
for j = 1:n
if i > j
L[i, j] = U[i, j]
U[i, j] = R()
elseif i == j
L[i, j] = R(1)
elseif j <= m
L[i, j] = R()
end
end
end
return rank, p, L, U
end

=#

################################################################################
#
#  Windowing
#
################################################################################

function Base.view(x::ZZModMatrix, r1::Int, c1::Int, r2::Int, c2::Int)

  _checkrange_or_empty(nrows(x), r1, r2) ||
  Base.throw_boundserror(x, (r1:r2, c1:c2))

  _checkrange_or_empty(ncols(x), c1, c2) ||
  Base.throw_boundserror(x, (r1:r2, c1:c2))

  if (r1 > r2)
    r1 = 1
    r2 = 0
  end
  if (c1 > c2)
    c1 = 1
    c2 = 0
  end

  z = ZZModMatrix()
  z.base_ring = x.base_ring
  z.view_parent = x
  ccall((:fmpz_mod_mat_window_init, libflint), Nothing,
        (Ref{ZZModMatrix}, Ref{ZZModMatrix}, Int, Int, Int, Int, Ref{fmpz_mod_ctx_struct}),
        z, x, r1 - 1, c1 - 1, r2, c2, base_ring(x).ninv)
  finalizer(_fmpz_mod_mat_window_clear_fn, z)
  return z
end

function Base.view(x::T, r::AbstractUnitRange{Int}, c::AbstractUnitRange{Int}) where T <: Zmod_fmpz_mat
  return Base.view(x, first(r), first(c), last(r), last(c))
end

function _fmpz_mod_mat_window_clear_fn(a::ZZModMatrix)
  ccall((:fmpz_mod_mat_window_clear, libflint), Nothing,
        (Ref{ZZModMatrix}, Ref{fmpz_mod_ctx_struct}), a, base_ring(a).ninv)
end

function sub(x::T, r1::Int, c1::Int, r2::Int, c2::Int) where T <: Zmod_fmpz_mat
  return deepcopy(Base.view(x, r1, c1, r2, c2))
end

function sub(x::T, r::AbstractUnitRange{Int}, c::AbstractUnitRange{Int}) where T <: Zmod_fmpz_mat
  return deepcopy(Base.view(x, r, c))
end

function getindex(x::T, r::AbstractUnitRange{Int}, c::AbstractUnitRange{Int}) where T <: Zmod_fmpz_mat
  sub(x, r, c)
end

################################################################################
#
#  Concatenation
#
################################################################################

function hcat(x::T, y::T) where T <: Zmod_fmpz_mat
  (base_ring(x) != base_ring(y)) && error("Matrices must have same base ring")
  (x.r != y.r) && error("Matrices must have same number of rows")
  z = similar(x, nrows(x), ncols(x) + ncols(y))
  ccall((:fmpz_mod_mat_concat_horizontal, libflint), Nothing,
        (Ref{T}, Ref{T}, Ref{T}, Ref{fmpz_mod_ctx_struct}),
        z, x, y, base_ring(x).ninv)
  return z
end

function vcat(x::T, y::T) where T <: Zmod_fmpz_mat
  (base_ring(x) != base_ring(y)) && error("Matrices must have same base ring")
  (x.c != y.c) && error("Matrices must have same number of columns")
  z = similar(x, nrows(x) + nrows(y), ncols(x))
  ccall((:fmpz_mod_mat_concat_vertical, libflint), Nothing,
        (Ref{T}, Ref{T}, Ref{T}, Ref{fmpz_mod_ctx_struct}),
        z, x, y, base_ring(x).ninv)
  return z
end

################################################################################
#
#  Lifting
#
################################################################################

#= Not implemented in Flint yet

@doc raw"""
    lift(a::T) where {T <: Zmod_fmpz_mat}

Return a lift of the matrix $a$ to a matrix over $\mathbb{Z}$, i.e. where the
entries of the returned matrix are those of $a$ lifted to $\mathbb{Z}$.
"""
function lift(a::T) where {T <: Zmod_fmpz_mat}
z = ZZMatrix(nrows(a), ncols(a))
ccall((:fmpz_mat_set_fmpz_mod_mat, libflint), Nothing,
(Ref{ZZMatrix}, Ref{T}), z, a)
return z
end

function lift!(z::ZZMatrix, a::T) where T <: Zmod_fmpz_mat
ccall((:fmpz_mat_set_fmpz_mod_mat, libflint), Nothing,
(Ref{ZZMatrix}, Ref{T}), z, a)
return z
end

=#

@doc raw"""
    lift(a::T) where {T <: Zmod_fmpz_mat}

Return a lift of the matrix $a$ to a matrix over $\mathbb{Z}$, i.e. where the
entries of the returned matrix are those of $a$ lifted to $\mathbb{Z}$.
"""
function lift(a::Zmod_fmpz_mat)
  z = zero_matrix(ZZ, nrows(a), ncols(a))
  GC.@preserve a z begin
    for i in 1:nrows(a)
      for j in 1:ncols(a)
        m = mat_entry_ptr(z, i, j)
        n = mat_entry_ptr(a, i, j)
        set!(m, n)
      end
    end
  end
  return z
end

################################################################################
#
#  Characteristic polynomial
#
################################################################################

#= Not implemented in Flint yet

function charpoly(R::ZZModPolyRing, a::ZZModMatrix)
m = deepcopy(a)
p = R()
ccall((:fmpz_mod_mat_charpoly, libflint), Nothing,
(Ref{ZZModPolyRingElem}, Ref{ZZModMatrix}), p, m)
return p
end

=#

################################################################################
#
#  Minimal polynomial
#
################################################################################

#= Not implemented in Flint yet

function minpoly(R::ZZModPolyRing, a::ZZModMatrix)
p = R()
ccall((:fmpz_mod_mat_minpoly, libflint), Nothing,
(Ref{ZZModPolyRingElem}, Ref{ZZModMatrix}), p, a)
return p
end

=#

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{ZZModMatrix}, ::Type{V}) where {V <: Integer} = ZZModMatrix

promote_rule(::Type{ZZModMatrix}, ::Type{ZZModRingElem}) = ZZModMatrix

promote_rule(::Type{ZZModMatrix}, ::Type{ZZRingElem}) = ZZModMatrix

################################################################################
#
#  Parent object overloading
#
################################################################################

function (a::ZZModMatrixSpace)()
  z = ZZModMatrix(base_ring(a), nrows(a), ncols(a))
  return z
end

function (a::ZZModMatrixSpace)(arr::AbstractMatrix{BigInt}, transpose::Bool = false)
  _check_dim(nrows(a), ncols(a), arr, transpose)
  z = ZZModMatrix(nrows(a), ncols(a), base_ring(a).ninv, arr, transpose)
  z.base_ring = a.base_ring
  return z
end

function (a::ZZModMatrixSpace)(arr::AbstractVector{BigInt})
  _check_dim(nrows(a), ncols(a), arr)
  z = ZZModMatrix(nrows(a), ncols(a), base_ring(a).ninv, arr)
  z.base_ring = a.base_ring
  return z
end

function (a::ZZModMatrixSpace)(arr::AbstractMatrix{ZZRingElem}, transpose::Bool = false)
  _check_dim(nrows(a), ncols(a), arr, transpose)
  z = ZZModMatrix(nrows(a), ncols(a), base_ring(a).ninv, arr, transpose)
  z.base_ring = a.base_ring
  return z
end

function (a::ZZModMatrixSpace)(arr::AbstractVector{ZZRingElem})
  _check_dim(nrows(a), ncols(a), arr)
  z = ZZModMatrix(nrows(a), ncols(a), base_ring(a).ninv, arr)
  z.base_ring = a.base_ring
  return z
end

function (a::ZZModMatrixSpace)(arr::AbstractMatrix{Int}, transpose::Bool = false)
  _check_dim(nrows(a), ncols(a), arr, transpose)
  z = ZZModMatrix(nrows(a), ncols(a), base_ring(a).ninv, arr, transpose)
  z.base_ring = a.base_ring
  return z
end

function (a::ZZModMatrixSpace)(arr::AbstractVector{Int})
  _check_dim(nrows(a), ncols(a), arr)
  z = ZZModMatrix(nrows(a), ncols(a), base_ring(a).ninv, arr)
  z.base_ring = a.base_ring
  return z
end

function (a::ZZModMatrixSpace)(arr::AbstractMatrix{ZZModRingElem}, transpose::Bool = false)
  _check_dim(nrows(a), ncols(a), arr, transpose)
  (length(arr) > 0 && (base_ring(a) != parent(arr[1]))) && error("Elements must have same base ring")
  z = ZZModMatrix(nrows(a), ncols(a), base_ring(a).ninv, arr, transpose)
  z.base_ring = a.base_ring
  return z
end

function (a::ZZModMatrixSpace)(arr::AbstractVector{ZZModRingElem})
  _check_dim(nrows(a), ncols(a), arr)
  (length(arr) > 0 && (base_ring(a) != parent(arr[1]))) && error("Elements must have same base ring")
  z = ZZModMatrix(nrows(a), ncols(a), base_ring(a).ninv, arr)
  z.base_ring = a.base_ring
  return z
end

function (a::ZZModMatrixSpace)(b::ZZMatrix)
  (ncols(a) != b.c || nrows(a) != b.r) && error("Dimensions do not fit")
  z = ZZModMatrix(base_ring(a), b.r, b.c)
  ccall((:fmpz_mod_mat_set_fmpz_mat, libflint), Nothing,
        (Ref{ZZModMatrix}, Ref{ZZMatrix}, Ref{fmpz_mod_ctx_struct}),
        z, b, base_ring(a).ninv)
  return z
end

###############################################################################
#
#   Matrix constructor
#
###############################################################################

function matrix(R::ZZModRing, arr::AbstractMatrix{<: Union{ZZModRingElem, ZZRingElem, Integer}})
  z = ZZModMatrix(size(arr, 1), size(arr, 2), R.ninv, arr)
  z.base_ring = R
  return z
end

function matrix(R::ZZModRing, r::Int, c::Int, arr::AbstractVector{<: Union{ZZModRingElem, ZZRingElem, Integer}})
  _check_dim(r, c, arr)
  z = ZZModMatrix(r, c, R.ninv, arr)
  z.base_ring = R
  return z
end

###############################################################################
#
#  Zero matrix
#
###############################################################################

function zero_matrix(R::ZZModRing, r::Int, c::Int)
  if r < 0 || c < 0
    error("dimensions must not be negative")
  end
  z = ZZModMatrix(R, r, c)
  return z
end

################################################################################
#
#  Kernel
#
################################################################################

# kernel(::ZZModMatrix) is in src/flint/nmod_mat.jl

function nullspace(M::ZZModMatrix)
  # Apparently this only works correctly if base_ring(M) is a field
  N = similar(M, ncols(M), ncols(M))
  nullity = ccall((:fmpz_mod_mat_nullspace, libflint), Int,
                  (Ref{ZZModMatrix}, Ref{ZZModMatrix}, Ref{fmpz_mod_ctx_struct}), N, M, base_ring(N).ninv)
  return nullity, view(N, 1:nrows(N), 1:nullity)
end

################################################################################
#
#  Entry pointers
#
################################################################################

mat_entry_ptr(A::ZZModMatrix, i::Int, j::Int) = unsafe_load(A.rows, i) + (j-1)*sizeof(ZZRingElem)
