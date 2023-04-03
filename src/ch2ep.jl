# Compute Effective ELASTIC Properties
using JSON
using SparseArrays
using LinearAlgebra

# Model data struct:
struct Model
    nx::UInt64;
    ny::UInt64;
    voxelSize::Float64;
    refinement::UInt64;
    nMat::UInt16;
    rhsType::UInt8;
    solverType::UInt8;
    pcgTol::Float64;
    pcgIter::UInt64;
    matKeys::Vector{UInt16};
    matProp::Matrix{Float64}; 
    nNodes::UInt64;
    nElems::UInt64;
    nDOFs::UInt64;
    DOFMap::Vector{UInt64}; 
    elemMatMap::Vector{UInt16}; 
    function Model(_nx, _ny, _voxelSize, _refinement, _nMat, _rhsType, _solverType, _pcgTol, _pcgIter, _matKeys, _matProp, _nNodes, _nElems, _nDOFs, _DOFMap, _elemMatMap)
        new(_nx, _ny, _voxelSize, _refinement, _nMat, _rhsType, _solverType, _pcgTol, _pcgIter, _matKeys, _matProp, _nNodes, _nElems, _nDOFs, _DOFMap, _elemMatMap);
    end
end

# Build Model:
function buildModel(_JsonFile::String, _RawFile::String)
    println(".Building Model!")
    # Read Json file:
    nx, ny, voxelSize, refinement,  nMat, rhsType, solverType, pcgTol, pcgIter, matKeys, matProp = readJSON(_JsonFile)
    # Read Raw file:
    elemMatMap = zeros(UInt16, nx * ny * refinement * refinement)
    readRAW!(nx, ny, refinement, matKeys, elemMatMap, _RawFile)
    # Update the parameters based on the given refinement level:
    nx *= refinement
    ny *= refinement
    nNodes::UInt64 = (nx + 1) * (ny + 1)
    nElems::UInt64 = (nx) * (ny)
    DOFperNode::UInt64 = 2
    nDOFs::UInt64 = nElems * DOFperNode
    # Generate a map of Degree of Freedom:
    DOFMap = zeros(UInt64, nNodes)
    generateDOFMap!(nx, ny, DOFMap)
    # Build the Model:
    model = Model(nx, ny, voxelSize, refinement, nMat, rhsType, solverType, pcgTol, pcgIter, matKeys, matProp, nNodes, nElems, nDOFs, DOFMap, elemMatMap)
    println("---------------------------")
    return model
end

# Read JSON file:
function readJSON(_filename::String)
    println("   .Reading JSON!")
    # Open and read file:
    open(_filename, "r") do f
        data = JSON.parse(f)
        nx::UInt64 = data["image_dimensions"][1]
        ny::UInt64 = data["image_dimensions"][2]
        refinement::UInt64 = 1
        if haskey(data, "refinement"); refinement = data["refinement"]; end
        voxelSize::Float64 = 1.0
        if haskey(data, "voxel_size"); voxelSize = data["voxel_size"]; end
        rhsType::UInt8 = 0
        if haskey(data, "type_of_rhs"); rhsType = data["type_of_rhs"]; end
        solverType::UInt8 = 0
        if haskey(data, "type_of_solver"); solverType = data["type_of_solver"]; end
        pcgTol::Float64 = 0.000001
        if haskey(data, "solver_tolerance"); pcgTol = data["solver_tolerance"]; end
        pcgIter::UInt64 = nx * ny * refinement * refinement
        if haskey(data, "number_of_iterations"); pcgIter = data["number_of_iterations"]; end
        nMat::UInt16 = data["number_of_materials"]
        materials = data["properties_of_materials"]
        matKeys = zeros(UInt16, 256)
        matProp = zeros(Float64, 256, 2)
        for i = 1:nMat
            matKeys[convert(UInt8, materials[i][1]) + 1] = i
            matProp[convert(UInt8, materials[i][1]) + 1,1] = convert(Float64, materials[i][2])
            matProp[convert(UInt8, materials[i][1]) + 1,2] = convert(Float64, materials[i][3])
        end
        materials = nothing
        data = nothing
        return nx, ny, voxelSize, refinement, nMat, rhsType, solverType, pcgTol, pcgIter, matKeys, matProp
    end
end

# Read RAW file:
function readRAW!(_nx::UInt64, _ny::UInt64, _refinement::UInt64, _matKeys::Vector{UInt16}, _elemMatMap::Vector{UInt16}, _filename::String)
    println("   .Reading RAW!")
    # Initializations:
    nelem::UInt64 = _nx * _ny;  buffer::UInt64 = 0;    
    elref::UInt64 = 0; ix::UInt64 = 0; iy::UInt64 = 0;
    row::UInt64 = _ny * _refinement; rowref::UInt64 = _ny * _refinement * _refinement; 
    # Open and read file:
    open(_filename, "r") do io
        bin_array = read(io)      
        # Build the element material map based on the refinement level:
        for e = 1:nelem
            buffer = _matKeys[bin_array[e] + 1]
            ix = ((e-1) % _nx)
            iy = ((e-1) ÷ _nx)    
            # el = 1 + ix * _ny + iy 
            elref = 1 + (ix * rowref) + (iy*_refinement)  
            for i = 1:_refinement
                for j = 1:_refinement     
                    _elemMatMap[elref + (j-1) + (i-1)*row] = buffer
                end        
            end 
        end     
        bin_array = nothing
    end
end

# Generate the Degree of Freedom Map:
function generateDOFMap!(_nx::UInt64, _ny::UInt64, _DOFMap::Vector{UInt64})
    println("   .Generating the Map of DOFs (Degrees of Freedom)!")
    # Number the DOFs following the nodes from top to bottom and left to right:
    nElems::UInt64 = _nx * _ny;
    nNodes::UInt64 = (_nx + 1) * (_ny + 1);
    @fastmath @inbounds @simd for n = 1:nNodes
        i = (n - 1) % nNodes
        _DOFMap[n] = (n - 1 - div(n - 1, (_ny + 1)) - _ny * (div((n - 1) % (_ny + 1), _ny))) % nElems + 1 
    end
end

# Estimate memory consuption:
function estimateMemory(_model::Model)
    println(".Estimating memory!")
    # elemMatMap = 16 bits * nElems
    # DOFMap = 64 bits * nNodes
    # RHS = 64 bits * nDOFs
    # PCGM Solver   / solverType == 0 / M d x q =  * 64 bits * nDOFs
    # Direct Solver / solverType == 1 / K = 18 * 64 bits * nElems (rough sparse estimative)
    mem::Float64 = 0.0
    if (_model.solverType == 0)
        mem += (16 * _model.nElems + 64 * _model.nNodes + 5 * 64 * _model.nDOFs) / 8 / 1_000_000
    elseif (_model.solverType == 1)
        mem += (16 * _model.nElems + 64 * _model.nNodes + 2 * 64 * _model.nDOFs + 18 * 64 * _model.nElems) / 8 / 1_000_000
    end
    println("   $(_model.nDOFs) DOFs")
    println("   $mem MB")
    println("---------------------------")
end

# Compute the element stiffness matrix for each material:
function elementStiffnessMatrices!(_model::Model, _K::Array{Float64,3}, _B::Array{Float64,3})
    println("   .Computing each element stiffness matrix!")
    # Compute the matrices for each material:
    i::UInt64 = 0
    for j = 1:256
        if (_model.matKeys[j] != 0)
            i += 1
            elemProps = _model.matProp[j,:]
            _K[:,:,i], _B[:,:,i] = Q4ElementStiffness(elemProps)
        end
    end
end

# Element Q4 Stiffness - FEM:
function Q4ElementStiffness(_elemProps::Vector{Float64})
    # Initializations
    k  = zeros(Float64, 8, 8)
    BC = zeros(Float64, 3, 8)
    C  = zeros(Float64, 3, 3)
    E::Float64 = _elemProps[1]
    p::Float64 = _elemProps[2]
    # Element coords 
    #  4*___*3
    #   |   |
    #   |   |
    #  1*___*2
    x = [0.;1.;1.;0.]
    y = [0.;0.;1.;1.]
    # Constitutive matrix (Plane Strain State)
    E /= ((1. + p) * (1. - 2 * p))
    C[1,1] = 1 - p;  C[1,2] = p;
    C[2,1] = p;      C[2,2] = 1 - p;
    C[3,3] = (1 - 2 * p) / 2;
    C .*= E
    # Constitutive matrix (Plane Stress State)
    # E /= (1. - p*p);
    # C[1,1] = 1;    C[1,2] = p;
    # C[2,1] = p;    C[2,2] = 1;
    # C[3,3] = (1-p)/2;
    # C .*= E;
    # Gauss Points and Weights
    gp = [-1.0 / sqrt(3) 1.0 / sqrt(3)]
    # w = [1.0 1.0];
    for i = 1:2
        r = gp[1,i]
        for j = 1:2
            s = gp[1,j]
            B, detJ = Q4BMatrix(r, s, x, y)
            k  += B' * C * B * detJ# *w[1,i]*w[1,j];
            BC += C * B * detJ# *w[1,i]*w[1,j];
        end
    end
    return k, BC
end

# Q4BMatrix - FEM
function Q4BMatrix(r::Float64, s::Float64, x::Vector{Float64}, y::Vector{Float64})
    # Initializations
    X = [x'; y']
    # Compute B matrix and Jacobian
    dN1dr = -(1 - s) * .25; dN2dr =  (1 - s) * .25; dN3dr = (1 + s) * .25; dN4dr = -(1 + s) * .25;
    dN1ds = -(1 - r) * .25; dN2ds = -(1 + r) * .25; dN3ds = (1 + r) * .25; dN4ds =  (1 - r) * .25;
    dN = [ dN1dr dN2dr dN3dr dN4dr;
           dN1ds dN2ds dN3ds dN4ds ];
    J = dN * X'
    dNxy = J \ dN
    B = [ dNxy[1,1]        0  dNxy[1,2]        0  dNxy[1,3]        0  dNxy[1,4]        0  ;
                 0  dNxy[2,1]        0  dNxy[2,2]        0  dNxy[2,3]        0  dNxy[2,4] ;
          dNxy[2,1] dNxy[1,1] dNxy[2,2] dNxy[1,2] dNxy[2,3] dNxy[1,3] dNxy[2,4] dNxy[1,4] ];
    return B, det(J)
end

# Compute the RHS: Boundary or Domain, rhsType: Boundary = 0 || Domain = 1, axis 0 = X || axis 1 = Y
function computeRHS!(_model::Model, _RHS::Vector{Float64}, _axis::Int, _K::Array{Float64,3}, _B::Array{Float64,3})    
    # Initializations:
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0;  e::UInt64 = 0;
    pElemDOFNum = zeros(UInt64, 8)
    # Compute each RHS (_axis) based on boundary or domain data (_model.m_rhsType):
    if _model.rhsType == 1  # Boundary 
        println("   .Computing RHS - Boundary!")
        delta::Float64 = 0.0        
        if _axis == 0     # _axis 0 = X
            delta = _model.nx
            for e = (_model.nElems - _model.ny + 1):_model.nElems
                N1 = e + 1 + div(e - 1, _model.ny); N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;                
                pElemDOFNum[1] = _model.DOFMap[N1] * 2 - 1; pElemDOFNum[2] = _model.DOFMap[N1] * 2;
                pElemDOFNum[3] = _model.DOFMap[N2] * 2 - 1; pElemDOFNum[4] = _model.DOFMap[N2] * 2;
                pElemDOFNum[5] = _model.DOFMap[N3] * 2 - 1; pElemDOFNum[6] = _model.DOFMap[N3] * 2;
                pElemDOFNum[7] = _model.DOFMap[N4] * 2 - 1; pElemDOFNum[8] = _model.DOFMap[N4] * 2;
                for i = 1:8
                    for j in [3 5]
                        _RHS[pElemDOFNum[i]] -= (_K[i,j,_model.elemMatMap[e]]) * delta
                    end
                end
            end 
        elseif _axis == 1 # _axis 1 = Y
            delta = _model.ny
            for e = 1:_model.ny:(_model.nElems - _model.ny + 1)
                N1 = e + 1 + div(e - 1, _model.ny); N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;                
                pElemDOFNum[1] = _model.DOFMap[N1] * 2 - 1; pElemDOFNum[2] = _model.DOFMap[N1] * 2;
                pElemDOFNum[3] = _model.DOFMap[N2] * 2 - 1; pElemDOFNum[4] = _model.DOFMap[N2] * 2;
                pElemDOFNum[5] = _model.DOFMap[N3] * 2 - 1; pElemDOFNum[6] = _model.DOFMap[N3] * 2;
                pElemDOFNum[7] = _model.DOFMap[N4] * 2 - 1; pElemDOFNum[8] = _model.DOFMap[N4] * 2;
                for i = 1:8
                    for j in [6 8]
                        _RHS[pElemDOFNum[i]] -= (_K[i,j,_model.elemMatMap[e]]) * delta
                    end
                end
            end
        elseif _axis == 2 # _axis 2 = XY
            delta = _model.ny
            for e = 1:_model.ny:(_model.nElems - _model.ny + 1)
                N1 = e + 1 + div(e - 1, _model.ny); N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                pElemDOFNum[1] = _model.DOFMap[N1] * 2 - 1; pElemDOFNum[2] = _model.DOFMap[N1] * 2;
                pElemDOFNum[3] = _model.DOFMap[N2] * 2 - 1; pElemDOFNum[4] = _model.DOFMap[N2] * 2;
                pElemDOFNum[5] = _model.DOFMap[N3] * 2 - 1; pElemDOFNum[6] = _model.DOFMap[N3] * 2;
                pElemDOFNum[7] = _model.DOFMap[N4] * 2 - 1; pElemDOFNum[8] = _model.DOFMap[N4] * 2;
                for i = 1:8
                    for j in [5 7]
                        _RHS[pElemDOFNum[i]] -= (_K[i,j,_model.elemMatMap[e]]) * delta
                    end
                end
            end
        end
    elseif _model.rhsType == 0  # Domain  
        println("   .Computing RHS - Domain!")   
        for e = 1:_model.nElems
            N1 = e + ((e - 1) ÷ _model.ny) + 1; N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;            
            pElemDOFNum[1] = _model.DOFMap[N1] * 2 - 1; pElemDOFNum[2] = _model.DOFMap[N1] * 2;
            pElemDOFNum[3] = _model.DOFMap[N2] * 2 - 1; pElemDOFNum[4] = _model.DOFMap[N2] * 2;
            pElemDOFNum[5] = _model.DOFMap[N3] * 2 - 1; pElemDOFNum[6] = _model.DOFMap[N3] * 2;
            pElemDOFNum[7] = _model.DOFMap[N4] * 2 - 1; pElemDOFNum[8] = _model.DOFMap[N4] * 2;
            for i = 1:8
                _RHS[pElemDOFNum[i]] += _B[_axis + 1,i,_model.elemMatMap[e]]
            end
        end
    end
end

# Direct Solver: [K] 64 bits * m_nDOFs * m_nDOFs 
function directMethod!(_model::Model, _x1::Vector{Float64}, _x2::Vector{Float64}, _x3::Vector{Float64}, _RHS1::Vector{Float64}, _RHS2::Vector{Float64}, _RHS3::Vector{Float64}, _K::Array{Float64,3})
    println("   .Direct Solver!")
    # Initializations:
    K = spzeros(_model.nDOFs, _model.nDOFs)
    pElemDOFNum = zeros(UInt64, 8)
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0; 
    thisElemMat::UInt16 = 0
    # Assembly system matrix:
    for e = 1:_model.nElems
        thisElemMat = _model.elemMatMap[e]
        N1 = e + 1 + div(e - 1, _model.ny); N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;        
        pElemDOFNum[1] = _model.DOFMap[N1] * 2 - 1; pElemDOFNum[2] = _model.DOFMap[N1] * 2;
        pElemDOFNum[3] = _model.DOFMap[N2] * 2 - 1; pElemDOFNum[4] = _model.DOFMap[N2] * 2;
        pElemDOFNum[5] = _model.DOFMap[N3] * 2 - 1; pElemDOFNum[6] = _model.DOFMap[N3] * 2;
        pElemDOFNum[7] = _model.DOFMap[N4] * 2 - 1; pElemDOFNum[8] = _model.DOFMap[N4] * 2;
        for i = 1:8
            for j = 1:8
                K[pElemDOFNum[i],pElemDOFNum[j]] += _K[i,j,thisElemMat]
            end
        end
    end
    # Solve for three rhs:
    _x1 .= K \ _RHS1
    _x2 .= K \ _RHS2
    _x3 .= K \ _RHS3
end

# Jacobi Preconditioner: assembly || M
function jacobiPrecond!(_model::Model, _M::Vector{Float64}, _K::Array{Float64,3})
    println("   .Jacobi Preconditioner!")
    # Initializations:
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0; 
    thisElemMat::UInt16 = 0
    pElemDOFNum = zeros(UInt64, 8)
    # Compute the preconditioner: 
    for e = 1:_model.nElems
        thisElemMat = _model.elemMatMap[e]
        N1 = e + 1 + div(e - 1, _model.ny); N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;        
        pElemDOFNum[1] = _model.DOFMap[N1] * 2 - 1; pElemDOFNum[2] = _model.DOFMap[N1] * 2;
        pElemDOFNum[3] = _model.DOFMap[N2] * 2 - 1; pElemDOFNum[4] = _model.DOFMap[N2] * 2;
        pElemDOFNum[5] = _model.DOFMap[N3] * 2 - 1; pElemDOFNum[6] = _model.DOFMap[N3] * 2;
        pElemDOFNum[7] = _model.DOFMap[N4] * 2 - 1; pElemDOFNum[8] = _model.DOFMap[N4] * 2;
        for i = 1:8
            _M[pElemDOFNum[i]] += _K[i,i,thisElemMat]
        end
    end
    _M .= _M .\ 1
end 

# Preconditioned Conjugate Gradient Method:
function pcg_old!(_model::Model, _x::Vector{Float64}, _r::Vector{Float64}, _M::Vector{Float64}, _K::Array{Float64,3})
    println("   .PCG Solver!");
    # Initializations:    
    d = zeros(Float64, _model.nDOFs);
    q = zeros(Float64, _model.nDOFs);
    pElemDOFNum = zeros(UInt64, 8);
    N1 = 0; N2 = 0; N3 = 0; N4 = 0; 
    q_temp = 0;    
    # PCG Initialization:
    d .= _r;
    d .*= _M;
    delta_new = (_r' * d)[1,1];
    delta_0 = delta_new;
    i_max = _model.pcgIter;
    ii = 0;
    # PCG Iterations:
    while (ii < i_max) && (abs(delta_new) > _model.pcgTol * _model.pcgTol * abs(delta_0)) # (maximum(abs.(_r))>_pcgTol)
        @fastmath @inbounds @simd for e = 1:_model.nElems
            N1 = e + ((e - 1) ÷ _model.ny) + 1; N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            pElemDOFNum[1] = _model.DOFMap[N1] * 2 - 1; pElemDOFNum[2] = _model.DOFMap[N1] * 2;
            pElemDOFNum[3] = _model.DOFMap[N2] * 2 - 1; pElemDOFNum[4] = _model.DOFMap[N2] * 2;
            pElemDOFNum[5] = _model.DOFMap[N3] * 2 - 1; pElemDOFNum[6] = _model.DOFMap[N3] * 2;
            pElemDOFNum[7] = _model.DOFMap[N4] * 2 - 1; pElemDOFNum[8] = _model.DOFMap[N4] * 2;
            for i = 1:8
                q_temp = 0;
                for j = 1:8
                    q_temp += _K[i,j,_model.elemMatMap[e]] * d[pElemDOFNum[j]];
                end
                q[pElemDOFNum[i]] += q_temp;
            end
        end
        alfa = delta_new / (d' * q)[1,1];
        d .*= alfa;
        _x .+= d;
        q .*= alfa;
        _r .-= q;
        q .= _r;
        q .*= _M;
        delta_old = delta_new;
        delta_new = (_r' * q)[1,1];
        beta = delta_new / delta_old;
        d .*= beta / alfa;
        d .+= q;
        q .*= 0;
        ii += 1;
    end
    println("    $ii steps");
    println("    Residue = ", sqrt(abs(delta_new) / abs(delta_0)));
end

# Preconditioned Conjugate Gradient Method:
function pcg!(_model::Model, _x::Vector{Float64}, _r::Vector{Float64}, _M::Vector{Float64}, _K::Array{Float64,3})
    println("   .PCG Solver!")
    # Initializations:    
    d = zeros(Float64, _model.nDOFs)
    q = zeros(Float64, _model.nDOFs)
    pElemDOFNum = zeros(UInt64, 8)
    pElemDOFVar = zeros(Float64, 8)
    thisElemMat::UInt16 = 0
    q_temp::Float64 = 0.0
    alfa::Float64 = 0.0
    beta::Float64 = 0.0
    ri::Float64   = 0.0
    qi::Float64   = 0.0
    delta_new::Float64 = 0.0
    delta_old::Float64 = 0.0
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0;   
    # PCG Initialization:
    @inbounds for i=1:_model.nDOFs; d[i] = _r[i]*_M[i]; end
    delta_new = dot(_r,d)
    delta_0 = delta_new
    if (abs(delta_0)<1e-14); println("    x0 satisfied absolute tolerance criteria: delta < 1e-14"); return; end
    tolerance::Float64 = _model.pcgTol * _model.pcgTol * abs(delta_0)
    iteration_count::UInt64 = _model.pcgIter
    # PCG Iterations:
    for ii = 1:_model.pcgIter
        # (EbE) q = Kd
        @fastmath @inbounds @simd for e = 1:_model.nElems
            thisElemMat = _model.elemMatMap[e]
            N1 = e + ((e - 1) ÷ _model.ny) + 1; N2 = N1 + _model.ny + 1; N3 = N2 - 1; N4 = N1 - 1;            
            pElemDOFNum[1] = 2*_model.DOFMap[N1]-1; pElemDOFNum[2] = pElemDOFNum[1]+1; 
            pElemDOFNum[3] = 2*_model.DOFMap[N2]-1; pElemDOFNum[4] = pElemDOFNum[3]+1;
            pElemDOFNum[5] = 2*_model.DOFMap[N3]-1; pElemDOFNum[6] = pElemDOFNum[5]+1;
            pElemDOFNum[7] = 2*_model.DOFMap[N4]-1; pElemDOFNum[8] = pElemDOFNum[7]+1;
            pElemDOFVar[1] = d[pElemDOFNum[1]]; pElemDOFVar[2] = d[pElemDOFNum[2]]; pElemDOFVar[3] = d[pElemDOFNum[3]]; pElemDOFVar[4] = d[pElemDOFNum[4]];
            pElemDOFVar[5] = d[pElemDOFNum[5]]; pElemDOFVar[6] = d[pElemDOFNum[6]]; pElemDOFVar[7] = d[pElemDOFNum[7]]; pElemDOFVar[8] = d[pElemDOFNum[8]];
            @inbounds for i = 1:8
                q_temp = _K[i,1,thisElemMat] * pElemDOFVar[1]
                q_temp+= _K[i,2,thisElemMat] * pElemDOFVar[2]
                q_temp+= _K[i,3,thisElemMat] * pElemDOFVar[3]
                q_temp+= _K[i,4,thisElemMat] * pElemDOFVar[4]
                q_temp+= _K[i,5,thisElemMat] * pElemDOFVar[5]
                q_temp+= _K[i,6,thisElemMat] * pElemDOFVar[6]
                q_temp+= _K[i,7,thisElemMat] * pElemDOFVar[7]
                q_temp+= _K[i,8,thisElemMat] * pElemDOFVar[8]
                q[pElemDOFNum[i]] += q_temp;
            end
        end
        alfa = delta_new / dot(d,q)
        delta_old = delta_new
        delta_new = 0.0
        @inbounds for i=1:_model.nDOFs
            _x[i]+=d[i]*alfa
            _r[i]-=q[i]*alfa
            ri   =_r[i]
            qi   =  ri*_M[i]
            q[i]  =qi
            delta_new+=ri*qi       
        end
        if (abs(delta_new) <= tolerance); iteration_count = ii; break; end
        beta = delta_new / delta_old
        @inbounds for i=1:_model.nDOFs
            d[i]*=beta
            d[i]+=q[i]
            q[i] = 0.0
        end
    end
    println("    $iteration_count steps")
    println("    Residue = ", sqrt(abs(delta_new) / abs(delta_0)))
end

# Compute Stress-FEM Effective property:
function femEffective(_model::Model, _X::Vector{Float64}, _axis::Int, _B::Array{Float64,3})
    println("   .Updating Constitutive Matrix!")
    # Initializations:
    SX::Float64 = 0.0; SY::Float64 = 0.0; SXY::Float64 = 0.0;
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0;
    pElemDOFNum = zeros(UInt64, 8)
    C = zeros(Float64, 3, 3)

    # Compute the effective properties for each test: 
    if _model.rhsType == 1  # Boundary
        if _axis == 0
            stressSX = open("out/data/stressSXX.txt", "w")
            stressSY = open("out/data/stressSYY.txt", "w")
            stressSXY = open("out/data/stressSXY.txt", "w")
            deltaX = _model.nx;   
            for eb = _model.nElems - (_model.ny - 1):_model.nElems
                for j in [3 5]
                    SX   += (_B[1,j,_model.elemMatMap[eb]] * deltaX)
                    SY   += (_B[2,j,_model.elemMatMap[eb]] * deltaX)
                    SXY  += (_B[3,j,_model.elemMatMap[eb]] * deltaX)
                end
                write(stressSX, string(SX, '\n'))
                write(stressSY, string(SY, '\n'))
                write(stressSXY, string(SXY, '\n'))
            end
            close(stressSX)
            close(stressSY)
            close(stressSXY)  
        elseif _axis == 1
            deltaX = _model.ny;
            for eb = 1:(_model.ny):_model.nElems
                for j in [6 8]
                    SX   += (_B[1,j,_model.elemMatMap[eb]] * deltaX)
                    SY   += (_B[2,j,_model.elemMatMap[eb]] * deltaX)
                    SXY  += (_B[3,j,_model.elemMatMap[eb]] * deltaX)
                end           
            end 
        elseif _axis == 2
            deltaX = _model.ny;
            for eb = 1:(_model.ny):_model.nElems
                for j in [5 7]
                    SX   += (_B[1,j,_model.elemMatMap[eb]] * deltaX)
                    SY   += (_B[2,j,_model.elemMatMap[eb]] * deltaX)
                    SXY  += (_B[3,j,_model.elemMatMap[eb]] * deltaX)
                end           
            end 
        end
        for e = 1:_model.nElems
            N1 = e + ((e - 1) ÷ _model.ny) + 1; N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            pElemDOFNum[1] = _model.DOFMap[N1] * 2 - 1; pElemDOFNum[2] = _model.DOFMap[N1] * 2;
            pElemDOFNum[3] = _model.DOFMap[N2] * 2 - 1; pElemDOFNum[4] = _model.DOFMap[N2] * 2;
            pElemDOFNum[5] = _model.DOFMap[N3] * 2 - 1; pElemDOFNum[6] = _model.DOFMap[N3] * 2;
            pElemDOFNum[7] = _model.DOFMap[N4] * 2 - 1; pElemDOFNum[8] = _model.DOFMap[N4] * 2;
            for i = 1:8
                SX  += (_B[1,i,_model.elemMatMap[e]] * _X[pElemDOFNum[i]])
                SY  += (_B[2,i,_model.elemMatMap[e]] * _X[pElemDOFNum[i]])
                SXY += (_B[3,i,_model.elemMatMap[e]] * _X[pElemDOFNum[i]])
            end
        end        
    elseif _model.rhsType == 0  # Domain
        x = zeros(Float64, 8)
        if (_axis == 0) # horizontal deformation
            x[3] = 1; x[5] = 1;
            stressSX = open("out/data/stressSXX.txt", "w")
            stressSY = open("out/data/stressSYY.txt", "w")
            stressSXY = open("out/data/stressSXY.txt", "w")
            for e = 1:_model.nElems            
                N1 = e + ((e - 1) ÷ _model.ny) + 1; N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                pElemDOFNum[1] = _model.DOFMap[N1] * 2 - 1; pElemDOFNum[2] = _model.DOFMap[N1] * 2;
                pElemDOFNum[3] = _model.DOFMap[N2] * 2 - 1; pElemDOFNum[4] = _model.DOFMap[N2] * 2;
                pElemDOFNum[5] = _model.DOFMap[N3] * 2 - 1; pElemDOFNum[6] = _model.DOFMap[N3] * 2;
                pElemDOFNum[7] = _model.DOFMap[N4] * 2 - 1; pElemDOFNum[8] = _model.DOFMap[N4] * 2;
                for i = 1:8
                    SX  += (_B[1,i,_model.elemMatMap[e]] * (x[i] - _X[pElemDOFNum[i]]))
                    SY  += (_B[2,i,_model.elemMatMap[e]] * (x[i] - _X[pElemDOFNum[i]]))
                    SXY += (_B[3,i,_model.elemMatMap[e]] * (x[i] - _X[pElemDOFNum[i]]))
                end
                write(stressSX, string((_B[1,1,_model.elemMatMap[e]] * (x[1] - _X[pElemDOFNum[1]])), '\n'))
                write(stressSY, string((_B[2,1,_model.elemMatMap[e]] * (x[1] - _X[pElemDOFNum[1]])), '\n'))
                write(stressSXY, string((_B[3,1,_model.elemMatMap[e]] * (x[1] - _X[pElemDOFNum[1]])), '\n'))
                #write(stressSX, string(SX, '\n'))
                #write(stressSY, string(SY, '\n'))
                #write(stressSXY, string(SXY, '\n'))
            end
            close(stressSX)
            close(stressSY)
            close(stressSXY)  
        else
            if (_axis == 1); x[6] = 1; x[8] = 1;
            elseif (_axis == 2); x[5] = 1; x[7] = 1; end

            for e = 1:_model.nElems            
                N1 = e + ((e - 1) ÷ _model.ny) + 1; N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                pElemDOFNum[1] = _model.DOFMap[N1] * 2 - 1; pElemDOFNum[2] = _model.DOFMap[N1] * 2;
                pElemDOFNum[3] = _model.DOFMap[N2] * 2 - 1; pElemDOFNum[4] = _model.DOFMap[N2] * 2;
                pElemDOFNum[5] = _model.DOFMap[N3] * 2 - 1; pElemDOFNum[6] = _model.DOFMap[N3] * 2;
                pElemDOFNum[7] = _model.DOFMap[N4] * 2 - 1; pElemDOFNum[8] = _model.DOFMap[N4] * 2;
                for i = 1:8
                    SX  += (_B[1,i,_model.elemMatMap[e]] * (x[i] - _X[pElemDOFNum[i]]))
                    SY  += (_B[2,i,_model.elemMatMap[e]] * (x[i] - _X[pElemDOFNum[i]]))
                    SXY += (_B[3,i,_model.elemMatMap[e]] * (x[i] - _X[pElemDOFNum[i]]))
                end
            end
        end
    end
    C[1,_axis + 1] = SX / _model.nElems; C[2,_axis + 1] = SY / _model.nElems; C[3,_axis + 1] = SXY / _model.nElems;
    return C
end

# -----------------
function homogenize(_arg)
    println("---------------------------")
    # Build the Model data struct:
    m_model = buildModel(_arg * ".json", _arg * ".raw")
    # Estimate Memory Consumption:
    estimateMemory(m_model)
    # SOLVE:
    println(".Solving")
    # Compute the stiffness matrix for each Material:
    m_K = zeros(Float64, 8, 8, m_model.nMat)
    m_B = zeros(Float64, 3, 8, m_model.nMat)
    elementStiffnessMatrices!(m_model, m_K, m_B)
    if (m_model.solverType == 0) # Preconditioned Conjugate Gradient Method
        # Initialize the effective tensor, the right hand side, the inicial guess and the preconditioner:
        m_C = zeros(Float64, 3, 3)
        m_RHS = zeros(Float64, m_model.nDOFs)  
        m_X = zeros(Float64, m_model.nDOFs)
        m_M = zeros(Float64, m_model.nDOFs)
        # Compute the Jacobi preconditioner:
        jacobiPrecond!(m_model, m_M, m_K)
        for axis = 0:2
            println("\n  Case ", axis + 1)
            # Compute the RHS: Boundary or Domain, rhsType: Boundary = 1 || Domain = 0, axis 0 = X || axis 1 = Y  || axis 2 = XY
            computeRHS!(m_model, m_RHS, axis, m_K, m_B)
            # Solver (to ensure optimal RAM usage we call GC before and after the PCGM):    
            GC.gc()
            #pcg_old!(m_model, m_X, m_RHS, m_M, m_K);
            pcg!(m_model, m_X, m_RHS, m_M, m_K)
            GC.gc()
            # Compute Effective Property:
            m_C .+= femEffective(m_model, m_X, axis, m_B)
            m_RHS .*= 0
            m_X .*= 0
        end
        m_M = nothing;
    elseif (m_model.solverType == 1) # Direct Method
        # Compute the RHS: Boundary or Domain, rhsType: Boundary = 1 || Domain = 0, axis 0 = X || axis 1 = Y
        m_RHS1 = zeros(Float64, m_model.nDOFs); m_RHS2 = zeros(Float64, m_model.nDOFs); m_RHS3 = zeros(Float64, m_model.nDOFs);
        computeRHS!(m_model, m_RHS1, 0, m_K, m_B)
        computeRHS!(m_model, m_RHS2, 1, m_K, m_B)
        computeRHS!(m_model, m_RHS3, 2, m_K, m_B)
        # Solver
        m_X1 = zeros(Float64, m_model.nDOFs); m_X2 = zeros(Float64, m_model.nDOFs); m_X3 = zeros(Float64, m_model.nDOFs);
        directMethod!(m_model, m_X1, m_X2, m_X3, m_RHS1, m_RHS2, m_RHS3, m_K)
        m_RHS1 = nothing; m_RHS2 = nothing; m_RHS3 = nothing;
        # Compute Effective Property:
        m_C = zeros(Float64, 3, 3)
        m_C .+= femEffective(m_model, m_X1, 0, m_B)
        m_C .+= femEffective(m_model, m_X2, 1, m_B)
        m_C .+= femEffective(m_model, m_X3, 2, m_B)
        m_X1 = nothing; m_X2 = nothing; m_X3 = nothing;
    end
    println("---------------------------")
    println("Effective Properties:\n")
    for i = 1:3
        println("C[$i,:] = ", m_C[i,:])
    end
    println("\n--------------------------------------")
    println("Reverting matrix\n")
    m_D = inv(m_C)
    for i = 1:3
        println("D[$i,:] = ", m_D[i,:])
    end
    println("\nE = ", 1 / m_D[1][1])
end