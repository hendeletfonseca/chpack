# Compute Effective PERMEABILITY Properties
using JSON
using SparseArrays
using LinearAlgebra

# Model data struct:
struct Model
    nx::UInt64;
    ny::UInt64;
    voxelSize::Float64;
    refinement::UInt64;
    rhsType::UInt8;
    solverType::UInt8;
    pcgTol::Float64;
    pcgIter::UInt64;
    nNodes::UInt64;
    nElems::UInt64;
    nDOFs::UInt64;
    nVelocityNodes::UInt64;
    nInterfaceNodes::UInt64;
    DOFMap::Vector{UInt64};
    FluidElems::Vector{UInt64};
    function Model(_nx, _ny, _voxelSize, _refinement, _rhsType, _solverType, _pcgTol, _pcgIter, _nNodes, _nElems, _nDOFs, _nVelocityNodes, _nInterfaceNodes, _DOFMap, _FluidElems)
        new(_nx, _ny, _voxelSize, _refinement, _rhsType, _solverType, _pcgTol, _pcgIter, _nNodes, _nElems, _nDOFs, _nVelocityNodes, _nInterfaceNodes, _DOFMap, _FluidElems);
    end
end

# Build Model: 
function buildModel(_JsonFile::String, _RawFile::String)
    println(".Building Model!")
    # Read Json file:
    nx, ny, voxelSize, refinement, rhsType, solverType, pcgTol, pcgIter, matKeys = readJSON(_JsonFile)
    # Read Raw file:
    elemMatMap = zeros(UInt16, nx * ny * refinement * refinement)
    readRAW!(nx, ny, refinement, matKeys, elemMatMap, _RawFile)
    # Update the parameters based on the given refinement level:
    nx *= refinement
    ny *= refinement
    nNodes::UInt64 = (nx + 1) * (ny + 1)
    nElems::UInt64 = (nx) * (ny)
    # Generate a map of Degree of Freedom:
    DOFMap = zeros(UInt64, nNodes)
    nVelocityNodes, nInterfaceNodes, nDOFs, FluidElems = generateDOFMap!(nx, ny, elemMatMap, DOFMap) 
    # Build the Model:
    model = Model(nx, ny, voxelSize, refinement, rhsType, solverType, pcgTol, pcgIter, nNodes, nElems, nDOFs, nVelocityNodes, nInterfaceNodes, DOFMap, FluidElems)
    matKeys = nothing
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
        for i = 1:nMat
            matKeys[convert(UInt8, materials[i][1]) + 1] = i
        end
        materials = nothing
        data = nothing
        return nx, ny, voxelSize, refinement, rhsType, solverType, pcgTol, pcgIter, matKeys
    end
end

# Read RAW file:
function readRAW!(_nx::UInt64, _ny::UInt64, _refinement::UInt64, _matKeys::Vector{UInt16}, _elemMatMap::Vector{UInt16}, _filename::String)
    println("   .Reading RAW!");
    # Initializations:
    nelem::UInt64 = _nx * _ny;  line::UInt64 = _ny * _refinement * _refinement - _ny;
    el::UInt64 = 0;  buffer::UInt64 = 0; c::UInt64 = 0; r::UInt64 = 0;
    # Open and read file:
    open(_filename, "r") do io
        bin_array = read(io)
        bin_array = reshape(bin_array, convert(Int64, _nx), convert(Int64, _ny))         
        # Build the element material map based on the refinement level:
        for e = 1:nelem
            c = ((e - 1) ÷ _ny) + 1
            r = e - (c - 1) * _ny
            buffer = _matKeys[bin_array[c,r] + 1]
            for i = 1:_refinement
                for j = 1:_refinement
                    el = e + ((e - 1) ÷ _ny) * line + (j - 1) + ((i - 1) % _refinement) * _ny * _refinement + ((e - 1) % _ny) * (_refinement - 1)
                    _elemMatMap[el] = buffer * (1 - (buffer >= 2)) + 1 * (buffer >= 2) + -1 * (buffer == 1)
                end        
            end 
        end   
        bin_array = nothing
    end
end

# Generate the Degree of Freedom Map:
function generateDOFMap!(_nx::UInt64, _ny::UInt64, _elemMatMap::Vector{UInt16}, _DOFMap::Vector{UInt64})
    println("   .Generating the Map of DOFs (Degrees of Freedom)!")
    # Initializations
    nNodes::UInt64 = (_nx + 1) * (_ny + 1); c::UInt64 = 0; r::UInt64 = 0;
    flag_lastcol::UInt64 = 0; flag_lastrow::UInt64 = 0; flag_firstcol::UInt64 = 0; flag_firstrow::UInt64 = 0;
    elem_se::UInt64 = 0; elem_ne::UInt64 = 0; elem_nw::UInt64 = 0; elem_sw::UInt64 = 0;
    numVelocityNodes::UInt64 = 0; numInterfaceNodes::UInt64 = 0; numFluidElems::UInt64 = 0;
    value_CaseVelocity::UInt64 = 0; value_CaseInterface::UInt64 = 0;
    top_node::UInt64 = 0; left_node::UInt64 = 0;
    vFluidElems = zeros(UInt64, nNodes)
    @fastmath @inbounds @simd for N = 1:nNodes
        # Nodes's column and row
        c = ((N - 1) ÷ (_ny + 1)) + 1
        r = (N - (c - 1) * (_ny + 1))
        # Flags - Boundary nodes
        flag_lastcol = (c / (_nx + 1) == 1)
        flag_lastrow = (r / (_ny + 1) == 1)
        flag_firstcol = (c == 1)
        flag_firstrow = (r == 1)
        # Check Neighboring   |nw|ne|
        #                      - N - 
        #                     |sw|se|
        elem_se = N - c + 1 + flag_lastrow * (-_ny) + flag_lastcol * (-_nx * (_ny + 1) + c - 1)
        elem_ne = N - c + flag_firstrow * (_ny) + flag_lastcol * (-_nx * _ny)
        elem_nw = N - c - _ny + flag_firstrow * (_ny) + flag_firstcol * (_nx) * _ny
        elem_sw = N - c + 1 - _ny + flag_firstcol * (_nx) * _ny + flag_lastrow * (-_ny)
        # Flags - type of node (At least one fluid elem / At least one solid elem / Only fluid elems)
        flag_oneFluidElem =  - (_elemMatMap[elem_se] * _elemMatMap[elem_ne] * _elemMatMap[elem_nw] * _elemMatMap[elem_sw] - 1)
        flag_oneSolidElem =  (_elemMatMap[elem_se] + _elemMatMap[elem_ne] + _elemMatMap[elem_nw] + _elemMatMap[elem_sw]) > 0
        flag_onlyFluidElems =  (_elemMatMap[elem_se] + _elemMatMap[elem_ne] + _elemMatMap[elem_nw] + _elemMatMap[elem_sw]) == 0
        # Creation of DofMap
        numVelocityNodes += flag_onlyFluidElems * (1 - flag_lastrow) * (1 - flag_lastcol)
        numInterfaceNodes += flag_oneSolidElem * flag_oneFluidElem * (1 - flag_lastrow) * (1 - flag_lastcol)  
        value_CaseVelocity = numVelocityNodes
        value_CaseInterface = nNodes + numInterfaceNodes      
        _DOFMap[N] = flag_onlyFluidElems * value_CaseVelocity + flag_oneFluidElem * flag_oneSolidElem * value_CaseInterface
        # Application of PBC 
        top_node = flag_lastrow * (N - _ny) + (1 - flag_lastrow)
        left_node = flag_lastcol * (N - _nx * (_ny + 1)) + (1 - flag_lastcol)  
        _DOFMap[N] += flag_lastrow * (-_DOFMap[N] + _DOFMap[top_node]) + flag_lastcol * (-_DOFMap[N] + _DOFMap[left_node]) + flag_lastrow * flag_lastcol * (+_DOFMap[N] - _DOFMap[left_node])
        # Vector of Fluid Elements
        numFluidElems += (_elemMatMap[elem_se] == 0) * (1 - flag_lastrow) * (1 - flag_lastcol);
        vFluidElems[numFluidElems * (_elemMatMap[elem_se] == 0) + (_elemMatMap[elem_se] != 0)] += ((_elemMatMap[elem_se] == 0) * elem_se) * (1 - flag_lastrow) * (1 - flag_lastcol)
    end
    # Correction of the numbering of interface nodes and reduction of the vector of fluid elements
    vFluidElems_new = zeros(UInt64, numFluidElems);
    for N = 1:nNodes
        _DOFMap[N] += (_DOFMap[N] > numVelocityNodes) * (-nNodes + numVelocityNodes)
        vFluidElems_new[(N <= numFluidElems) * N + (N > numFluidElems)] += vFluidElems[N] * (N <= numFluidElems)
    end
    # Total number of degrees of freedom
    numDOFs::UInt64 = 3 * numVelocityNodes + numInterfaceNodes
    vFluidElems = nothing
    return numVelocityNodes, numInterfaceNodes, numDOFs, vFluidElems_new
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
        mem = (16 * _model.nElems + 64 * _model.nNodes + 5 * 64 * _model.nDOFs) / 8 / 1_000_000
    elseif (_model.solverType == 1)
        mem = (16 * _model.nElems + 64 * _model.nNodes + 2 * 64 * _model.nDOFs + 18 * 64 * _model.nElems) / 8 / 1_000_000
    end
    println("   $(_model.nDOFs) DOFs")
    println("   $mem MB")
    println("---------------------------")
end

# Compute the element matrices for the fluid:
function finiteElementMatrices!(_K::Matrix{Float64}, _G::Matrix{Float64}, _Pe::Matrix{Float64}, _F::Vector{Float64}, _SN::Matrix{Float64})
    println("   .Computing FEM matrices!")
    # Compute the matrices
    delta::Float64 = 1.0
    coordsElem = [0. 0.; delta 0.; delta delta; 0. delta]
    rr = [-1.0 / sqrt(3) 1.0 / sqrt(3)]
    ww = [1.0 1.0]
    ss = rr
    C = [2.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 1.0]
    # h2 = ((dx)^2)+((dy)^2);
    h2::Float64 = 2.0
    stab::Float64 = h2 / 12 / 1
    for i = 1:2
        r = rr[1,i]
        for j = 1:2
            s = ss[1,j]
            N = (1.0 / 4.0) * [(1 - s) * (1 - r) 0 (1 - s) * (1 + r) 0 (1 + s) * (1 + r) 0 (1 - r) * (1 + s) 0;
                     0 (1 - s) * (1 - r) 0 (1 - s) * (1 + r) 0 (1 + s) * (1 + r) 0 (1 - r) * (1 + s)];  
            dN1dr = -1 / 4 * (1 - s)
            dN2dr = +1 / 4 * (1 - s)
            dN3dr = +1 / 4 * (1 + s)
            dN4dr = -1 / 4 * (1 + s)  
            dN1ds = -1 / 4 * (1 - r)
            dN2ds = -1 / 4 * (1 + r)
            dN3ds = +1 / 4 * (1 + r)
            dN4ds = +1 / 4 * (1 - r)
            DN = [dN1dr dN2dr dN3dr dN4dr;
                  dN1ds dN2ds dN3ds dN4ds];
            J = DN * coordsElem
            invJ = 1.0 / det(J) * [J[2,2] -J[1,2]; -J[2,1] J[1,1]]
            DNxy = invJ * DN
            B = [DNxy[1,1]         0 DNxy[1,2]         0 DNxy[1,3]         0 DNxy[1,4]         0;
                         0 DNxy[2,1]         0 DNxy[2,2]         0 DNxy[2,3]         0 DNxy[2,4];
                 DNxy[2,1] DNxy[1,1] DNxy[2,2] DNxy[1,2] DNxy[2,3] DNxy[1,3] DNxy[2,4] DNxy[1,4]];
            _K .= _K .+ B' * C * B * det(J) * ww[1,i] * ww[1,j]
            _G .= _G .+ B' * [1.0;1.0;0.0] * [N[1,1] N[1,3] N[1,5] N[1,7]] * (det(J)) * ww[1,i] * ww[1,j]
            Bp = invJ * DN
            _Pe .= _Pe .+ stab * Bp' * Bp * det(J) * ww[1,i] * ww[1,j]
            _F .= _F .+ [N[1,1]; N[1,3]; N[1,5]; N[1,7]] * det(J) * ww[1,i] * ww[1,j]
            _SN .= _SN .+ N * det(J) * ww[1,i] * ww[1,j]
        end
    end 
end

# Compute the RHS: Boundary or Domain, rhsType: Boundary = 0 || Domain = 1, axis 0 = X || axis 1 = Y
function computeRHS!(_model::Model, _RHS::Vector{Float64}, _axis::UInt64, _fe::Vector{Float64}, _k::Matrix{Float64}, _g::Matrix{Float64}, _p::Matrix{Float64})
    # Initializations
    numVelocityNodes::UInt64 = _model.nVelocityNodes
    ny::UInt64 = _model.ny
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0;
    nN1::UInt64 = 0; nN2::UInt64 = 0; nN3::UInt64 = 0; nN4::UInt64 = 0;
    # Compute each RHS (_axis) based on boundary or domain data (_rhsType)
    if _model.rhsType == 0  # Domain
        println("   .Computing RHS - Domain!")
        c::UInt64 = 0
        if _axis == 0
            for e in _model.FluidElems
                c = ((e - 1) ÷ ny) + 1
                N1 = e + c
                N3 = N1 + ny
                N2 = N3 + 1
                N4 = N1 - 1
                nN1 = _model.DOFMap[N1]
                nN2 = _model.DOFMap[N2]
                nN3 = _model.DOFMap[N3]
                nN4 = _model.DOFMap[N4]
                _RHS[(numVelocityNodes >= nN1) * (nN1 * 2 - 2) + 1] += _fe[1] * (numVelocityNodes >= nN1)
                _RHS[(numVelocityNodes >= nN2) * (nN2 * 2 - 2) + 1] += _fe[2] * (numVelocityNodes >= nN2)
                _RHS[(numVelocityNodes >= nN3) * (nN3 * 2 - 2) + 1] += _fe[3] * (numVelocityNodes >= nN3)
                _RHS[(numVelocityNodes >= nN4) * (nN4 * 2 - 2) + 1] += _fe[4] * (numVelocityNodes >= nN4)
            end
        elseif _axis == 1
            for e in _model.FluidElems
                c = ((e - 1) ÷ ny) + 1
                N1 = e + c
                N3 = N1 + ny
                N2 = N3 + 1
                N4 = N1 - 1
                nN1 = _model.DOFMap[N1]
                nN2 = _model.DOFMap[N2]
                nN3 = _model.DOFMap[N3]
                nN4 = _model.DOFMap[N4]
                _RHS[(numVelocityNodes >= nN1) * (nN1 * 2 - 1) + 1] += _fe[1] * (numVelocityNodes >= nN1)
                _RHS[(numVelocityNodes >= nN2) * (nN2 * 2 - 1) + 1] += _fe[2] * (numVelocityNodes >= nN2)
                _RHS[(numVelocityNodes >= nN3) * (nN3 * 2 - 1) + 1] += _fe[3] * (numVelocityNodes >= nN3)
                _RHS[(numVelocityNodes >= nN4) * (nN4 * 2 - 1) + 1] += _fe[4] * (numVelocityNodes >= nN4)
            end
        end
    elseif _model.rhsType == 1  # Boundary
        println("   .Computing RHS - Boundary!")      
        isOnlyFluid = zeros(UInt8, 8)
        pElemDOFNum = zeros(UInt64, 8)
        if _axis == 0     # _axis 0 = X
            deltaX::Float64 = _model.nx
            numVDOFs::UInt64 = 2 * numVelocityNodes
            for e in _model.FluidElems
                if (e >= (_model.nElems - _model.ny + 1))
                    N1 = e + ((e - 1) ÷ ny) + 1
                    N3 = N1 + ny
                    N2 = N3 + 1
                    N4 = N1 - 1
                    nN1 = _model.DOFMap[N1]
                    nN2 = _model.DOFMap[N2]
                    nN3 = _model.DOFMap[N3]
                    nN4 = _model.DOFMap[N4]
                    # Applying pressure on the left border
                    _RHS[(nN1 * 2 - 2) + 1] += -_g[1,2] * deltaX - _g[1,3] * deltaX
                    _RHS[(nN1 * 2 - 1) + 1] += -_g[2,2] * deltaX - _g[2,3] * deltaX
                    _RHS[(nN2 * 2 - 2) + 1] += -_g[3,2] * deltaX - _g[3,3] * deltaX
                    _RHS[(nN2 * 2 - 1) + 1] += -_g[4,2] * deltaX - _g[4,3] * deltaX
                    _RHS[(nN3 * 2 - 2) + 1] += -_g[5,2] * deltaX - _g[5,3] * deltaX
                    _RHS[(nN3 * 2 - 1) + 1] += -_g[6,2] * deltaX - _g[6,3] * deltaX
                    _RHS[(nN4 * 2 - 2) + 1] += -_g[7,2] * deltaX - _g[7,3] * deltaX
                    _RHS[(nN4 * 2 - 1) + 1] += -_g[8,2] * deltaX - _g[8,3] * deltaX
                    _RHS[(nN1 + numVDOFs - 1) + 1] += -_p[1,2] * deltaX - _p[1,3] * deltaX
                    _RHS[(nN2 + numVDOFs - 1) + 1] += -_p[2,2] * deltaX - _p[2,3] * deltaX
                    _RHS[(nN3 + numVDOFs - 1) + 1] += -_p[3,2] * deltaX - _p[3,3] * deltaX
                    _RHS[(nN4 + numVDOFs - 1) + 1] += -_p[4,2] * deltaX - _p[4,3] * deltaX
                end
            end
        end
    end
end

# Jacobi Preconditioner: assembly || M
function jacobiPrecond!(_model::Model, _M::Vector{Float64}, _k::Matrix{Float64}, _p::Matrix{Float64})
    println("   .Jacobi Preconditioner!")
    # Initializations:
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0;
    nN1::UInt64 = 0; nN2::UInt64 = 0; nN3::UInt64 = 0; nN4::UInt64 = 0;
    numVelocityNodes::UInt64 = _model.nVelocityNodes
    numVDOFs::UInt64 = 2 * numVelocityNodes
    c::UInt64 = 0;
    # Compute the preconditioner: 
    @fastmath @inbounds @simd for e in _model.FluidElems
        c = ((e - 1) ÷ _model.ny) + 1
        N1 = e + c
        N3 = N1 + _model.ny
        N2 = N3 + 1
        N4 = N1 - 1
        nN1 = _model.DOFMap[N1]
        nN2 = _model.DOFMap[N2]
        nN3 = _model.DOFMap[N3]
        nN4 = _model.DOFMap[N4]
        _M[(numVelocityNodes >= nN1) * (nN1 * 2 - 2) + 1] += _k[1,1] * (numVelocityNodes >= nN1)
        _M[(numVelocityNodes >= nN2) * (nN2 * 2 - 2) + 1] += _k[3,3] * (numVelocityNodes >= nN2)
        _M[(numVelocityNodes >= nN3) * (nN3 * 2 - 2) + 1] += _k[5,5] * (numVelocityNodes >= nN3)
        _M[(numVelocityNodes >= nN4) * (nN4 * 2 - 2) + 1] += _k[7,7] * (numVelocityNodes >= nN4)
        _M[(numVelocityNodes >= nN1) * (nN1 * 2 - 1) + 1] += _k[2,2] * (numVelocityNodes >= nN1)
        _M[(numVelocityNodes >= nN2) * (nN2 * 2 - 1) + 1] += _k[4,4] * (numVelocityNodes >= nN2)
        _M[(numVelocityNodes >= nN3) * (nN3 * 2 - 1) + 1] += _k[6,6] * (numVelocityNodes >= nN3)
        _M[(numVelocityNodes >= nN4) * (nN4 * 2 - 1) + 1] += _k[8,8] * (numVelocityNodes >= nN4)
        _M[numVDOFs + nN1] -= _p[1,1]
        _M[numVDOFs + nN2] -= _p[2,2]
        _M[numVDOFs + nN3] -= _p[3,3]
        _M[numVDOFs + nN4] -= _p[4,4]
    end
    _M .= _M .\ 1
end 

# Preconditioned Conjugate Gradient Method:
function pcg_old!(_model::Model, _x::Vector{Float64}, _r::Vector{Float64}, _M::Vector{Float64}, _k::Matrix{Float64}, _g::Matrix{Float64}, _p::Matrix{Float64})
    println("   .PCG Solver!")
    # Initializations
    numVelocityNodes = _model.nVelocityNodes;
    numVDOFs = 2 * numVelocityNodes;
    d = zeros(Float64, _model.nDOFs);
    q = zeros(Float64, _model.nDOFs);
    N1  = 0; N2  = 0; N3  = 0; N4  = 0; 
    nN1 = 0; nN2 = 0; nN3 = 0; nN4 = 0;
    c = 0; q_temp = 0;
    pElemDOFNumP = zeros(UInt64, 4);
    pElemDOFNumV = zeros(UInt64, 8);
    isOnlyFluid = zeros(UInt8, 8);
    # PCG Initialization:
    d .= _r;
    d .*= _M;
    delta_new = (_r' * d)[1,1];
    delta_0 = delta_new;
    i_max = _model.pcgIter;
    ii = 0;
    # PCG Iterations:
    while (ii < i_max) && (abs(delta_new) > _model.pcgTol * _model.pcgTol * abs(delta_0))
    # while (ii<=i_max) && (maximum(abs.(_r))>_pcgTol)
        @fastmath @inbounds @simd for e in _model.FluidElems
            c = ((e - 1) ÷ _model.ny) + 1;
            N1 = e + c;
            N3 = N1 + _model.ny;
            N2 = N3 + 1;
            N4 = N1 - 1;
            nN1 = _model.DOFMap[N1];
            nN2 = _model.DOFMap[N2];
            nN3 = _model.DOFMap[N3];
            nN4 = _model.DOFMap[N4];
            pElemDOFNumP[1] = numVDOFs + nN1;
            pElemDOFNumP[2] = numVDOFs + nN2;
            pElemDOFNumP[3] = numVDOFs + nN3;
            pElemDOFNumP[4] = numVDOFs + nN4;
            pElemDOFNumV[1] = ((numVelocityNodes >= nN1) * (nN1 * 2 - 2)) + 1;
            pElemDOFNumV[2] = ((numVelocityNodes >= nN1) * (nN1 * 2 - 1)) + 1;
            pElemDOFNumV[3] = ((numVelocityNodes >= nN2) * (nN2 * 2 - 2)) + 1;
            pElemDOFNumV[4] = ((numVelocityNodes >= nN2) * (nN2 * 2 - 1)) + 1;
            pElemDOFNumV[5] = ((numVelocityNodes >= nN3) * (nN3 * 2 - 2)) + 1;
            pElemDOFNumV[6] = ((numVelocityNodes >= nN3) * (nN3 * 2 - 1)) + 1;
            pElemDOFNumV[7] = ((numVelocityNodes >= nN4) * (nN4 * 2 - 2)) + 1;
            pElemDOFNumV[8] = ((numVelocityNodes >= nN4) * (nN4 * 2 - 1)) + 1;
            isOnlyFluid[1] = (numVelocityNodes >= nN1);
            isOnlyFluid[2] = (numVelocityNodes >= nN1);
            isOnlyFluid[3] = (numVelocityNodes >= nN2);
            isOnlyFluid[4] = (numVelocityNodes >= nN2);
            isOnlyFluid[5] = (numVelocityNodes >= nN3);
            isOnlyFluid[6] = (numVelocityNodes >= nN3);
            isOnlyFluid[7] = (numVelocityNodes >= nN4);
            isOnlyFluid[8] = (numVelocityNodes >= nN4);
            for i = 1:8
                q_temp = 0;
                for j = 1:8
                    q_temp += _k[i,j] * d[pElemDOFNumV[j]] * isOnlyFluid[j];
                end
                q[pElemDOFNumV[i]] += q_temp * isOnlyFluid[i];
                q_temp = 0;
                for j = 1:4
                    q_temp += _g[i,j] * d[pElemDOFNumP[j]];
                end
                q[pElemDOFNumV[i]] += q_temp * isOnlyFluid[i];
            end
            for i = 1:4
                q_temp = 0;
                for j = 1:4
                    q_temp -= _p[i,j] * d[pElemDOFNumP[j]];
                end
                q[pElemDOFNumP[i]] += q_temp;
                q_temp = 0;
                for j = 1:8
                    q_temp += _g[j,i] * d[pElemDOFNumV[j]] * isOnlyFluid[j];
                end
                q[pElemDOFNumP[i]] += q_temp;
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
    println("    $ii steps")
    println("    Residue = ", sqrt(abs(delta_new) / abs(delta_0)))
end

# Preconditioned Conjugate Gradient Method:
function pcg!(_model::Model, _x::Vector{Float64}, _r::Vector{Float64}, _M::Vector{Float64}, _k::Matrix{Float64}, _g::Matrix{Float64}, _p::Matrix{Float64})
    println("   .PCG Solver!")
    # Initializations:    
    d = zeros(Float64, _model.nDOFs)
    q = zeros(Float64, _model.nDOFs)
    pElemDOFNumP = zeros(UInt64, 4)
    pElemDOFNumV = zeros(UInt64, 8)
    pElemDOFNumV_local = zeros(UInt8, 8)
    pElemDOFVarP = zeros(Float64, 4)
    pElemDOFVarV = zeros(Float64, 8)
    isNodeFluid::Bool = false
    dofs = zeros(UInt64, 4)
    count_elemDOFV::UInt64 = 0
    q_temp::Float64 = 0.0
    alfa::Float64 = 0.0
    beta::Float64 = 0.0
    ri::Float64   = 0.0
    qi::Float64   = 0.0
    delta_new::Float64 = 0.0
    delta_old::Float64 = 0.0
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0;
    numVelocityNodes::UInt64 = _model.nVelocityNodes
    numVDOFs::UInt64 = 2 * numVelocityNodes
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
        @fastmath @inbounds @simd for e in _model.FluidElems
            N1 = e + ((e - 1) ÷ _model.ny) + 1; N2 = N1 + _model.ny + 1; N3 = N2 - 1; N4 = N1 - 1;
            dofs[1] = _model.DOFMap[N1]; dofs[2] = _model.DOFMap[N2]; dofs[3] = _model.DOFMap[N3]; dofs[4] = _model.DOFMap[N4];
            count_elemDOFV = 0            
            @inbounds for i = 1:4
                pElemDOFNumP[i] = numVDOFs + dofs[i]
                pElemDOFVarP[i] = d[pElemDOFNumP[i]]
                isNodeFluid = (numVelocityNodes >= dofs[i])
                pElemDOFNumV[count_elemDOFV+1] = 2*dofs[i]-1; pElemDOFNumV[count_elemDOFV+2] = 2*dofs[i];
                pElemDOFNumV_local[count_elemDOFV+1] = 2*i-1; pElemDOFNumV_local[count_elemDOFV+2] = 2*i;
                count_elemDOFV += 2*isNodeFluid
                pElemDOFVarV[2*i-1] = isNodeFluid ? d[2*dofs[i]-1] : 0.0
                pElemDOFVarV[2*i]   = isNodeFluid ?   d[2*dofs[i]] : 0.0
            end
            @inbounds for i = 1:count_elemDOFV
                iii::UInt16 = pElemDOFNumV_local[i]
                q_temp = _k[iii,1] * pElemDOFVarV[1]
                @inbounds for j=2:8; q_temp += _k[iii,j] * pElemDOFVarV[j]; end
                @inbounds for j=1:4; q_temp += _g[iii,j] * pElemDOFVarP[j]; end
                q[pElemDOFNumV[i]] += q_temp
            end
            @inbounds for i = 1:4
                q_temp = _g[1,i] * pElemDOFVarV[1]
                @inbounds for j=2:8; q_temp += _g[j,i] * pElemDOFVarV[j]; end
                @inbounds for j=1:4; q_temp -= _p[i,j] * pElemDOFVarP[j]; end
                q[pElemDOFNumP[i]] += q_temp
            end
        end
        alfa = delta_new / dot(d,q)
        delta_old = delta_new
        delta_new= 0.0
        @inbounds for i=1:_model.nDOFs
            _x[i] += d[i]*alfa
            _r[i] -= q[i]*alfa
            ri = _r[i]
            qi = ri*_M[i]
            q[i] = qi
            delta_new += ri*qi;       
        end
        if (abs(delta_new) <= tolerance); iteration_count = ii; break; end
        beta = delta_new / delta_old
        @inbounds for i=1:_model.nDOFs
            d[i] *= beta
            d[i] += q[i]
            q[i]  = 0.0
        end
    end
    println("    $iteration_count steps")
    println("    Residue = ", sqrt(abs(delta_new) / abs(delta_0)))
end

# Direct Solver: [K] 64 bits * m_nDOFs * m_nDOFs 
function directMethod!(_model::Model, _x1::Vector{Float64}, _x2::Vector{Float64}, _RHS1::Vector{Float64}, _RHS2::Vector{Float64}, _K::Matrix{Float64}, _G::Matrix{Float64}, _Pe::Matrix{Float64})
    println("   .Direct Solver!")
    # Initializations:
    numVelocityNodes::UInt64 = _model.nVelocityNodes
    numVDOFs::UInt64 = 2 * numVelocityNodes
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0;
    nN1::UInt64 = 0; nN2::UInt64 = 0; nN3::UInt64 = 0; nN4::UInt64 = 0;
    pElemDOFNum = zeros(UInt64, 12)
    isOnlyFluid = zeros(UInt8, 8)
    A = spzeros(_model.nDOFs, _model.nDOFs)
    # Assembly system matrix:
    @fastmath @inbounds @simd for e in _model.FluidElems
        N1 = e + 1 + div(e - 1, _model.ny);
        N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
        nN1 = _model.DOFMap[N1]
        nN2 = _model.DOFMap[N2]
        nN3 = _model.DOFMap[N3]
        nN4 = _model.DOFMap[N4]
        pElemDOFNum[1] = ((numVelocityNodes >= nN1) * (nN1 * 2 - 2)) + 1
        pElemDOFNum[2] = ((numVelocityNodes >= nN1) * (nN1 * 2 - 1)) + 1
        pElemDOFNum[3] = ((numVelocityNodes >= nN2) * (nN2 * 2 - 2)) + 1
        pElemDOFNum[4] = ((numVelocityNodes >= nN2) * (nN2 * 2 - 1)) + 1
        pElemDOFNum[5] = ((numVelocityNodes >= nN3) * (nN3 * 2 - 2)) + 1
        pElemDOFNum[6] = ((numVelocityNodes >= nN3) * (nN3 * 2 - 1)) + 1
        pElemDOFNum[7] = ((numVelocityNodes >= nN4) * (nN4 * 2 - 2)) + 1
        pElemDOFNum[8] = ((numVelocityNodes >= nN4) * (nN4 * 2 - 1)) + 1
        pElemDOFNum[9] = numVDOFs + nN1
        pElemDOFNum[10] = numVDOFs + nN2
        pElemDOFNum[11] = numVDOFs + nN3
        pElemDOFNum[12] = numVDOFs + nN4
        isOnlyFluid[1] = (numVelocityNodes >= nN1)
        isOnlyFluid[2] = (numVelocityNodes >= nN1)
        isOnlyFluid[3] = (numVelocityNodes >= nN2)
        isOnlyFluid[4] = (numVelocityNodes >= nN2)
        isOnlyFluid[5] = (numVelocityNodes >= nN3)
        isOnlyFluid[6] = (numVelocityNodes >= nN3)
        isOnlyFluid[7] = (numVelocityNodes >= nN4)
        isOnlyFluid[8] = (numVelocityNodes >= nN4)
        for i = 1:8
            for j = 1:8
                A[pElemDOFNum[i],pElemDOFNum[j]] += _K[i,j] * isOnlyFluid[i] * isOnlyFluid[j]
            end
        end
        for i = 1:4
            for j = 1:4
                A[pElemDOFNum[i + 8],pElemDOFNum[j + 8]] += -_Pe[i,j]
            end
        end
        for i = 1:8
            for j = 1:4
                A[pElemDOFNum[i],pElemDOFNum[j + 8]] += _G[i,j] * isOnlyFluid[i]
            end
        end
        for i = 1:4
            for j = 1:8
                A[pElemDOFNum[i + 8],pElemDOFNum[j]] += _G[j,i] * isOnlyFluid[j]
            end
        end
    end
    # Solve for three rhs:
    _x1 .= A \ _RHS1
    _x2 .= A \ _RHS2
end

# Compute Velocity-FEM Effective property:
function femEffective(_model::Model, _X::Vector{Float64}, _axis::UInt64, _SN::Matrix{Float64})
    println("   .Updating Constitutive Matrix!")
    # Initializations
    nx::UInt64 = _model.nx
    ny::UInt64 = _model.ny
    numVelocityNodes::UInt64 = _model.nVelocityNodes
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0;
    nN1::UInt64 = 0; nN2::UInt64 = 0; nN3::UInt64 = 0; nN4::UInt64 = 0;
    c::UInt64 = 0; q_temp = 0;
    # loop over the fluid elements
    QXx::Float64 = 0; QXy::Float64  = 0; QYx::Float64 = 0; QYy::Float64 = 0 ; V::Float64 = 0;
    for e in _model.FluidElems
        c = ((e - 1) ÷ ny) + 1
        N1 = e + c
        N3 = N1 + ny
        N2 = N3 + 1
        N4 = N1 - 1
        nN1 = _model.DOFMap[N1]
        nN2 = _model.DOFMap[N2]
        nN3 = _model.DOFMap[N3]
        nN4 = _model.DOFMap[N4]
        if _axis == 0
            QXx += _SN[1,1] * _X[(numVelocityNodes >= nN1) * (nN1 * 2 - 2) + 1] * (numVelocityNodes >= nN1)
            QXx += _SN[1,3] * _X[(numVelocityNodes >= nN2) * (nN2 * 2 - 2) + 1] * (numVelocityNodes >= nN2)
            QXx += _SN[1,5] * _X[(numVelocityNodes >= nN3) * (nN3 * 2 - 2) + 1] * (numVelocityNodes >= nN3)
            QXx += _SN[1,7] * _X[(numVelocityNodes >= nN4) * (nN4 * 2 - 2) + 1] * (numVelocityNodes >= nN4)
            QYx += _SN[2,2] * _X[(numVelocityNodes >= nN1) * (nN1 * 2 - 1) + 1] * (numVelocityNodes >= nN1)
            QYx += _SN[2,4] * _X[(numVelocityNodes >= nN2) * (nN2 * 2 - 1) + 1] * (numVelocityNodes >= nN2)
            QYx += _SN[2,6] * _X[(numVelocityNodes >= nN3) * (nN3 * 2 - 1) + 1] * (numVelocityNodes >= nN3)
            QYx += _SN[2,8] * _X[(numVelocityNodes >= nN4) * (nN4 * 2 - 1) + 1] * (numVelocityNodes >= nN4)
        elseif _axis == 1
            QXy += _SN[1,1] * _X[(numVelocityNodes >= nN1) * (nN1 * 2 - 2) + 1] * (numVelocityNodes >= nN1)
            QXy += _SN[1,3] * _X[(numVelocityNodes >= nN2) * (nN2 * 2 - 2) + 1] * (numVelocityNodes >= nN2)
            QXy += _SN[1,5] * _X[(numVelocityNodes >= nN3) * (nN3 * 2 - 2) + 1] * (numVelocityNodes >= nN3)
            QXy += _SN[1,7] * _X[(numVelocityNodes >= nN4) * (nN4 * 2 - 2) + 1] * (numVelocityNodes >= nN4)
            QYy += _SN[2,2] * _X[(numVelocityNodes >= nN1) * (nN1 * 2 - 1) + 1] * (numVelocityNodes >= nN1)
            QYy += _SN[2,4] * _X[(numVelocityNodes >= nN2) * (nN2 * 2 - 1) + 1] * (numVelocityNodes >= nN2)
            QYy += _SN[2,6] * _X[(numVelocityNodes >= nN3) * (nN3 * 2 - 1) + 1] * (numVelocityNodes >= nN3)
            QYy += _SN[2,8] * _X[(numVelocityNodes >= nN4) * (nN4 * 2 - 1) + 1] * (numVelocityNodes >= nN4)
        end    
    end
    V = nx * ny 
    C = [ QXx / V QXy / V; QYx / V QYy / V]   
    return C * _model.voxelSize * _model.voxelSize / (_model.refinement  * _model.refinement )
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
    # Compute the conductivity matrix for each Material:
    m_K = zeros(Float64, 8, 8);       m_G = zeros(Float64, 8, 4);       m_Pe = zeros(Float64, 4, 4);
    m_F = zeros(Float64, 4);          m_SN = zeros(Float64, 2, 8);      m_C = zeros(Float64, 2, 2);
    finiteElementMatrices!(m_K, m_G, m_Pe, m_F, m_SN)
    if (m_model.solverType == 0) # Preconditioned Conjugate Gradient Method
        # Initialize the effective tensor, the right hand side, the inicial guess and the preconditioner:
        m_M = zeros(Float64, m_model.nDOFs) 
        m_X = zeros(Float64, m_model.nDOFs)  
        m_RHS = zeros(Float64, m_model.nDOFs)
        # Compute the Jacobi preconditioner:
        jacobiPrecond!(m_model, m_M, m_K, m_Pe)
        #m_M = ones(Float64, m_model.nDOFs);
        for axis = 0:1
            println("\n  Case ", axis + 1)
            # Compute the RHS: Boundary or Domain, rhsType: Boundary = 1 || Domain = 0, axis 0 = X || axis 1 = Y 
            computeRHS!(m_model, m_RHS, axis, m_F, m_K, m_G, m_Pe)  
            # Solver (to ensure optimal RAM usage we call GC before and after the PCGM):    
            GC.gc()
            #pcg_old!(m_model, m_X, m_RHS, m_M, m_K, m_G, m_Pe);
            pcg!(m_model, m_X, m_RHS, m_M, m_K, m_G, m_Pe)
            GC.gc()     
            # Compute Effective Property:
            m_C .+= femEffective(m_model, m_X, axis, m_SN)
            m_RHS .*= 0
            m_X .*= 0
        end
        M = nothing
    elseif (m_model.solverType == 1) # Direct Method
        # Compute the RHS: Boundary or Domain, rhsType: Boundary = 1 || Domain = 0, axis 0 = X || axis 1 = Y
        m_RHS1 = zeros(Float64, m_model.nDOFs); m_RHS2 = zeros(Float64, m_model.nDOFs);
        computeRHS!(m_model, m_RHS1, 0, m_F, m_K, m_G, m_Pe) 
        computeRHS!(m_model, m_RHS2, 1, m_F, m_K, m_G, m_Pe)
        # Solver
        m_X1 = zeros(Float64, m_model.nDOFs); m_X2 = zeros(Float64, m_model.nDOFs);
        directMethod!(m_model, m_X1, m_X2, m_RHS1, m_RHS2, m_K, m_G, m_Pe)
        m_RHS1 = nothing; m_RHS2 = nothing;
        # Compute Effective Property:       
        m_C .+= femEffective(m_model, m_X1, 0, m_SN)
        m_C .+= femEffective(m_model, m_X2, 1, m_SN)
    end
    println("---------------------------")
    println("Effective Properties:\n")
    for i = 1:2
        println("C[$i,:] = ", m_C[i,:])
    end
    println("\n--------------------------------------")
    return m_C
end
