# Compute Effective CONDUCTIVITY Properties
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
    matKeys::Array{UInt16};
    matProp::Array{Float64}; 
    nNodes::UInt64;
    nElems::UInt64;
    nDOFs::UInt64;
    DOFMap::Array{UInt64}; 
    elemMatMap::Array{UInt64}; 
    function Model(_nx, _ny, _voxelSize, _refinement, _nMat, _rhsType, _solverType, _pcgTol, _pcgIter, _matKeys, _matProp, _nNodes, _nElems, _nDOFs, _DOFMap, _elemMatMap)
        new(_nx, _ny, _voxelSize, _refinement, _nMat, _rhsType, _solverType, _pcgTol, _pcgIter, _matKeys, _matProp, _nNodes, _nElems, _nDOFs, _DOFMap, _elemMatMap);
    end
end

# Build Model:
function buildModel(_JsonFile::String, _RawFile::String)
    println(".Building Model!")
    # Read Json file:
    nx, ny, voxelSize, refinement, nMat, rhsType, solverType, pcgTol, pcgIter, matKeys, matProp = readJSON(_JsonFile);
    # Read Raw file:
    elemMatMap = zeros(UInt64, nx * ny * refinement * refinement);
    readRAW!(nx, ny, refinement, matKeys, elemMatMap, _RawFile);
    # Update the parameters based on the given refinement level:
    nx *= refinement;
    ny *= refinement;
    nNodes = (nx + 1) * (ny + 1);
    nElems = (nx) * (ny);
    DOFperNode = 1;
    nDOFs = nElems * DOFperNode;
    # Generate a map of Degree of Freedom:
    DOFMap = zeros(UInt64, nNodes);
    generateDOFMap!(nx, ny, DOFMap);    
    # Build the Model:
    model = Model(nx, ny, voxelSize, refinement, nMat, rhsType, solverType, pcgTol, pcgIter, matKeys, matProp, nNodes, nElems, nDOFs, DOFMap, elemMatMap);
    println("---------------------------");
    return model
end

# Read JSON file:
function readJSON(_filename::String)
    println("   .Reading JSON!")
    # Open and read file:
    open(_filename, "r") do f
        data = JSON.parse(f);
        nx = data["image_dimensions"][1];
        ny = data["image_dimensions"][2];
        refinement = 1;
        if haskey(data, "refinement"); refinement = data["refinement"]; end
        voxelSize = 1;
        if haskey(data, "voxel_size"); voxelSize = data["voxel_size"]; end
        rhsType = 0;
        if haskey(data, "type_of_rhs"); rhsType = data["type_of_rhs"]; end
        solverType = 0;
        if haskey(data, "type_of_solver"); solverType = data["type_of_solver"]; end
        pcgTol = 0.000001;
        if haskey(data, "solver_tolerance"); pcgTol = data["solver_tolerance"]; end
        pcgIter = nx * ny * refinement * refinement;
        if haskey(data, "number_of_iterations"); pcgIter = data["number_of_iterations"]; end
        nMat = data["number_of_materials"];
        materials = data["properties_of_materials"];
        matKeys = zeros(UInt16, 256);
        matProp = zeros(Float64, 256);
        for i = 1:nMat
            matKeys[convert(UInt8, materials[i][1]) + 1] = i;
            matProp[convert(UInt8, materials[i][1]) + 1] = convert(Float64, materials[i][2]);
        end
        materials = nothing;
        data = nothing;
        return nx, ny, voxelSize, refinement, nMat, rhsType, solverType, pcgTol, pcgIter, matKeys, matProp
    end
end

# Read RAW file:
function readRAW!(_nx::Int, _ny::Int, _refinement::Int, _matKeys::Array{UInt16,1}, _elemMatMap::Array{UInt64,1}, _filename::String)
    println("   .Reading RAW!");
    # Initializations:
    nelem = _nx * _ny; el = 0; line = _ny * _refinement * _refinement - _ny; buffer = 0;
    # Open and read file:
    open(_filename, "r") do io
        bin_array = read(io);        
        # Build the element material map based on the refinement level:
        for e = 1:nelem
            buffer = _matKeys[bin_array[e] + 1];
            for i = 1:_refinement
                for j = 1:_refinement
                    el = e + ((e - 1) ÷ _ny) * line + (j - 1) + ((i - 1) % _refinement) * _ny * _refinement + ((e - 1) % _ny) * (_refinement - 1);
                    _elemMatMap[el] = buffer;
                end        
            end 
        end        
        bin_array = nothing;
    end
end

# Generate the Degree of Freedom Map:
function generateDOFMap!(_nx::Int, _ny::Int, _DOFMap::Array{UInt64,1})
    println("   .Generating the Map of DOFs (Degrees of Freedom)!");
    # Number the DOFs following the nodes from top to bottom and left to right:
    nElems = _nx * _ny;
    nNodes = (_nx + 1) * (_ny + 1);
    for n = 1:nNodes
        i = (n - 1) % nNodes;
        _DOFMap[n] = (i - (i ÷ (_ny + 1)) - _ny * ((i % (_ny + 1) ÷ _ny))) % nElems + 1;
    end
end

# Estimate memory consuption:
function estimateMemory(_model::Model)
    println(".Estimating memory!");
    # elemMatMap = 16 bits * nElems
    # DOFMap = 64 bits * nNodes
    # RHS = 64 bits * nDOFs
    # PCGM Solver   / solverType == 0 / M d x q =  * 64 bits * nDOFs
    # Direct Solver / solverType == 1 / K = 18 * 64 bits * nElems (rough sparse estimative)
    mem = 0;
    if (_model.solverType == 0)
        mem += (16 * _model.nElems + 64 * _model.nNodes + 5 * 64 * _model.nDOFs) / 8 / 1_000_000;
    elseif (_model.solverType == 1)
        mem += (16 * _model.nElems + 64 * _model.nNodes + 2 * 64 * _model.nDOFs + 18 * 64 * _model.nElems) / 8 / 1_000_000;
    end
    println("   $(_model.nDOFs) DOFs");
    println("   $mem MB");
    println("---------------------------");
end

# Compute the element conductivity matrix for each material:
function elementConductivityMatrices!(_model::Model, _K::Array{Float64,3}, _B::Array{Float64,3})
    println("   .Computing each element conductivity matrix!");
    # Compute the matrices for each material:
    i = 0;
    for elemProps in _model.matProp[_model.matKeys .!= 0]
        i += 1;
        _K[:,:,i], _B[:,:,i] = Q4ElementConductivity(elemProps);
    end
end

# Element Q4 Conductivity - FEM:
function Q4ElementConductivity(_elemProps::Float64)
    # Initializations:
    K = zeros(Float64, 4, 4);
    B = zeros(Float64, 2, 4);
    # Analytical k and B for a Pixel-Based FEM: 
    a = 2 / 3 * _elemProps; b = 1 / 3 * _elemProps; c = 1 / 6 * _elemProps;
    q = 1 / 2 * _elemProps;
    K[1,1] = +a; K[1,2] = -c; K[1,3] = -b; K[1,4] = -c;
    K[2,1] = -c; K[2,2] = +a; K[2,3] = -c; K[2,4] = -b;
    K[3,1] = -b; K[3,2] = -c; K[3,3] = +a; K[3,4] = -c;
    K[4,1] = -c; K[4,2] = -b; K[4,3] = -c; K[4,4] = +a;
    B[1,1] = -q; B[1,2] = +q; B[1,3] = +q; B[1,4] = -q;
    B[2,1] = -q; B[2,2] = -q; B[2,3] = +q; B[2,4] = +q;
    return K, B
end

# Compute the RHS: Boundary or Domain, rhsType: Boundary = 0 || Domain = 1, axis 0 = X || axis 1 = Y
function computeRHS!(_model::Model, _RHS::Array{Float64,1}, _axis::Int, _K::Array{Float64,3}, _B::Array{Float64,3})
    println("   .Computing RHS!");
    # Initializations:
    N1 = 0; N2 = 0; N3 = 0; N4 = 0; e = 0; 
    # Compute each RHS (_axis) based on boundary or domain data (_model.m_rhsType):
    if _model.rhsType == 1  # Boundary         
        if _axis == 0
            deltaT = _model.nx;
            c = _model.nx;
            for r = 1:_model.ny
                e = r + (c - 1) * _model.ny;
                N1 = e + c; N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                _RHS[_model.DOFMap[N1]] -= (_K[1,2,_model.elemMatMap[e]] + _K[1,3,_model.elemMatMap[e]]) * deltaT;
                _RHS[_model.DOFMap[N2]] -= (_K[2,2,_model.elemMatMap[e]] + _K[2,3,_model.elemMatMap[e]]) * deltaT;
                _RHS[_model.DOFMap[N3]] -= (_K[3,2,_model.elemMatMap[e]] + _K[3,3,_model.elemMatMap[e]]) * deltaT;
                _RHS[_model.DOFMap[N4]] -= (_K[4,2,_model.elemMatMap[e]] + _K[4,3,_model.elemMatMap[e]]) * deltaT;
            end
        elseif _axis == 1
            deltaT = _model.ny;
            r = 1;
            for c = 1:_model.nx
                e = r + (c - 1) * _model.ny;
                N1 = e + c; N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                _RHS[_model.DOFMap[N1]] -= (_K[1,3,_model.elemMatMap[e]] + _K[1,4,_model.elemMatMap[e]]) * deltaT;
                _RHS[_model.DOFMap[N2]] -= (_K[2,3,_model.elemMatMap[e]] + _K[2,4,_model.elemMatMap[e]]) * deltaT;
                _RHS[_model.DOFMap[N3]] -= (_K[3,3,_model.elemMatMap[e]] + _K[3,4,_model.elemMatMap[e]]) * deltaT;
                _RHS[_model.DOFMap[N4]] -= (_K[4,3,_model.elemMatMap[e]] + _K[4,4,_model.elemMatMap[e]]) * deltaT;
            end
        end
    elseif _model.rhsType == 0  # Domain      
        for e = 1:_model.nElems
            N1 = e + ((e - 1) ÷ _model.ny) + 1; N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            _RHS[_model.DOFMap[N1]] += _B[_axis + 1,1,_model.elemMatMap[e]];
            _RHS[_model.DOFMap[N2]] += _B[_axis + 1,2,_model.elemMatMap[e]];
            _RHS[_model.DOFMap[N3]] += _B[_axis + 1,3,_model.elemMatMap[e]];
            _RHS[_model.DOFMap[N4]] += _B[_axis + 1,4,_model.elemMatMap[e]];
        end    
    end
end

# Direct Solver: [K] 64 bits * m_nDOFs * m_nDOFs 
function directMethod!(_model::Model, _x1::Array{Float64,1}, _x2::Array{Float64,1}, _RHS1::Array{Float64,1}, _RHS2::Array{Float64,1}, _K::Array{Float64,3})
    println("   .Direct Solver!");
    # Initializations:
    K = spzeros(_model.nDOFs, _model.nDOFs);
    pElemDOFNum = zeros(UInt64, 4);
    N1 = 0; N2 = 0; N3 = 0; N4 = 0;
    # Assembly system matrix:
    for e = 1:_model.nElems
        N1 = e + ((e - 1) ÷ _model.ny) + 1; N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
        pElemDOFNum[1] = _model.DOFMap[N1]; pElemDOFNum[2] = _model.DOFMap[N2]; pElemDOFNum[3] = _model.DOFMap[N3]; pElemDOFNum[4] = _model.DOFMap[N4];
        for i = 1:4
            for j = 1:4
                K[pElemDOFNum[i],pElemDOFNum[j]] += _K[i,j,_model.elemMatMap[e]];
            end
        end
    end
    # Solve for two rhs:
    _x1 .= K \ _RHS1;
    _x2 .= K \ _RHS2;
end

# Jacobi Preconditioner: assembly || M
function jacobiPrecond!(_model::Model, _M::Array{Float64,1}, _K::Array{Float64,3})
    println("   .Jacobi Preconditioner!");
    # Initializations:
    N1 = 0; N2 = 0; N3 = 0; N4 = 0;
    # Compute the preconditioner: 
    for e = 1:_model.nElems
        N1 = e + ((e - 1) ÷ _model.ny) + 1; N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
        _M[_model.DOFMap[N1]] += _K[1,1,_model.elemMatMap[e]];
        _M[_model.DOFMap[N2]] += _K[2,2,_model.elemMatMap[e]];
        _M[_model.DOFMap[N3]] += _K[3,3,_model.elemMatMap[e]];
        _M[_model.DOFMap[N4]] += _K[4,4,_model.elemMatMap[e]];
    end
    _M .= _M .\ 1;
end 

# Preconditioned Conjugate Gradient Method:
function pcg!(_model::Model, _x::Array{Float64,1}, _r::Array{Float64,1}, _M::Array{Float64,1}, _K::Array{Float64,3})
    println("   .PCG Solver!");
    # Initializations:    
    d = zeros(Float64, _model.nDOFs);
    q = zeros(Float64, _model.nDOFs);
    pElemDOFNum = zeros(UInt64, 4);
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
            pElemDOFNum[1] = _model.DOFMap[N1]; pElemDOFNum[2] = _model.DOFMap[N2]; pElemDOFNum[3] = _model.DOFMap[N3]; pElemDOFNum[4] = _model.DOFMap[N4];
            for i = 1:4
                q_temp = 0;
                for j = 1:4
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

# Compute Flux-FEM Effective property:
function femEffective(_model::Model, _T::Array{Float64,1}, _axis::Int, _B::Array{Float64,3})
    println("   .Updating Constitutive Matrix!");
    # Initializations:
    QX = 0; QY = 0;
    N1 = 0; N2 = 0; N3 = 0; N4 = 0;
    pElemDOFNum = zeros(UInt64, 4);
    C = zeros(Float64, 2, 2);
    qx = 0; qy = 0;

    # opening txt file:
    if (_axis == 0)
        t_txt = open("pyview/values/t0.txt", "w");
        qx_txt = open("pyview/values/qx0.txt", "w");
        qy_txt = open("pyview/values/qy0.txt", "w");
    elseif (_axis == 1)
        t_txt = open("pyview/values/t1.txt", "w");
        qx_txt = open("pyview/values/qx1.txt", "w");
        qy_txt = open("pyview/values/qy1.txt", "w");
    end

    # Compute the effective properties for each test: 
    if _model.rhsType == 1  # Boundary
        if _axis == 0
            deltaT = _model.nx;   
            for eb = _model.nElems - (_model.ny - 1):_model.nElems  
                QX += (_B[1,2,_model.elemMatMap[eb]] * deltaT); QX += (_B[1,3,_model.elemMatMap[eb]] * deltaT);
                QY += (_B[2,2,_model.elemMatMap[eb]] * deltaT); QY += (_B[2,3,_model.elemMatMap[eb]] * deltaT);             
            end    
        elseif _axis == 1
            deltaT = _model.ny;
            for eb = 1:(_model.ny):_model.nElems
                QX += (_B[1,3,_model.elemMatMap[eb]] * deltaT); QX += (_B[1,4,_model.elemMatMap[eb]] * deltaT);
                QY += (_B[2,3,_model.elemMatMap[eb]] * deltaT); QY += (_B[2,4,_model.elemMatMap[eb]] * deltaT);             
            end 
        end
        for e = 1:_model.nElems
            qx = 0; qy = 0;
            N1 = e + ((e - 1) ÷ _model.ny) + 1; N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            pElemDOFNum[1] = N1; pElemDOFNum[2] = N2; pElemDOFNum[3] = N3; pElemDOFNum[4] = N4;
            for i = 1:4
                QX += (_B[1,i,_model.elemMatMap[e]] * _T[_model.DOFMap[pElemDOFNum[i]]]);
                QY += (_B[2,i,_model.elemMatMap[e]] * _T[_model.DOFMap[pElemDOFNum[i]]]);
                qx += (_B[1,i,_model.elemMatMap[e]] * _T[_model.DOFMap[pElemDOFNum[i]]]);
                qy += (_B[2,i,_model.elemMatMap[e]] * _T[_model.DOFMap[pElemDOFNum[i]]]);
            end            
            write(qx_txt, string(qx, "\n"));
            write(qy_txt, string(qy, "\n"));
            #write(qx_txt, string((_B[1,4,_model.elemMatMap[e]] * _T[_model.DOFMap[pElemDOFNum[4]]]), '\n'));
            #write(qy_txt, string((_B[2,4,_model.elemMatMap[e]] * _T[_model.DOFMap[pElemDOFNum[4]]]), '\n'));
            if _axis == 0 write(t_txt, string(_T[_model.DOFMap[pElemDOFNum[4]]], "\n"));
            elseif _axis == 1 write(t_txt, string(_T[_model.DOFMap[pElemDOFNum[2]]], "\n"));
            end

        end        
    elseif _model.rhsType == 0  # Domain
        t = zeros(Float64, 4);
        if (_axis == 0);     t[1] = 0; t[2] = 1; t[3] = 1; t[4] = 0;
        elseif (_axis == 1); t[1] = 0; t[2] = 0; t[3] = 1; t[4] = 1; end

        for e = 1:_model.nElems            
            N1 = e + ((e - 1) ÷ _model.ny) + 1; N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            pElemDOFNum[1] = N1; pElemDOFNum[2] = N2; pElemDOFNum[3] = N3; pElemDOFNum[4] = N4;
            for i = 1:4
                QX += (_B[1,i,_model.elemMatMap[e]] * (t[i] - _T[_model.DOFMap[pElemDOFNum[i]]]));
                QY += (_B[2,i,_model.elemMatMap[e]] * (t[i] - _T[_model.DOFMap[pElemDOFNum[i]]]));
            end
            write(t_txt, string(_T[_model.DOFMap[pElemDOFNum[4]]], "\n"));
            write(qx_txt, string((_B[1,4,_model.elemMatMap[e]] * (t[4] - _T[_model.DOFMap[pElemDOFNum[4]]])), '\n'));
            write(qy_txt, string((_B[2,4,_model.elemMatMap[e]] * (t[4] - _T[_model.DOFMap[pElemDOFNum[4]]])), '\n'));
        end
    end
    close(t_txt);
    close(qx_txt);
    close(qy_txt);
    C[1,_axis + 1] = QX / _model.nElems; C[2,_axis + 1] = QY / _model.nElems;
    return C
end

# -----------------
function main(_arg)
    println("---------------------------");
    # Build the Model data struct:
    m_model = buildModel(_arg * ".json", _arg * ".raw");
    # Estimate Memory Consumption:
    estimateMemory(m_model);
    # SOLVE:
    println(".Solving");
    # Compute the conductivity matrix for each Material:
    m_K = zeros(Float64, 4, 4, m_model.nMat);
    m_B = zeros(Float64, 2, 4, m_model.nMat);
    elementConductivityMatrices!(m_model, m_K, m_B);
    if (m_model.solverType == 0) # Preconditioned Conjugate Gradient Method
        # Initialize the effective tensor, the right hand side, the inicial guess and the preconditioner:
        m_C = zeros(Float64, 2, 2);
        m_RHS = zeros(Float64, m_model.nDOFs);   
        m_X = zeros(Float64, m_model.nDOFs);
        m_M = zeros(Float64, m_model.nDOFs);
        # Compute the Jacobi preconditioner:
        jacobiPrecond!(m_model, m_M, m_K);
        for axis = 0:1
            println("\n  Case ", axis + 1);
            # Compute the RHS: Boundary or Domain, rhsType: Boundary = 1 || Domain = 0, axis 0 = X || axis 1 = Y
            computeRHS!(m_model, m_RHS, axis, m_K, m_B);            
            # Solver (to ensure optimal RAM usage we call GC before and after the PCGM):    
            GC.gc();
            pcg!(m_model, m_X, m_RHS, m_M, m_K);
            GC.gc();        
            # Compute Effective Property:
            m_C .+= femEffective(m_model, m_X, axis, m_B);
            m_RHS .*= 0;
            m_X .*= 0;
        end
        m_M = nothing;
    elseif (m_model.solverType == 1) # Direct Method
        # Compute the RHS: Boundary or Domain, rhsType: Boundary = 1 || Domain = 0, axis 0 = X || axis 1 = Y
        m_RHS1 = zeros(Float64, m_model.nDOFs); m_RHS2 = zeros(Float64, m_model.nDOFs);
        computeRHS!(m_model, m_RHS1, 0, m_K, m_B);
        computeRHS!(m_model, m_RHS2, 1, m_K, m_B);
        # Solver
        m_X1 = zeros(Float64, m_model.nDOFs); m_X2 = zeros(Float64, m_model.nDOFs);
        directMethod!(m_model, m_X1, m_X2, m_RHS1, m_RHS2, m_K);
        m_RHS1 = nothing; m_RHS2 = nothing;
        # Compute Effective Property:
        m_C = zeros(Float64, 2, 2);
        m_C .+= femEffective(m_model, m_X1, 0, m_B);
        m_C .+= femEffective(m_model, m_X2, 1, m_B);
        m_X1 = nothing; m_X2 = nothing;
    end
    println("---------------------------");
    println("Effective Properties:\n")
    for i = 1:2
        println("C[$i,:] = ", m_C[i,:]);
    end
    println("--------------------------------------");
end

# Starts application
if length(ARGS) > 0
    @time main(ARGS[1])
end
