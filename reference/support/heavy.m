function heavy_prefix_indices(n)
    value := GetEnv("HEAVY_PREFIX_INDICES");
    if value eq "" then return [1..n]; end if;
    indices := [ StringToInteger(s) : s in Split(value,",") ];
    error if IsEmpty(indices) or #Seqset(indices) ne #indices or
        &or[ i lt 1 or i gt n : i in indices ],
        "MISMATCH: HEAVY_PREFIX_INDICES must contain distinct valid indices";
    return indices;
end function;

function heavy_processes()
    value := GetEnv("REFERENCE_PROCESSES");
    n := value eq "" select 8 else StringToInteger(value);
    error if n lt 1, "MISMATCH: REFERENCE_PROCESSES must be positive";
    return n;
end function;

function heavy_target_limit(n)
    value := GetEnv("HEAVY_TARGET_LIMIT");
    if value eq "" then return n; end if;
    m := StringToInteger(value);
    error if m lt 1, "MISMATCH: HEAVY_TARGET_LIMIT must be positive";
    return Min(m,n);
end function;

function heavy_shard(processes)
    value := GetEnv("REFERENCE_SHARD");
    n := StringToInteger(value);
    error if n lt 1 or n gt processes, "MISMATCH: invalid reference shard";
    return n;
end function;
