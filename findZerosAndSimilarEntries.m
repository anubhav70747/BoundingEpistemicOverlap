function [indicesList, zeroList,matrix] = findZerosAndSimilarEntries(matrix, tolerance)
% FINDZEROSANDSIMILARENTRIES_CPLX_UT
% Same behaviour as the complex-aware routine, but it only scans entries
% (r,c) with  c ≥ r  (upper triangular, diagonal included).

    [rows, cols] = size(matrix);

    indicesList = {};
    zeroList    = {};
    nZeros      = 0;

    for r = 1:rows
        r
        % ---- iterate only over the upper-triangular columns ----------
        for c = r : cols          % <-- key change
            z = matrix(r,c);

            % ---------- 1. “Zero” test (rectangular region) ------------
            if abs(real(z)) <= tolerance && abs(imag(z)) <= tolerance
                nZeros          = nZeros + 1;
                zeroList{nZeros}= [r, c];
                continue
            end

            % ---------- 2. Check existing groups -----------------------
            found = false;
            for k = 1:numel(indicesList)
                idx        = indicesList{k};
                zExisting  = matrix(sub2ind([rows,cols], idx(:,1), idx(:,2)));

                sameReal = abs(real(zExisting) - real(z)) <= tolerance;
                sameImag = abs(abs(imag(zExisting)) - abs(imag(z))) <= tolerance;

                if any(sameReal & sameImag)
                    indicesList{k} = [idx; r, c];
                    found          = true;
                    break
                end
            end

            % ---------- 3. Create a new group ---------------------------
            if ~found
                indicesList{end+1} = [r, c];
            end
        end
    end
end