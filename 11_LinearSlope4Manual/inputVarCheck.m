%#codegen
function ErrCode = inputVarCheck(fvar,ErrCodeValue)
         ErrCode=0;
        if isempty(fvar) 
            ErrCode  = ErrCodeValue;
        end
        for ii=1:length(fvar)
            if  isnan(fvar(ii))
                ErrCode  = ErrCodeValue;
            end
        end
end