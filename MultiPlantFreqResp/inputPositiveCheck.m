%#codegen
%%previous is zero;
%ver02 @ Dec-06-2016 ,chg nErrCode from 1 to 0,if no err;
function nErrCode = inputPositiveCheck(fvar,ErrCodeValue)
         nErrCode=0;%%0 value means no err, keep the same as others @Dec-06-2016
       
        for ii=1:length(fvar)
            if fvar<=0
                nErrCode = ErrCodeValue;
            end
        end
end